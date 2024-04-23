!
!   SCID-TDSE: Simple 1-electron atomic TDSE solver
!   Copyright (C) 2015-2024 Serguei Patchkovskii, Serguei.Patchkovskii@mbi-berlin.de
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!   Common data for the spherical_tdse program driver
!
module spherical_tdse_data
  use accuracy
  use bicg_tools
  use cap_tools
  use checkpoint_tools
  use composition_analysis
  use constants
  use math
  use node_tools
  use potential_tools
  use propagator_tools
  use rotation_tools
  use spherical_data
  use spherical_data_initialize
  use spherical_tsurf
  use spherical_tsurf_data
  use test_tools
  use timer
  use tridiagonal_tools
  use vectorpotential_tools
  use wavefunction_tools
  implicit none
  public
  public rcsid_spherical_tdse_data
  !
  character(len=clen), save :: rcsid_spherical_tdse_data = "$Id: spherical_tdse_data.f90,v 1.1 2024/02/13 14:22:14 ps Exp $"
  !
  integer(ik)              :: verbose               = 2_ik         ! How verbose do we need to be?
  integer                  :: omp_num_threads       = 0            ! Non-zero value will cause number of OpenMP threads
                                                                   ! to be set explicitly. It looks like some libraries
                                                                   ! mess with the environment variables, and stop
                                                                   ! OMP_NUM_THREADS from working.
  character(len=clen)      :: comment               = ' '          ! Descriptive string, to be copied to the output
  integer(ik)              :: ensemble_size         = 1_ik         ! Number of wavefunctions propagated in lockstep.
  integer(ik)              :: ensemble_index        = 1_ik         ! Index of the current wavefunction in the ensemble
                                                                   ! Only relevant during the input, and not used otherwise
  character(len=clen)      :: initial_wfn           = 'atomic'     ! Choice of the initial wavefunction. One of:
                                                                   ! 'random' = Use random initial guess
                                                                   ! 'unit'   = Use unit vector as the initial guess
                                                                   ! 'atomic' = Use atomic field-free solution (full matrix)
                                                                   !            Must supply initial_wfn_index.
                                                                   ! 'single' = Use atomic filed-free solution (single solution)
                                                                   !            Must supply initial_wfn_index and initial_wfn_energy
                                                                   ! 'read'   = Read initial wfn from a file
  integer(ik)              :: initial_wfn_index(3)  = (/0,0,1/)    ! Initial wavefunction quantum numbers; L,M,I
  complex(rk)              :: initial_wfn_energy    = -1._rk       ! Initial guess for wavefunction energy
  character(len=clen)      :: initial_wfn_file      = ' '          ! File containing initial wavefunction
  character(len=clen)      :: task                  = 'real time'  ! What to do after the initialization. One of:
                                                                   ! 'real time'      = Propagate in real time
                                                                   ! 'imaginary time' = Propagate in imaginary time
  logical                  :: skip_tests            = .true.       ! Skip all non-essential tests, and assume user knows
                                                                   ! where she is going
  logical                  :: do_dipole(3)          = (/ .true., .true., .true. /)
                                                                   ! do_dipole(1) = .True. if the expectation of dipole moment should be calculated
                                                                   ! do_dipole(2) = .True. if the dipole velocity should be calculated
                                                                   ! do_dipole(3) = .True. if the dipole acceleration should be calculated
  logical                  :: do_dipole_plasma      = .true.       ! .True. will include the "plasma" term in the dipole velocity and acceleration.
                                                                   ! Effectively, we simply add the free-electron contribution for
                                                                   ! the part of electron density which already left the box.
                                                                   ! WARNING: The resulting dipole acceleration will differ from
                                                                   ! WARNING: the numerical derivative of the dipole expectation
                                                                   ! WARNING: or dipole velocity
  logical                  :: skip_left_propagation = .false.      ! Do not propagate the left wavefunction. This makes the calculation almost
                                                                   ! twice as fast, but causes incorrect results for almost all observables.
                                                                   ! The only observables guaranteed to work correctly are the final right
                                                                   ! wavefunction, and the photoelectron spectrum for sts_atend_mode='direct'.
                                                                   ! This parameter exists almost exclusively for benchmark bragging rights.
  logical                  :: fake_left_propagation = .true.       ! For evaluation of properties during time propagation, replace left
                                                                   ! wavefunction by the complex conjugate of the right wavefunction.
                                                                   ! This is nearly always an approximation to the exact right wavefunction,
                                                                   ! but it may produce reasonable observables when the Hamiltonian matrix
                                                                   ! is nearly Hermitian.
                                                                   ! fake_left_propagation has no effect unless skip_left_propagation==.true.
  logical                  :: reconstruct_left      = .false.      ! If set to .True., the left wavefunction will be reconstructed from the
                                                                   ! right at the end of time propagation.
  real(xk)                 :: dt                    = 0.01_xk      ! Time step, in atomic units of time
  integer(ik)              :: timesteps             = 50000_ik     ! Number of time steps to perform
  character(len=clen)      :: dt_subdivision        = 'on'         ! Time step subdivision conrol. Can be one of:
                                                                   ! 'off'     - Do not subdivide time steps
                                                                   ! 'on'      - Subdivide the time steps. Same as 'lmax'.
                                                                   ! 'lmax'    - Keep (DT_EFF*(LMAX+1)<DT_MAX_LA).
                                                                   ! 'lmax-a'  - Keep (DT_EFF*(LMAX+1)*APOT<DT_MAX_LA).
                                                                   ! 'lmax-a2' - Keep (DT_EFF*(LMAX+1)*APOT**2<DT_MAX_LA).
  real(xk)                 :: dt_max_la             = 2.00_xk      ! Criterion for time-step subdivision.
                                                                   ! Smaller thresholds are generally needed for grids dense close to
                                                                   ! the origin. This value has no effect unless dt_subdivision=='on'
  integer(ik)              :: dt_interpolant_width  = 5_ik         ! Width of the interpolant used to evaluate vector-potential at the subdivided
                                                                   ! timestep points.
  character(len=clen)      :: rotation_mode         = 'auto'       ! Choice of the frame rotation implementation. One of:
                                                                   ! 'none'        - Use fixed frame. Vector-potential must be along Z;
                                                                   !                 other vector-potentials will be accepted, but 
                                                                   !                 (silently) produce incorrect results.
                                                                   !                 Choosing 'none' will also force field_unwrap=.false.
                                                                   ! 'brute force' - Explicit finite-angle rotation, using full-matrix
                                                                   !                 implementation through Wigner rotation matrices.
                                                                   !                 This version becomes numerically unstable for large
                                                                   !                 angular momenta, and is intended primarily for 
                                                                   !                 debugging.
                                                                   ! 'sparse'      - Small-angle rotation using sparse matrices. Will
                                                                   !                 implement large-angle rotation by breaking it up
                                                                   !                 in small steps (see rt_max_rot below, and the code)
                                                                   ! 'auto'        - Switches between 'none' and 'sparse' depending 
                                                                   !                 on the vector potential.
  logical                  :: field_unwrap          = .true.       ! Try to unwrap vector-potential in spherical coordinates
  real(rk)                 :: unwrap_threshold      = 1e-8_rk      ! Do not try to manipulate vector potential magnitude smaller 
                                                                   ! than unwrap_threshold * max(vp)
  character(len=clen)      :: field_preview         = 'field.table'! File containing laser field during the simulation.
  integer(ik)              :: output_each           = 100_ik       ! Reduce summary output by a factor output_each
  integer(ik)              :: detail_frequency      = 1_ik         ! Frequency of output on detail_output
  character(len=clen)      :: detail_output         = 'detail.table'
                                                                   ! File containing full-accuracy results from the simulation
  character(len=clen)      :: ensemble_output       = 'ensemble.table'
                                                                   ! Cross-w.f. output from an ensemble calculation. Only
                                                                   ! relevant if ensemble_size>1
  real(rk)                 :: composition_threshold = 1e-10_rk     ! Threshold for reporting the field-free amplitudes
                                                                   ! Set to a negative number to disable composition analysis
  real(rk)                 :: composition_max_energy= huge(1._rk)  ! Largest real part of eigenvalue to include in the analysis
  character(len=clen)      :: initial_wf_dump_prefix= ' '          ! Initial wavefunction dump in human-readable from
  character(len=clen)      :: final_wf_dump_prefix  = ' '          ! Final radial wavefunction in human-readable form is dumped in
                                                                   ! files with final_wf_dump_prefix used as a prefix. Empty
                                                                   ! string suppresses the dump
  character(len=clen)      :: visualize_prefix      = ' '          ! Prefix of the filename for wavefunction snapshots during the
                                                                   ! simulation. The actual filename will consist of the prefix,
                                                                   ! followed by a 10-digit integer (the timestep index) and ".data" 
                                                                   ! Blank disables visualization.
  integer(ik)              :: visualize_each        = 1000_ik      ! Produce a snapshot each (visualize_each) timesteps. Additionally,
                                                                   ! the initial and final wavefunctions are saved.
  integer(ik)              :: visualize_1stseq      = 1_ik         ! Index number of the first snapshot. Per-ensemble W.F. values
                                                                   ! are maintained in ens_visualize_1stseq(:); this value is used
                                                                   ! purely for input.
  logical                  :: vp_as_is              = .false.      ! If true, do not reset the vector-potential to zero at the beginning
                                                                   ! and the end of the simulation. Please do not use unless you know
                                                                   ! EXACTLY what you are doing.
  real(xk), allocatable    :: vpot_table(:,:)                      ! Vector-potential parameters table for the entire simulation
                                                                   ! First index: (0) = time, (1) = signed magnitude, (2) = theta, (3) = phi
                                                                   ! Second index: timestep, from -1 to 2*timesteps+1; even numbers are whole
                                                                   ! time steps; odd numbers are half-steps. 
                                                                   ! These data must be handled in at least double precision, or bad things
                                                                   ! will happen in time propagation
  real(xk), allocatable    :: efield_table(:,:)                    ! -(d A/d t) in the laboratory frame, otherwise known as the electric
                                                                   ! field. This quantity is needed for computing dipole acceleration.
                                                                   ! First index: 1=X, 2=Y, 3=Z components of the electric field
                                                                   ! Second index: timestep, from 0 to 2*timesteps; note that this is NOT
                                                                   ! the same as in vpot_table.
                                                                   ! efield_table is derived from the contents of vpot_table in
                                                                   ! fill_efield_table() below.
  integer(ik)              :: max_subdivision       = 1_ik         ! Largest timestep subdivision factor seen
  integer(ik)              :: timesteps_subdivided  = 0_ik         ! Number of time steps where subdivision was necessary
  integer(ik)              :: timesteps_micro       = 0_ik         ! Number of time steps after the subdivision
                                                                   ! The fields below are per-wavefunction in an ensemble.
                                                                   ! The index runs from 1 to ensemble_size.
  type(sd_wfn), allocatable    :: wfns_r(:), wfns_l(:)             ! Our wavefunctions (right and left).
  type(sts_data), allocatable  :: tsurfs(:)                        ! Running data for photoelectron spectrum calculation with t-SURF
  type(ckpt_data), allocatable :: ckpts(:)                         ! Global checkpointer state, one per wavefunction
  !
  !  Parameters which are per-ensemble. These arrays are allocated and initialized
  !  in spherical_tdse%prepare_initial_wavefunction. Their meaning is the same as
  !  that of the scalar parametrs with no ens_ prefix.
  !
  character(len=clen), allocatable :: &
      ens_detail_output(:), ens_final_wf_dump_prefix(:), ens_initial_wf_dump_prefix(:), &
      ens_visualize_prefix(:), ens_ckpt_save_basename(:), ens_ckpt_load_filename(:), &
      ens_volkov_opendx(:), ens_volkov_table(:), ens_coulomb_waves(:), &
      ens_coulomb_table(:), ens_coulomb_opendx(:)
  integer(ik), allocatable         :: iu_detail(:)                 ! Units for detailed output; remains open during the entire run
                                                                   ! Negative values indicate that the unit wasn't opened.
  integer(ik), allocatable         :: ens_visualize_1stseq(:)      ! Sequence number for visualization snapshot
  integer(ik)                      :: iu_ens = -1_ik               ! Unit for ensemble output; remains open during the entire run
  !
  !  Simulation parameters; we are collecting variables from many modules.
  !
  namelist /sph_tdse/ &
                      ! Parameters defined locally
                      verbose, comment, &
                      omp_num_threads, &
                      ensemble_size, ensemble_index, &
                      initial_wfn, initial_wfn_index, initial_wfn_energy, initial_wfn_file, &
                      task, dt, dt_subdivision, dt_max_la, dt_interpolant_width, timesteps, rotation_mode, &
                      field_unwrap, unwrap_threshold, field_preview, skip_tests, &
                      output_each, detail_frequency, detail_output, ensemble_output, &
                      composition_threshold, composition_max_energy, final_wf_dump_prefix, &
                      initial_wf_dump_prefix, &
                      vp_as_is, &
                      do_dipole, do_dipole_plasma, &
                      visualize_prefix, visualize_each, visualize_1stseq, &
                      skip_left_propagation, fake_left_propagation, reconstruct_left, &
                      ! Parameters from spherical_data
                      sd_nradial, sd_nspin, sd_lmax, sd_mmin, sd_mmax, &
                      sd_adaptive_l, sd_tolerance_l, sd_adaptive_r, sd_tolerance_r, &
                      sd_nradial_scl, sd_radial_edge, &
                      sd_rgrid, sd_rgrid_zeta, sd_rgrid_rmin, sd_rgrid_dr, sd_rgrid_r0, sd_rgrid_scale, sd_rgrid_npad, &
                      sd_rgrid_file, sd_rgrid_report, sd_rgrid_grad_delta, &
                      sd_operators, &
                      ! Parameters from potential_tools
                      pot_name, pot_param, pot_mask, pot_mask_r0, pot_mask_rx, pot_shift, &
                      pot_lmax, pot_tong05, pot_input, &
                      ! Parameters from propagator_tools
                      pt_mix_solver, pt_chunk, pt_sense_real, pt_eps_dt, pt_force_par_l, &
                      ! Parameters from vectorpotential_tools
                      vp_scale, vp_scale_x, vp_scale_x2, vp_shape, vp_param, vp_param_x, &
                      vp_param_x2, vp_table, &
                      ! Parameters from timer
                      timer_disable, &
                      ! Parameters from bicg_tools
                      bicg_failtrace, bicg_epsilon, bicg_maxiter, &
                      ! Parameters from rotation_tools
                      rt_max_rot, rt_blocksize, rt_sense_real, rt_nstep_max, &
                      ! Parameters from wavefunction_tools
                      wt_atomic_cache_prefix, wt_iterative_improvement, &
                      wt_disable_orthogonalization, wt_max_solution_iterations, &
                      wt_enable_memory_caches, wt_tme, wt_tme_file, &
                      ! Parameters from cap_tools
                      cap_name, cap_param, &
                      ! Parameters from spherical_tsurf
                      sts_volkov, sts_volkov_atend, sts_coulomb_atend, &
                      sts_step_frequency, sts_verbose, &
                      sts_rmatch, sts_kgrid, sts_dgrid, &
                      sts_kgrid_min, sts_kgrid_max, sts_kgrid_count, sts_kgrid_file, &
                      sts_dgrid_ntheta, sts_dgrid_nphi, sts_dgrid_count, sts_dgrid_file, &
                      sts_asymptotic_q, &
                      sts_volkov_table, sts_volkov_opendx, &
                      sts_coulomb_waves, sts_coulomb_table, sts_coulomb_opendx, &
                      sts_atend_block, sts_r2r_scale, sts_atend_mode, &
                      ! Parameters from tridiagonal_tools
                      m3d_solver, m3d_iterations, &
                      ! Parameters from checkpoint_tools
                      ckpt_save_basename, ckpt_load_filename, ckpt_max_checkpoints, &
                      ckpt_interval, &
                      ! Parameters from composition_analysis
                      ca_maxram, &
                      ! Parameters from node_tools
                      nt_node_output, nt_rebalance_interval, nt_use_multinode, nt_verbose, &
                      nt_max_requests
  !
  namelist /sph_multi/ &
                      ensemble_index, &
                      initial_wfn, initial_wfn_index, initial_wfn_energy, initial_wfn_file, &
                      detail_output, final_wf_dump_prefix, initial_wf_dump_prefix, &
                      visualize_1stseq, &
                      visualize_prefix, ckpt_save_basename, ckpt_load_filename, &
                      sts_volkov_opendx, sts_volkov_table, sts_coulomb_waves, &
                      sts_coulomb_table, sts_coulomb_opendx
  !
end module spherical_tdse_data
