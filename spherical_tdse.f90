!
!   SCID-TDSE: Simple 1-electron atomic TDSE solver
!   Copyright (C) 2015-2021 Serguei Patchkovskii, Serguei.Patchkovskii@mbi-berlin.de
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
!
!  Solving atomic TDSE in the presence of a strong laser field.
!  The structure of this code is strongly inspired by the paper:
!
!   H.G. Muller, An Efficient Propagation Scheme for the Time-Dependent 
!   Schroedinger Equation in the Velocity Gauge, Laser Physics, 9, 138-148 (1999)
!
!  However, the capabilities of the code, and some of the numerical 
!  details are quite different from HGM's. They are described in:
!
!   S. Patchkovskii and H.G. Muller, "Simple, accurate, and efficient 
!   implementation of 1-electron atomic time-dependent Schrödinger equation 
!   in spherical coordinates", Comp. Phys. Comm. 199, 153-169 (2016). 
!   doi:10.1016/j.cpc.2015.10.014
!
!  In parallel runs, this code may benefit significantly from large page support,
!  especially on systems with 2K default pages and small TLB (x86_64!). Please
!  refer to your system documentation on how to enable large page support. Most
!  recent Linux systems will come with anonymous huge page support enabled.
!
module spherical_tdse
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
  private
  public start
  public rcsid_spherical_tdse
  !
  character(len=clen), save :: rcsid_spherical_tdse = "$Id: spherical_tdse.f90,v 1.137 2023/12/23 11:09:20 ps Exp $"
  !
  integer, parameter       :: iu_detail             = 29           ! Unit for detailed output; remains open during the entire run
  integer, parameter       :: iu_temp               = 22           ! An arbitrary unit number, which can be used here
  !                                                 
  integer(ik)              :: verbose               = 2_ik         ! How verbose do we need to be?
  integer(ik)              :: omp_num_threads       = 0_ik         ! Non-zero value will cause number of OpenMP threads
                                                                   ! to be set explicitly. It looks like some libraries
                                                                   ! mess with the environment variables, and stop
                                                                   ! OMP_NUM_THREADS from working.
  character(len=clen)      :: comment               = ' '          ! Descriptive string, to be copied to the output
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
  character(len=clen)      :: detail_output         = 'detail.table'! File containing full-accuracy results from the simulation
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
  integer(ik)              :: visualize_1stseq      = 1_ik         ! Index number of the first snapshot
  logical                  :: vp_as_is              = .false.      ! If true, do not reset the vector-potential to zero at the beginning
                                                                   ! and the end of the simulation. Please do not use unless you know
                                                                   ! EXACTLY what you are doing.
  type(sd_wfn)             :: wfn_r                                ! Our wavefunction (right)
  type(sd_wfn)             :: wfn_l                                ! Ditto (left)
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
  type(sts_data)           :: tsurf                                ! Running data for photoelectron spectrum calculation with t-SURF
  type(ckpt_data)          :: ckpt                                 ! Global checkpointer state. 
  integer(ik)              :: max_subdivision       = 1_ik         ! Largest timestep subdivision factor seen
  integer(ik)              :: timesteps_subdivided  = 0_ik         ! Number of time steps where subdivision was necessary
  integer(ik)              :: timesteps_micro       = 0_ik         ! Number of time steps after the subdivision
  !
  !  Simulation parameters; we are collecting variables from many modules.
  !
  namelist /sph_tdse/ &
                      ! Parameters defined locally
                      verbose, comment, &
                      omp_num_threads, &
                      initial_wfn, initial_wfn_index, initial_wfn_energy, initial_wfn_file, &
                      task, dt, dt_subdivision, dt_max_la, dt_interpolant_width, timesteps, rotation_mode, &
                      field_unwrap, unwrap_threshold, field_preview, skip_tests, &
                      output_each, detail_frequency, detail_output, &
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
  contains
  !
  subroutine random_initial_wfn(tag,wfn)
    character(len=*), intent(in) :: tag
    type(sd_wfn), intent(inout)  :: wfn
    integer(ik)                  :: lval, mval, sval, ipt
    real(rk)                     :: hr(2)
    !
    call random_seed()
    loop_m: do mval=sd_mmin,sd_mmax
      loop_l: do lval=abs(mval),sd_lmax
        loop_s: do sval=1,sd_nspin
          loop_pt: do ipt=1,sd_nradial
            call random_number(hr)
            wfn%wfn(:,sval,lval,mval) = cmplx(hr(1),hr(2),kind=rk)
          end do loop_pt
        end do loop_s
      end do loop_l
    end do loop_m
    write (out,"('Initial ',a,' wavefunction set to random values')") trim(tag)
  end subroutine random_initial_wfn
  !
  subroutine atomic_initial_wfn
    integer(ik)              :: ipt, lval, mval, ind, nvec, alloc
    complex(rk), allocatable :: evec(:,:,:)   ! Eigenvectors
    complex(rk), allocatable :: eval(:)       ! Eigenvalues
    !
    lval = initial_wfn_index(1)
    mval = initial_wfn_index(2)
    ind  = initial_wfn_index(3)
    nvec = sd_nradial
    if (initial_wfn=='single') then
      ind  = 1
      nvec = 1
    end if
    !
    write (out,"('Choosing atomic solution with L=',i0,' M=',i0,' I=',i0)") lval, mval, ind
    if (lval<0 .or. lval>sd_lmax .or. mval<sd_mmin .or. mval>sd_mmax .or. &
        abs(mval)>lval .or. ind<1 .or. ind>sd_nradial) then
      stop 'spherical_tdse%atomic_initial_wfn - bad initial_wfn_index'
    end if
    !
    allocate (evec(sd_nradial,nvec,2),eval(nvec),stat=alloc)
    if (alloc/=0) stop 'spherical_tdse%atomic_initial_wfn - allocation failed'
    select case (initial_wfn)
      case default
        stop 'spherical_tdse%atomic_initial_wfn - bad initial_wfn'
      case ('atomic')
        call wt_atomic_solutions(verbose,lval,eval,evec)
      case ('single')
        eval(ind) = initial_wfn_energy
        call wt_one_atomic_solution(verbose,lval,eval(ind),evec(:,ind,:))
    end select
    !
    write (out,"('Energy of the atomic solution is ',2(g28.16e3,1x),'Hartree')") eval(ind)
    wfn_l%wfn = 0
    wfn_l%wfn(:,1,lval,mval) = evec(:,ind,1)
    wfn_r%wfn = 0
    wfn_r%wfn(:,1,lval,mval) = evec(:,ind,2)
    !
    if (verbose>=2) then 
      write (out,"(/'Initial radial wavefunction was:'/)")
      write (out,"((1x,a6,5(1x,a24)))") &
             ' I ', '  R,Bohr ', '  Re(psi_l)  ', '  Im(psi_l)  ', '  Re(psi_r)  ', '  Im(psi_r)  ', &
             '---', '---------', '-------------', '-------------', '-------------', '-------------'
      print_radial: do ipt=1,sd_nradial
        write (out,"(1x,i6,5(1x,g24.13e3))") ipt, sd_rtab(ipt), evec(ipt,ind,:)
      end do print_radial
      write (out,"()")
    end if
    !
    deallocate (evec,eval)
  end subroutine atomic_initial_wfn
  !
  subroutine prepare_initial_wavefunction
    real(rk)    :: ram
    complex(rk) :: norm(2)
    integer(ik) :: alloc
    !
    call TimerStart('Initial wavefunction')
    ram = 2*2*rk_bytes()*real(sd_nradial,kind=rk)*sd_nspin*real(sd_lmax,kind=rk)*real(sd_mmax-sd_mmin+1,kind=rk)
    write (out,"('Wavefunction requires ',f12.3,' Mbytes of memory')") ram/(1024._rk**2)
    call flush_wrapper(out)
    allocate (wfn_l%wfn(sd_nradial,sd_nspin,0:sd_lmax,sd_mmin:sd_mmax), &
              wfn_r%wfn(sd_nradial,sd_nspin,0:sd_lmax,sd_mmin:sd_mmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('Allocation of wavefunction array failed, code = ',i0)") alloc
      stop 'spherical_tdse%prepare_initial_wavefunction - out of memory'
    end if
    !
    !  Wavefunction initialization happens on the master node. The wavefunction is broadcast
    !  to other nodes if necessary.
    !
    if (nts%this_node==1) then
      select case (initial_wfn)
        case default
          write (out,"('Initial wavefunction choice ',a,' is not recognized')") trim(initial_wfn)
          stop 'spherical_tdse%prepare_initial_wavefunction - bad initial_wfn'
        case ('random')
          call random_initial_wfn('left',wfn_l)
          call random_initial_wfn('right',wfn_r)
        case ('unit')
          wfn_l%wfn = 1
          wfn_r%wfn = 1
        case ('atomic','single')
          call atomic_initial_wfn
        case ('read')
          call fetch_wavefunctions(initial_wfn_file)
      end select
    end if
    call nt_broadcast(wfn_l)
    call nt_broadcast(wfn_r)
    !
    !  If we are running adaptive Lmax, we do not yet know the appropriate lmax value here.
    !  Before we can figure it out, we need to normalize
    !
    wfn_l%lmax     = sd_lmax
    wfn_r%lmax     = sd_lmax
    ! sd_tolerance_* arrays are in the right/left order
    wfn_l%sd_tol_l = sd_tolerance_l(2)
    wfn_l%sd_tol_r = sd_tolerance_r(2)
    wfn_r%sd_tol_l = sd_tolerance_l(1)
    wfn_r%sd_tol_r = sd_tolerance_r(1)
    !
    call wt_normalize(wfn_l,wfn_r,norm)
    !
    !  The next subroutine is not parallelized, and runs on the master!
    !  This is inefficient, but is only done once.
    !
    call nt_merge_all(wfn_l)
    call nt_merge_all(wfn_r)
    if (nts%this_node==1) then
      call wt_reset_lrmax(wfn_l)
      call wt_reset_lrmax(wfn_r)
      wfn_l%lmax_top = wfn_l%lmax
      wfn_r%lmax_top = wfn_r%lmax
    end if
    call nt_broadcast(wfn_l)
    call nt_broadcast(wfn_r)
    !
    write (out,"('        Initial wavefunction norm was ',g24.13e3,1x,g24.13e3)") norm(1)
    write (out,"('Right wavefunction Cartesian norm was ',g24.13e3)") real(norm(2),kind=rk)
    write (out,"('                  Left/Right Lmax was ',2i4)") wfn_l%lmax, wfn_r%lmax
    write (out,"('               Left/Right nradial was ',2i8)") wfn_l%nradial, wfn_r%nradial
    call TimerStop('Initial wavefunction')
  end subroutine prepare_initial_wavefunction
  !
  subroutine fill_vpot_table
    integer(ik) :: its, alloc
    real(xk)    :: time, vp, th, ph
    !
    call TimerStart('Prepare VP table')
    allocate (vpot_table(0:3,-1:2*timesteps+1),stat=alloc)
    if (alloc/=0) stop 'spherical_tdse%fill_vpot_table - allocation failed'
    if (vp_shape=='table') then
      open (iu_temp,form='formatted',recl=256,action='read',position='rewind',status='old',file=trim(vp_table))
      read (iu_temp,*) vpot_table
      close (iu_temp)
    else
      vp = vp_apot(0._xk) ! Initialize internal structures of vp_apot
      !
      !  Start by filling the table 
      !
      time_steps: do its=-1,2*timesteps+1
        time = dt*0.5_xk*its
        vp   = vp_apot(time,theta=th,phi=ph)
        vpot_table(:,its) = (/ time, vp, th, ph /)
      end do time_steps
    end if
    !
    !  Issue a warning if vpot_table at its=-1,0 and 2*timesteps is not zero
    !
    if (any(vpot_table(1:3,(/-1,0,2*timesteps,2*timesteps+1/))/=0._rk)) then
      write (out,"(/'WARNING: Vector potential at time steps <= 0 and/or the end of the simulation is not zero')")
      if (.not. vp_as_is) then
        write (out,"( 'WARNING: Resetting vector-potential to zero, and initial/final orientation to lab.'/)")
      else
        write (out,"( 'WARNING: USING VECTOR-POTENTIAL AS IS. THE RESULTS ARE LIKELY UNPHYSICAL.'/)")
      end if
    end if
    !
    !  Force vector-potential to be zero at its=-1,0 and 2*timesteps
    !  Force local coordinate system at these time steps to be the laboratory system
    !
    if (.not. vp_as_is) then
      vpot_table(1:3,(/-1,0,2*timesteps,2*timesteps+1/)) = 0._rk
    end if
    !
    call TimerStop('Prepare VP table')
  end subroutine fill_vpot_table
  !
  subroutine unwrap_vpot_table
    integer(ik) :: its
    real(xk)    :: vp, th, ph, vp_extrap, th_ref, ph_ref
    real(xk)    :: th_v1, ph_v1, th_v2, ph_v2, r2_v1, r2_v2
    real(xk)    :: vp_max
    !
    if (.not.field_unwrap .or. rotation_mode=='none') then
      write (out,"(/'Skipping vector-potential unwrapping'/)")
      return
    end if
    call TimerStart('Unwrap VP table')
    !
    !  Spherical representation is not necessarily continuous in the normal parameter domain
    !  vp [0:inf), th [0:pi], ph [0:2pi)
    !  Our propagator is much happier when all three parameters to be continuos though the 
    !  first order (althogh it should be able to handle rapid change in theta and phi), 
    !  so that we need to unwrap the v.p.
    !
    !  Step 1: make sure V.P. magnitude's derivative is continuous. Points -1 and 0 are
    !          assumed to be "good" already, and will be used to start the process.
    ! 
    vp_max = maxval(abs(vpot_table(1,:)))
    unwrap_magnitude: do its=1,2*timesteps+1
      vp_extrap = 2*vpot_table(1,its-1) - vpot_table(1,its-2)
      vp        =   vpot_table(1,its)
      ! Try to avoid messing with vector-potential when it is nearly zero
      if (abs(vp_extrap+vp)>abs(vp_extrap-vp) .or. (abs(vp)+abs(vp_extrap))<unwrap_threshold*vp_max) cycle unwrap_magnitude
      ! We get smoother vector-potential by flipping the magnitude; do it!
      vpot_table(1,its) = -vpot_table(1,its)
      vpot_table(2,its) = -pi + vpot_table(2,its)
    end do unwrap_magnitude
    !
    !  Step 2: keeping VP magnitude constant, try to choose the smoothest possible version of (theta,phi)
    !
    !   We have two main possibilities:
    !   1.  th =  th0 + 2pi*n ; ph = ph0 + 2pi*k
    !   2   th = -th0 + 2pi*n ; ph = ph0 + 2pi*k + pi
    !
    unwrap_theta_phi: do its=0,2*timesteps+1
      th     = vpot_table(2,its) 
      ph     = vpot_table(3,its) 
      th_ref = vpot_table(2,its-1)
      ph_ref = vpot_table(3,its-1)
      th_v1  = nearest_2pi( th,   th_ref)
      ph_v1  = nearest_2pi( ph,   ph_ref)
      th_v2  = nearest_2pi(-th,   th_ref)
      ph_v2  = nearest_2pi( ph+pi,ph_ref)
      r2_v1  = (th_v1-th_ref)**2 + (ph_v1-ph_ref)**2
      r2_v2  = (th_v2-th_ref)**2 + (ph_v2-ph_ref)**2
      if (r2_v1<=r2_v2) then
        vpot_table(2,its) = th_v1
        vpot_table(3,its) = ph_v1
      else
        vpot_table(2,its) = th_v2
        vpot_table(3,its) = ph_v2
      end if
    end do unwrap_theta_phi
    call TimerStop('Unwrap VP table')
    contains
    function nearest_2pi(x,ref) result(v)
      real(xk), intent(in) :: x    
      real(xk), intent(in) :: ref 
      real(xk)             :: v
      !
      v = x + 2*pi_xk*nint((ref-x)/(2*pi_xk),kind=ik)
    end function nearest_2pi
  end subroutine unwrap_vpot_table
  !
  !  We calculate the electric field using numerical differentiation of the 
  !  vector-potential. Since we require high numerical accuracy, we use
  !  implicit derivative expressions.
  !
  !  WARNING: This routine implicitly assumes that the vector-potential starts
  !  WARNING: at zero and ends at zero. It will produce garbage at the ends of
  !  WARNING: the interval if this condition is violated
  !
  subroutine fill_efield_table
    integer(ik)            :: its, alloc
    real(xk), allocatable  :: vpot_xyz(:,:)              ! Vector-potential in Cartesian laboratory coordinates
    real(xk), allocatable  :: dt(:)                      ! dt(i) is the timestep leading to the time at point (i)
    real(xk), allocatable  :: d1(:,:), m1(:,:), m1f(:,:) ! tri-diagonal matrices defining the implict gradient
    logical, allocatable   :: m1p(:)                     ! pivot table 
    real(xk), allocatable  :: tmp1(:,:), tmp2(:,:)       ! Temporary arrays
    real(xk), allocatable  :: scr(:,:)                   ! Scratch for the linear solver
    !
    call TimerStart('Prepare efield table')
    !
    !  efield_table will remain until the end of the run. Remaining arrays must be
    !  deallocated before leaving fill_efied_table.
    !
    allocate (efield_table(3,0:2*timesteps), &
              vpot_xyz(0:2*timesteps,3), dt(0:2*timesteps+1), &
              d1(0:2*timesteps,3), m1(0:2*timesteps,3), m1f(0:2*timesteps,m3d_dc_size), m1p(0:2*timesteps), &
              tmp1(0:2*timesteps,3), tmp2(0:2*timesteps,3), scr(0:2*timesteps,3*m3d_sc_size), &
              stat=alloc)
    if (alloc/=0) stop 'spherical_tdse%fill_efield_table - allocation failed'
    !
    !  Calculate Cartesian vector-potential using spherical vector-potential in vp_apot()
    !  Also build the table of the timesteps and the rotation matrices.
    !
    fill_vpot_xyz: do its=0,2*timesteps
      vpot_xyz(its,:) = vpot_sph2xyz(vpot_table(1:3,its))
      dt(its)         = vpot_table(0,its) - vpot_table(0,its-1)
    end do fill_vpot_xyz
    dt(2*timesteps+1) = vpot_table(0,2*timesteps+1)-vpot_table(0,2*timesteps)
    !
    !  We need to prepare the Delta_1 and M_1 matrices for the first derivative.
    !  The code snipped below is lifted from initialize_radial_gradient() in spherical_data.f90
    !  We can't just reuse initialize_radial_gradient, since the boundary conditions are different
    !
    grad_tables: do its=0,2*timesteps
      ! Diagonal
      d1(its,1) = 1/dt(its) - 1/dt(its+1) + (-dt(its) + dt(its+1))/ (dt(its)**2 + dt(its)*dt(its+1) + dt(its+1)**2)
      m1(its,1) = (dt(its) + dt(its+1))**2/ (2*(dt(its)**2 + dt(its)*dt(its+1) + dt(its+1)**2))
      if (its>=2*timesteps) cycle grad_tables
      ! Sub-diagonal
      d1(its,2) = -((dt(its+2)**2*(2*dt(its+1) + dt(its+2)))/ (dt(its+1)*(dt(its+1) + dt(its+2))* &
                         (dt(its+1)**2 + dt(its+1)*dt(its+2) + dt(its+2)**2)))
      m1(its,2) = dt(its+2)**2/(2*(dt(its+1)**2 + dt(its+1)*dt(its+2) + dt(its+2)**2))
      ! Super-diagonal
      d1(its,3) = (dt(its)**2*(dt(its) + 2*dt(its+1)))/(dt(its+1)* (dt(its) + dt(its+1))*(dt(its)**2 + &
                         dt(its)*dt(its+1) + dt(its+1)**2))
      m1(its,3) = dt(its)**2/(2*(dt(its)**2 + dt(its)*dt(its+1) + dt(its+1)**2))
    end do grad_tables
    call m3d_decompose_x(m1,m1f,m1p)
    !
    !  We are done with all the matrices; the derivative is now given by:
    !
    !    M1^-1 Delta_1 A
    !
    call m3d_multiply_x(d1,vpot_xyz,tmp1)
    call m3d_solve_x(m1,m1f,m1p,tmp1,tmp2,scr)
    efield_table = -transpose(tmp2)
    !
    deallocate (vpot_xyz,dt,d1,m1,m1f,m1p,tmp1,tmp2,scr)
    call TimerStop('Prepare efield table')
  end subroutine fill_efield_table
  !
  !  Differentiate vector-potential without making assumptions about its
  !  values outside the time domain. Less accurate than fill_efield_table()
  !
  subroutine fill_efield_table_single_sided
    integer(ik)            :: its, alloc
    integer(ik)            :: il, ih
    real(xk), allocatable  :: vpot_xyz(:,:)              ! Vector-potential in Cartesian laboratory coordinates
    !
    call TimerStart('Prepare efield table (special)')
    !
    !  efield_table will remain until the end of the run. Remaining arrays must be
    !  deallocated before leaving fill_efied_table.
    !
    allocate (efield_table(3,0:2*timesteps), vpot_xyz(3,-1:2*timesteps+1), &
              stat=alloc)
    if (alloc/=0) stop 'spherical_tdse%fill_efield_table_single_sided - allocation failed'
    !
    !  Calculate Cartesian vector-potential using spherical vector-potential in vp_apot()
    !  Also build the table of the timesteps.
    !
    fill_vpot_xyz: do its=-1,2*timesteps+1
      vpot_xyz(:,its) = vpot_sph2xyz(vpot_table(1:3,its))
    end do fill_vpot_xyz
    !
    differentiate_vpot: do its=0,2*timesteps
      il = max(its-2,-1)
      ih = min(its+2,2*timesteps+1)
      efield_table(:,its) = - poly_gradient(vpot_table(0,il:ih),vpot_xyz(:,il:ih),its-il+1)
    end do differentiate_vpot
    !
    deallocate (vpot_xyz)
    call TimerStop('Prepare efield table (special)')
  end subroutine fill_efield_table_single_sided
  !
  !  Evaluate gradient at a grid point from the first derivative of the Lagrange interpolant
  !
  function poly_gradient(x,f,pos) result (g)
    real(xk), intent(in)    :: x(  :)             ! Time grid
    real(xk), intent(in)    :: f(:,:)             ! (Vector) function on a time frid
    integer(ik), intent(in) :: pos                ! Time point where we need the gradient
    real(xk)                :: g(size(f,dim=1))   ! Time derivative of the function
    !
    integer(ik) :: npts
    integer(ik) :: i, j
    real(xk)    :: term(size(f,dim=1))
    !
    npts = size(x)
    if (size(f,dim=1)<=0)    stop 'spherical_tdse%poly_gradient - bad vector length'
    if (size(f,dim=2)/=npts) stop 'spherical_tdse%poly_gradient - inconsistent sizes'
    if (pos<1 .or. pos>npts) stop 'spherical_tdse%poly_gradient - bad pos'
    !
    g(:) = 0
    accumulate_gradient: do i=1,npts
      if (i==pos) cycle accumulate_gradient
      term = f(:,i) - f(:,pos)
      accumulate_weights: do j=1,npts
        if (j==i) cycle accumulate_weights
        term = term / (x(i)-x(j))
        if (j==pos) cycle accumulate_weights
        term = term * (x(pos)-x(j))
      end do accumulate_weights
      g = g + term
    end do accumulate_gradient
  end function poly_gradient
  !
  function vpot_sph2xyz(atp)
    real(xk), intent(in) :: atp(3)
    real(xk)             :: vpot_sph2xyz(3)
    !
    vpot_sph2xyz(1) = atp(1)*sin(atp(2))*cos(atp(3))
    vpot_sph2xyz(2) = atp(1)*sin(atp(2))*sin(atp(3))
    vpot_sph2xyz(3) = atp(1)*cos(atp(2))
  end function vpot_sph2xyz
  !
  subroutine preview_laser_field
    integer(ik) :: its
    real(xk)    :: tatp(0:3) ! (time, a, theta, phi)
    real(xk)    :: efield(3) ! (Ex,Ey,Ez)
    !
    if (field_preview==' ') return
    call TimerStart('Dump VP table')
    !
    open (iu_temp,form='formatted',action='write',position='rewind',status='replace',file=trim(field_preview))
    !
    write (iu_temp,"(('#',a11,1x,a26,3(1x,a30,4x),1x,3(1x,a20,4x)))") &
           ' step ', ' Time, au[t] ', ' V.P. magnitude ', ' V.P. theta ', ' V.P. phi ', ' Ex ', ' Ey ', ' Ez ', &
           '------', '-------------', '----------------', '------------', '----------', '----', '----', '----'
    time_steps: do its=0,2*timesteps
      tatp   = vpot_table(:,its)
      efield = efield_table(:,its)
      write (iu_temp,"(1x,f11.1,1x,f26.16,3(1x,g30.19e3),1x,3(1x,g24.13e3))") 0.5_rk*its, tatp, efield
    end do time_steps
    !
    close (iu_temp)
    call TimerStop('Dump VP table')
  end subroutine preview_laser_field
  !
  subroutine visualize_wavefunctions(prefix,timestep,th,ph)
    character(len=*), intent(in) :: prefix
    integer(ik), intent(in)      :: timestep
    real(xk), intent(in)         :: th, ph    ! Rotation angles for the current localm coordinate system
    !
    integer(ik)          :: lval, mval, sval
    character(len=clen)  :: filename
    real(rk)             :: rm(3,3)
    complex(rk)          :: l_lap(sd_nradial), l_grad(sd_nradial)
    complex(rk)          :: r_lap(sd_nradial), r_grad(sd_nradial)
    complex(rk)          :: tmp(sd_nradial), scr(sd_nradial,m3d_sc_size)
    !
    if (prefix==' ') return
    !
    call TimerStart('Visualize wavefunction')
    !
    if (skip_left_propagation .and. fake_left_propagation) then
      call wt_fake_left(wfn_l,wfn_r)
    end if
    call nt_merge_all(wfn_l)
    call nt_merge_all(wfn_r)
    if (nts%this_node/=1) then
      call TimerStop('Visualize wavefunction')
      return
    end if
    write (filename,"(a,i10.10,'.data')") trim(prefix), visualize_1stseq
    visualize_1stseq = visualize_1stseq + 1
    open(iu_temp,form='formatted',recl=256,action='write',position='rewind',status='replace',file=trim(filename))
    write (iu_temp,"(6(1x,i0))") sd_nradial, sd_nspin, sd_lmax, sd_mmin, sd_mmax, timestep
    !
    call MathRotationMatrix((/real(ph,kind=rk),real(th,kind=rk),0._rk/),rm)
    write (iu_temp,"(3(1x,g27.18e3))") rm(1,:)
    write (iu_temp,"(3(1x,g27.18e3))") rm(2,:)
    write (iu_temp,"(3(1x,g27.18e3))") rm(3,:)
    !
    write (iu_temp,"(5(1x,g22.13e3))") sd_rtab(1:sd_nradial)
    !
    dump_m_channels: do mval=sd_mmin,sd_mmax
      dump_l_channels: do lval=0,sd_lmax
        dump_s_channels: do sval=1,sd_nspin
          if (lval==0) then
            !$omp parallel sections default(none) &
            !$omp&  shared(sd_d2n_l0,sd_m2n_l0,sd_m2np_l0,sd_m2nf_l0) &
            !$omp&  shared(sd_d1n_l0,sd_m1n_l0,sd_m1np_l0,sd_m1nf_l0) &
            !$omp&  shared(wfn_l,wfn_r,sval,lval,mval,l_lap,r_lap,l_grad,r_grad) &
            !$omp&  private(tmp,scr)
            ! Left laplacian
            !$omp section
            call m3d_multiply(sd_d2n_l0,wfn_l%wfn(:,sval,lval,mval),tmp)
            call m3d_solve(sd_m2n_l0,sd_m2nf_l0,sd_m2np_l0,tmp,l_lap,scr)
            ! Right laplacian
            !$omp section
            call m3d_multiply(sd_d2n_l0,wfn_r%wfn(:,sval,lval,mval),tmp)
            call m3d_solve(sd_m2n_l0,sd_m2nf_l0,sd_m2np_l0,tmp,r_lap,scr)
            ! Left gradient
            !$omp section
            call m3d_multiply(sd_d1n_l0,wfn_l%wfn(:,sval,lval,mval),tmp)
            call m3d_solve(sd_m1n_l0,sd_m1nf_l0,sd_m1np_l0,tmp,l_grad,scr)
            ! Right gradient
            !$omp section
            call m3d_multiply(sd_d1n_l0,wfn_r%wfn(:,sval,lval,mval),tmp)
            call m3d_solve(sd_m1n_l0,sd_m1nf_l0,sd_m1np_l0,tmp,r_grad,scr)
            !$omp end parallel sections
          else ! lval /= 0
            !$omp parallel sections default(none) &
            !$omp&  shared(sd_d2n_lx,sd_m2n_lx,sd_m2np_lx,sd_m2nf_lx) &
            !$omp&  shared(sd_d1n_lx,sd_m1n_lx,sd_m1np_lx,sd_m1nf_lx) &
            !$omp&  shared(wfn_l,wfn_r,sval,lval,mval,l_lap,r_lap,l_grad,r_grad) &
            !$omp&  private(tmp,scr)
            ! Left laplacian
            !$omp section
            call m3d_multiply(sd_d2n_lx,wfn_l%wfn(:,sval,lval,mval),tmp)
            call m3d_solve(sd_m2n_lx,sd_m2nf_lx,sd_m2np_lx,tmp,l_lap,scr)
            ! Right laplacian
            !$omp section
            call m3d_multiply(sd_d2n_lx,wfn_r%wfn(:,sval,lval,mval),tmp)
            call m3d_solve(sd_m2n_lx,sd_m2nf_lx,sd_m2np_lx,tmp,r_lap,scr)
            ! Left gradient
            !$omp section
            call m3d_multiply(sd_d1n_lx,wfn_l%wfn(:,sval,lval,mval),tmp)
            call m3d_solve(sd_m1n_lx,sd_m1nf_lx,sd_m1np_lx,tmp,l_grad,scr)
            ! Right gradient
            !$omp section
            call m3d_multiply(sd_d1n_lx,wfn_r%wfn(:,sval,lval,mval),tmp)
            call m3d_solve(sd_m1n_lx,sd_m1nf_lx,sd_m1np_lx,tmp,r_grad,scr)
            !$omp end parallel sections
          end if
          write (iu_temp,"(5(1x,g22.13e3))") (1._rk/tsurf%wfn_scale) * wfn_l%wfn(1:sd_nradial,sval,lval,mval)
          write (iu_temp,"(5(1x,g22.13e3))") (1._rk/tsurf%wfn_scale) * l_grad
          write (iu_temp,"(5(1x,g22.13e3))") (1._rk/tsurf%wfn_scale) * l_lap
          write (iu_temp,"(5(1x,g22.13e3))") tsurf%wfn_scale * wfn_r%wfn(1:sd_nradial,sval,lval,mval)
          write (iu_temp,"(5(1x,g22.13e3))") tsurf%wfn_scale * r_grad
          write (iu_temp,"(5(1x,g22.13e3))") tsurf%wfn_scale * r_lap
        end do dump_s_channels
      end do dump_l_channels
    end do dump_m_channels
    close(iu_temp)
    call TimerStop('Visualize wavefunction')
  end subroutine visualize_wavefunctions
  !
  subroutine dump_wavefunctions(prefix)
    character(len=*), intent(in) :: prefix
    integer(ik)                  :: lval, mval, sval, ir, ios
    character(len=clen)          :: filename
    !
    call TimerStart('Dump wavefunction')
    dump_l_channels: do lval=0,sd_lmax
      dump_m_channels: do mval=max(sd_mmin,-lval),min(lval,sd_mmax)
        if (mval>=0) then
          write (filename,"(a,'-L',i0.3,'-M+',i0.4)") trim(prefix), lval, mval
        else
          write (filename,"(a,'-L',i0.3,'-M-',i0.4)") trim(prefix), lval, -mval
        end if
        open(iu_temp,form='formatted',recl=256,action='write',position='rewind',status='replace',file=trim(filename),iostat=ios)
        if (ios/=0) then
          write (out,"('WARNING: File ',a,' cannot be opened for writing. Code = ',i0)") trim(filename), ios
          cycle dump_m_channels
        end if
        write (iu_temp,"('#',a2,1x,a24,2(2x,a34,1x,a34))") &
              'S', ' R, Bohr ', ' Re[left wfn] ', ' Im[left wfn] ', ' Re[right wfn] ', ' Im[right wfn] '
        dump_s_channels: do sval=1,sd_nspin
          dump_radial: do ir=1,sd_nradial
            write (iu_temp,"(1x,i2,1x,g24.13,2(2x,g34.21e3,1x,g34.21e3))") &
                   sval, sd_rtab(ir), wfn_l%wfn(ir,sval,lval,mval), wfn_r%wfn(ir,sval,lval,mval)
          end do dump_radial
        end do dump_s_channels
        close(iu_temp)
      end do dump_m_channels
    end do dump_l_channels
    call TimerStop('Dump wavefunction')
  end subroutine dump_wavefunctions
  !
  subroutine fetch_wavefunctions(prefix)
    character(len=*), intent(in) :: prefix
    integer(ik)                  :: lval, mval, sval, ir
    character(len=clen)          :: filename
    character(len=1)             :: buf
    integer(ik)                  :: line, ios
    integer(ik)                  :: stmp
    real(rk)                     :: rtmp(5)
    !
    call TimerStart('Fetch wavefunction')
    fetch_l_channels: do lval=0,sd_lmax
      fetch_m_channels: do mval=max(sd_mmin,-lval),min(lval,sd_mmax)
        if (mval>=0) then
          write (filename,"(a,'-L',i0.3,'-M+',i0.4)") trim(prefix), lval, mval
        else
          write (filename,"(a,'-L',i0.3,'-M-',i0.4)") trim(prefix), lval, -mval
        end if
        open(iu_temp,form='formatted',recl=256,action='read',position='rewind',status='old',file=trim(filename))
        line = 1
        read (iu_temp,"(a1)",iostat=ios) buf
        if (ios/=0) then
          write (out,"('Error ',i0,' skipping line ',i0,' of file ',a)") ios, line, trim(filename)
          stop 'spherical_tdse%fetch_wavefunctions - skip error'
        end if
        fetch_s_channels: do sval=1,sd_nspin
          fetch_radial: do ir=1,sd_nradial
            line = line + 1
            read (iu_temp,*,iostat=ios) stmp, rtmp
            if (ios/=0) then
              write (out,"('Error ',i0,' reading line ',i0,' of file ',a)") ios, line, trim(filename)
              stop 'spherical_tdse%fetch_wavefunctions - read error'
            end if
            if (stmp/=sval .or. abs(sd_rtab(ir)-rtmp(1))>1e-10_rk*max(1._rk,sd_rtab(ir))) then
              write (out,"('Grid mismatch reading line ',i0,' of file ',a)") line, trim(filename)
              stop 'spherical_tdse%fetch_wavefunctions - data error'
            end if
            wfn_l%wfn(ir,sval,lval,mval) = cmplx(rtmp(2),rtmp(3),kind=rk)
            wfn_r%wfn(ir,sval,lval,mval) = cmplx(rtmp(4),rtmp(5),kind=rk)
          end do fetch_radial
        end do fetch_s_channels
        close(iu_temp)
      end do fetch_m_channels
    end do fetch_l_channels
    call TimerStop('Fetch wavefunction')
  end subroutine fetch_wavefunctions
  !
  subroutine lab_dipole(wfn_l,wfn_r,rot,afield,efield,norm,dipole,velocity,acceleration)
    type(sd_wfn), intent(in) :: wfn_l           ! Left wavefunction 
    type(sd_wfn), intent(in) :: wfn_r           ! Right wavefunction 
    real(rk), intent(in)     :: rot(3,3)        ! Rotation matrix (lab->local)
    real(rk), intent(in)     :: afield(3)       ! Vector-potential, needed for the dipole velocity
    real(rk), intent(in)     :: efield(3)       ! Electric field, needed for the dipole acceleration
    complex(rk), intent(in)  :: norm            ! Wavefunction norm, ditto
    complex(rk), intent(out) :: dipole(3)       ! <L|q r|R> expectation value in the lab frame
    complex(rk), intent(out) :: velocity(3)     ! (d/d t) <L|q r|R> in the lab frame
    complex(rk), intent(out) :: acceleration(3) ! (d^2/d t^2) <L|q r|R> in the lab frame
    !
    complex(rk) :: loc_dipole(3,3)     ! Dipole and its derivatives in the local coordinate system
    !
    call wt_dipole(wfn_l,wfn_r,do_dipole,loc_dipole)
    !
    !  Rotate everything back into the lab system
    !
    !  gfortran will generate a temporary array below, promotiong rot to complex(rk).
    !  that matmul implementation really sucks ...
    ! 
    dipole       = matmul(transpose(rot),loc_dipole(:,1))
    velocity     = matmul(transpose(rot),loc_dipole(:,2))
    acceleration = matmul(transpose(rot),loc_dipole(:,3))
    !
    !  Add operator-derivative terms
    !
    if (do_dipole_plasma) then
      if (do_dipole(2)) velocity     = velocity     - (electron_charge**2/electron_mass) * afield
      if (do_dipole(3)) acceleration = acceleration + (electron_charge**2/electron_mass) * efield
    else
      if (do_dipole(2)) velocity     = velocity     - (electron_charge**2/electron_mass) * afield * norm
      if (do_dipole(3)) acceleration = acceleration + (electron_charge**2/electron_mass) * efield * norm
    end if
  end subroutine lab_dipole
  !
  subroutine subdivide_timestep(its,lmax,nsub,vsub)
    integer(ik), intent(in)  :: its       ! "Macro" time step
    integer(ik), intent(in)  :: lmax      ! Largest angular momentum on this time step
    integer(ik), intent(out) :: nsub      ! Number of "small" time steps after subdivision 
    real(xk), allocatable    :: vsub(:,:) ! Subdivided time step parameters
    !
    real(xk)                 :: cdt       ! The "big" time (half-)step
    real(xk)                 :: sdt       ! The "small" time (half-)step
    real(xk)                 :: vp        ! Vector-potential at the current time step
    real(xk)                 :: ladt      ! The value of (lmax+1)*vector_potential*dt for the "big" step
    integer(ik)              :: alloc
    integer(ik)              :: isub, ic
    !
    cdt = 0.5_rk*(vpot_table(0,2*its+2)-vpot_table(0,2*its))
    vp  = abs(vpot_table(1,2*its+1))
    !
    select case (dt_subdivision)
      case default
        write (out,"('subdivide_timestep: Unrecognized value of dt_subdivision: ',a)") trim(dt_subdivision)
        stop 'spherical_tdse%subdivide_timestep - bad dt_subdivision'
      case ('off')
        ladt = dt_max_la
      case ('lmax','on')
        ladt = abs((lmax+1)*dt)
      case ('lmax-a')
        ladt = abs((lmax+1)*vp**dt)
      case ('lmax-a2')
        ladt = abs((lmax+1)*vp**2*dt)
    end select
    if (ladt<=dt_max_la) then
      nsub = 1
    else
      nsub = ceiling(ladt/dt_max_la)
    end if
    if (nsub>max_subdivision) then 
      max_subdivision = nsub
      if (verbose>=0) write (out,"('Maximum of the time step subdivision factor increased to ',i0)") max_subdivision
    end if
    if (nsub>1) timesteps_subdivided = timesteps_subdivided + 1
    timesteps_micro = timesteps_micro + nsub
    !
    !  Prepare vector-potential array as needed.
    !
    if (allocated(vsub)) then
      if (ubound(vsub,dim=2)<2*nsub) deallocate (vsub) ! Spurious gfortran warning here
    end if
    if (.not.allocated(vsub)) then
      allocate (vsub(0:3,0:2*nsub),stat=alloc)
      if (alloc/=0) then
        write (out,"('spherical_tdse%subdivide_timestep: allocation for nsub = ',i0,' failed. Error = ',i0)") nsub, alloc
        call flush_wrapper(out)
        stop 'spherical_tdse%subdivide_timestep - alloc'
      end if
    end if
    !
    !  Copy vector-potential points at the beginning and end of the interval.
    !  For the remaining points, we'll need to interpolate.
    !
    sdt = cdt / nsub
    vsub(:,0) = vpot_table(:,2*its)
    fill_vsub: do isub=1,2*nsub-1
      vsub(0,isub) = vsub(0,0) + sdt*isub
      fill_components: do ic=1,3
        vsub(ic,isub) = interpolate_vpot(ic,vsub(0,isub),its)
      end do fill_components
    end do fill_vsub
    vsub(:,2*nsub) = vpot_table(:,2*its+2)
    !
    if (verbose>=3) then
      write (out,"(/(2x,a8,4(1x,a20)))") &
                       ' TS ', ' T ', ' A ', ' Theta ', ' Phi ', &
                       '----', '---', '---', '-------', '-----'
      debug_print_subdivide: do isub=0,2*nsub
        if (isub==0*nsub) write (out,"('R ',i8,4(1x,g20.12e3))") 2*its+0, vpot_table(:,2*its+0)
        if (isub==1*nsub) write (out,"('R ',i8,4(1x,g20.12e3))") 2*its+1, vpot_table(:,2*its+1)
        if (isub==2*nsub) write (out,"('R ',i8,4(1x,g20.12e3))") 2*its+2, vpot_table(:,2*its+2)
                          write (out,"('S ',i8,4(1x,g20.12e3))") isub, vsub(:,isub)
      end do debug_print_subdivide
      write (out,"()")
    end if
  end subroutine subdivide_timestep
  !
  function interpolate_vpot(ic,tim,its) result(v)
    integer(ik), intent(in) :: ic  ! Component of vpot_table to interpolate: vpot_table(ic,:)
    real(xk), intent(in)    :: tim ! Time point to interpolate at
    integer(ik), intent(in) :: its ! Reference point in vpot_table: vpot_table(:,its)
    real(xk)                :: v   ! Interpolated value
    !
    integer(ik) :: lp, rp ! Left and right points in vpot_table
    !
    lp = max(lbound(vpot_table,dim=2),2*its-dt_interpolant_width/2)
    rp = min(ubound(vpot_table,dim=2),lp+dt_interpolant_width-1)
    v  = MathInterpolate(tim,vpot_table(0,lp:rp),vpot_table(ic,lp:rp))
  end function interpolate_vpot
  !
  subroutine write_detail_header
    if (detail_output==' ') return
    !
    write (out,"(/'Saving detailed output to ',a/)") trim(detail_output)
    open(iu_detail,file=trim(detail_output),form='formatted',status='replace',position='rewind',recl=1050,pad='no')
    write (iu_detail,"('# Field Columns Data')")
    write (iu_detail,"('#  1    2,13    Timestep')")
    write (iu_detail,"('#  2   15,46    Time, au[t]')")
    write (iu_detail,"('#  3   48,79    Vector-potential magnitude')")
    write (iu_detail,"('#  4   81,112   Vector-potential, lab theta')")
    write (iu_detail,"('#  5  114,145   Vector-potential, lab phi')")
    write (iu_detail,"('#  6  147,178   Re[<Psi_L|Psi_R>]')")
    write (iu_detail,"('#  7  180,211   Im[<Psi_L|Psi_R>]')")
    write (iu_detail,"('#  8  213,244   Re[<Psi_L|H_at+H_L+V_cap|Psi_R>], Hartree')")
    write (iu_detail,"('#  9  246,277   Im[<Psi_L|H_at+H_L+V_cap|Psi_R>], Hartree')")
    write (iu_detail,"('# 10  279,310   Re[<Psi_L|H_at+H_L|Psi_R>], Hartree')")
    write (iu_detail,"('# 11  312,343   Im[<Psi_L|H_at+H_L|Psi_R>], Hartree')")
    write (iu_detail,"('# 12  345,376   Re[<Psi_L|e . x|Psi_R>], e-Bohr')")
    write (iu_detail,"('# 13  378,409   Im[<Psi_L|e . x|Psi_R>], e-Bohr')")
    write (iu_detail,"('# 14  411,442   Re[<Psi_L|e . y|Psi_R>], e-Bohr')")
    write (iu_detail,"('# 15  444,475   Im[<Psi_L|e . y|Psi_R>], e-Bohr')")
    write (iu_detail,"('# 16  477,508   Re[<Psi_L|e . z|Psi_R>], e-Bohr')")
    write (iu_detail,"('# 17  510,541   Im[<Psi_L|e . z|Psi_R>], e-Bohr')")
    write (iu_detail,"('# 18  543,574   Re[(d^2/d t^2) <Psi_L|e . x|Psi_R>], e-Bohr/jiffy^2')")
    write (iu_detail,"('# 19  576,607   Im[(d^2/d t^2) <Psi_L|e . x|Psi_R>], e-Bohr/jiffy^2')")
    write (iu_detail,"('# 20  609,640   Re[(d^2/d t^2) <Psi_L|e . y|Psi_R>], e-Bohr/jiffy^2')")
    write (iu_detail,"('# 21  642,673   Im[(d^2/d t^2) <Psi_L|e . y|Psi_R>], e-Bohr/jiffy^2')")
    write (iu_detail,"('# 22  675,706   Re[(d^2/d t^2) <Psi_L|e . z|Psi_R>], e-Bohr/jiffy^2')")
    write (iu_detail,"('# 23  708,739   Im[(d^2/d t^2) <Psi_L|e . z|Psi_R>], e-Bohr/jiffy^2')")
    write (iu_detail,"('# 24  741,772   Re[(d/d t) <Psi_L|e . x|Psi_R>], e-Bohr/jiffy')")
    write (iu_detail,"('# 25  774,805   Im[(d/d t) <Psi_L|e . x|Psi_R>], e-Bohr/jiffy')")
    write (iu_detail,"('# 26  807,838   Re[(d/d t) <Psi_L|e . y|Psi_R>], e-Bohr/jiffy')")
    write (iu_detail,"('# 27  840,871   Im[(d/d t) <Psi_L|e . y|Psi_R>], e-Bohr/jiffy')")
    write (iu_detail,"('# 28  873,904   Re[(d/d t) <Psi_L|e . z|Psi_R>], e-Bohr/jiffy')")
    write (iu_detail,"('# 29  906,937   Im[(d/d t) <Psi_L|e . z|Psi_R>], e-Bohr/jiffy')")
    write (iu_detail,"('# 30  939,970   Electric field, lab X, atomic units')")
    write (iu_detail,"('# 31  972,1003  Electric field, lab Y, atomic units')")
    write (iu_detail,"('# 32 1005,1036  Electric field, lab Z, atomic units')")
    if (sd_pot_nonlocal) then
      write (iu_detail,"('# WARNING: Non-local potential detected. Dipole acceleration results are incorrect')")
      write (out,      "(/ 'WARNING: Non-local potential detected. Dipole acceleration results are incorrect'/)")
    end if
    write (iu_detail,"('#',1x,a12,31(1x,a32))") &
           ' i ', ' time ', ' vp ', ' theta ', ' phi ', ' re(norm) ', ' im(norm) ', &
           ' re(energy) ', ' im(energy) ', ' re(en-nocap) ', ' im(en-nocap) ', &
           ' re(dip_x) ', ' im(dip_x) ', ' re(dip_y) ', ' im(dip_y) ', ' re(dip_z) ', ' im(dip_z) ', &
           ' re(acc_x) ', ' im(acc_x) ', ' re(acc_y) ', ' im(acc_y) ', ' re(acc_z) ', ' im(acc_z) ', &
           ' re(vel_x) ', ' im(vel_x) ', ' re(vel_y) ', ' im(vel_y) ', ' re(vel_z) ', ' im(vel_z) ', &
           ' F_x ', ' F_y ', ' F_z '
    write (iu_detail,"('#',1x,a12,31(1x,a32))") &
           ' 1 ', ' 2 ', ' 3 ', ' 4 ', ' 5 ', ' 6 ', ' 7 ', ' 8 ', ' 9 ', ' 10 ', ' 11 ', &
           ' 12 ', ' 13 ', ' 14 ', ' 15 ', ' 16 ', ' 17 ', &
           ' 18 ', ' 19 ', ' 20 ', ' 21 ', ' 22 ', ' 23 ', &
           ' 24 ', ' 25 ', ' 26 ', ' 27 ', ' 28 ', ' 29 ', &
           ' 30 ', ' 31 ', ' 32 '
  end subroutine write_detail_header
  !
  subroutine propagation
    integer(ik)           :: its
    integer(ik)           :: its_from          ! Time step (possibly half-step) from which we want to step with t-SURF
    integer(ik)           :: its_to            ! Time step (possibly half-step) to which we want to step with t-SURF
    integer(ik)           :: nsub, isub        ! Number of micro time steps after subdivision, and the counter
    real(xk), allocatable :: vsub(:,:)         ! Same as vpot_table(), but for the subdivided time step.
                                               ! The first index in 0=time, 1=signed magnitude, 2=theta, 3=phi
                                               ! The last index runs from 0 to 2*nsub.
                                               ! The array is maintaned by subdivide_timestep
    real(xk)              :: time              ! Time at the beginning of the macro time step
    real(rk)              :: afield(3)         ! The vector potential
    real(rk)              :: efield(3)         ! The electric field
    real(xk)              :: rdt, rdt2
    real(xk)              :: vp1, vp2          ! Magnitude of the vector-potential for the first and second half-steps
    real(xk)              :: th1, th2          ! vector-potential direction, ditto
    real(xk)              :: ph1, ph2          ! vector-potential direction, ditto
    complex(xk)           :: cdt1, cdt2        ! Time steps for the first and second half-steps
    logical               :: checkpoint_go
    !
    !  Propagation routine can handle varying time steps, at least in principle
    !
    if (detail_output==' ') then
      if (any(do_dipole(2:3))) then
        write (out,"(/'Detailed output is not enabled; will not calculate dipole velocity or acceleration'/)")
      end if
      do_dipole(2:3) = .false.
    end if
    if (skip_left_propagation) then
      if (fake_left_propagation) then
        write (out,"(/'WARNING: Left wavefunction approximated by conjugate of the right wavefunction.')")
        write (out,"( 'WARNING: Expect reduced accuracy for the observables during propagation.'/)")
      else
        write (out,"(/'WARNING: Left wavefunction is not evaluated. Observables are incorrect.'/)")
      end if
    end if
    if (nts%this_node==1) then
      call write_detail_header
    end if
    !
    call grid_resize(max(wfn_l%nradial,wfn_r%nradial))
    call nt_merge_all(wfn_r)
    if (.not.skip_left_propagation) then
      call nt_merge_all(wfn_l)
    end if
    call nt_rebalance(max(wfn_l%lmax,wfn_r%lmax))
    !
    time_steps: do its=0,timesteps-1
      time   = vpot_table(0,2*its)
      !
      !  vector-potential at the beginning of the macro step
      !
      vp1    = vpot_table(1,2*its)
      th1    = vpot_table(2,2*its)
      ph1    = vpot_table(3,2*its)
      !
      call ckpt_do_checkpoint(ckpt,checkpoint_go,its,vpot_table(0:3,2*its),(/out,iu_detail/),wfn_l,wfn_r,tsurf)
      if (.not.checkpoint_go) cycle time_steps ! Restart has been requested, but we are not at the point
                                               ! indicated in the checkpoint file yet.
      !
      !  Cartesian vector-potential and electric field are at the beginning of the time step
      !
      afield = real(vp1*(/sin(th1)*cos(ph1),sin(th1)*sin(ph1),cos(th1)/),kind=kind(afield))
      efield = real(efield_table(:,2*its),kind=kind(efield))
      call write_output(force=.false.)
      if (mod(its,visualize_each)==0) then
        call visualize_wavefunctions(visualize_prefix,its,th1,ph1)
      end if
      !
      !  See whether this time step needs to be subdivided to keep things stable
      !
      call subdivide_timestep(its,max(wfn_l%lmax,wfn_r%lmax),nsub,vsub)
      !
      !  The subdivided timestep.
      !
      subdivided_timestep: do isub=0,nsub-1
        cdt1 = 0.5_xk*(vsub(0,2*isub+1) - vsub(0,2*isub+0))
        cdt2 = 0.5_xk*(vsub(0,2*isub+2) - vsub(0,2*isub+1))
        !  vector-potential at the beginning of the micro step
        vp1    = vsub(1,2*isub)
        th1    = vsub(2,2*isub)
        ph1    = vsub(3,2*isub)
        !  vector-potential at the end of the micro step
        vp2    = vsub(1,2*isub+2)
        th2    = vsub(2,2*isub+2)
        ph2    = vsub(3,2*isub+2)
        !
        call nt_merge_borders(wfn_r)
        call pt_fwd_laser_n(wfn_r,vp1,cdt1)
        if (.not.skip_left_propagation) then
          call nt_merge_borders(wfn_l)
          call pt_fwd_laser_t(wfn_l,vp1,cdt1)
        end if
        !
        !  We are at the mid-step; apply atomic propagator for the first half-step
        !
        call pt_fwd_atomic_n(wfn_r,cdt1)
        call pt_rev_atomic_n(wfn_r,cdt1)
        if (.not.skip_left_propagation) then
          call pt_fwd_atomic_t(wfn_l,cdt1)
          call pt_rev_atomic_t(wfn_l,cdt1)
        end if
        !
        !  Rotation; nominally a full step. It will be broken into still smaller sub-steps
        !  if this turns out to be necessary.
        !
        call rotate(from=(/th1,ph1/),to=(/th2,ph2/))
        !
        !  Atomic propagator for the second half-step
        !
        call pt_fwd_atomic_n(wfn_r,cdt2)
        call pt_rev_atomic_n(wfn_r,cdt2)
        if (.not.skip_left_propagation) then
          call pt_fwd_atomic_t(wfn_l,cdt2)
          call pt_rev_atomic_t(wfn_l,cdt2)
        end if
        !
        !  Second half of the time step; vector potential is at the <i>next</i> time-step
        !
        call nt_merge_borders(wfn_r)
        call pt_rev_laser_n(wfn_r,vp2,cdt2)
        if (.not.skip_left_propagation) then
          call nt_merge_borders(wfn_l)
          call pt_rev_laser_t(wfn_l,vp2,cdt2)
        end if
        !
        call grid_adapt
        !
      end do subdivided_timestep
      !
      !  Node rebalancing is potentially an expensive operation!
      !
      if (nt_rebalance_needed(max(wfn_l%lmax,wfn_r%lmax))) then
        call nt_merge_all(wfn_r)
        if (.not.skip_left_propagation) then
          call nt_merge_all(wfn_l)
        end if
        call nt_rebalance(max(wfn_l%lmax,wfn_r%lmax))
      end if
      !
      !  Photoelectron spectrum accumulation with t-SURF. 
      !
      !  We do not necessarily wish to run t-SURF at each time step (it is expensive,
      !  compared to the propagator). However, we also do not want to overshoot the
      !  end of the simulation. As the result, we need to do a bit of fiddling with
      !  the effective t-SURF time step here.
      !
      !  The test below is deliberately for the time step index "its", not (its+1)
      !  where the actual time step will occur. We wish to perform (an asymmetric) 
      !  t-SURF time step at the very beginning of the simulation as well.
      !
      if (modulo(its,sts_step_frequency)==0) then
         !
         !  We are good to do a t-SURF time step centered at the (full) time step (its+1)
         !  The corresponding time and vector potential are at the offset (2*its+2) in
         !  vpot_table. 
         !
         its_from = max(          0,(2*its+2)-sts_step_frequency)  ! Don't step beyond beginning of the simulation
         its_to   = min(2*timesteps,(2*its+2)+sts_step_frequency)  ! Don't step beyond end of the simulation
         !
         !  If this is the last time step to trigger t-SURF, there is a possibility of missing 
         !  a partial time-step at the very end. Extend the time step to the end of the simulation.
         !
         if (its+sts_step_frequency>timesteps-1) its_to = 2*timesteps
         !
         !  Now we can calculate the duration of the t-SURF time step, and finally do it. 
         !  Because we try to use large time steps in t-SURF, it is not a good idea to assume
         !  that vector-potential remains constant; thankfully, we know the time derivative
         !  of the vector-potential (the electric field; except for the sign change).
         !
         rdt  = vpot_table(0,2*its+2)-vpot_table(0,its_from)  ! Time step to the expansion point
         rdt2 = vpot_table(0,its_to) -vpot_table(0,2*its+2)   ! Time step to the end of the interval
         call sts_timestep(wfn_l,wfn_r,tsurf,rdt,rdt2,vpot_table(1:3,2*its+2),efield_table(1:3,2*its+2))
      end if
    end do time_steps
    if (allocated(vsub)) deallocate (vsub)
    !
    time   = vpot_table(0,2*its)
    vp1    = vpot_table(1,2*its)
    th1    = vpot_table(2,2*its)
    ph1    = vpot_table(3,2*its)
    afield = real(vp1*(/sin(th1)*cos(ph1),sin(th1)*sin(ph1),cos(th1)/),kind=kind(afield))
    efield = real(efield_table(:,2*its),kind=kind(efield))
    !
    call write_output(force=.true.)
    call visualize_wavefunctions(visualize_prefix,its,th1,ph1)
    !
    if (detail_output/=' ' .and. nts%this_node==1) then
      close(iu_detail)
    end if
    !
    !  Rotate wavefunction back into the lab orientation for analysis
    !
    call rotate(from=vpot_table(2:3,2*its),to=(/0._xk,0._xk/))
    if (rotation_mode/='none') then
      write (out,"(/'Wavefunction restored back to laboratory frame'/)")
    end if
    !
    call grid_resize(sd_nradial_max)
    !
    contains 
    !
    subroutine write_output(force)
      logical, intent(in) :: force
      !
      real(rk)            :: abs_dipole
      complex(rk)         :: energy(2), norm
      complex(rk)         :: dipole(3)
      complex(rk)         :: velocity(3)
      complex(rk)         :: acceleration(3)
      real(rk)            :: rm(3,3)
      logical             :: normal_out     ! Produce output on the standard output
      logical             :: detail_out     ! Produce output on the detailed output
      !
      normal_out = force .or. mod(its,output_each)==0
      detail_out = (detail_output/=' ') .and. (force .or. mod(its,detail_frequency)==0)
      if (.not.normal_out .and. .not.detail_out) return
      !
      if (skip_left_propagation .and. fake_left_propagation) then
        call wt_fake_left(wfn_l,wfn_r)
      end if
      call nt_merge_borders(wfn_l)
      call nt_merge_borders(wfn_r)
      call wt_energy(wfn_l,wfn_r,vp1,energy,norm)
      call MathRotationMatrix(real((/ph1,th1,0._xk/),kind=kind(rm)),rm)
      call lab_dipole(wfn_l,wfn_r,rm,afield,efield,norm,dipole,velocity,acceleration) 
      !
      if (detail_out .and. nts%this_node==1) then
        write (iu_detail,"(1x,i12,31(1x,g32.22e4))") &
               its, time, vp1, th1, ph1, norm, energy, dipole, acceleration, velocity, efield
      end if
      !
      abs_dipole = sqrt(sum(abs(dipole)**2))
      if (normal_out) then
        write (out,"('@ ',i9,' t= ',f12.4,' a= ',f14.6,' th= ',f10.3,' ph= ',f10.3,' <l|r>= ',g23.12,1x,g17.6," // &
                   "' <l|h|r>= ',g23.12,1x,g17.6,' |d|= ',g20.9)") &
               its, time, vp1, th1, ph1, norm, energy(1), abs_dipole
        call flush_wrapper(out)
      end if
    end subroutine write_output
    !
    subroutine rotate(from,to)
      real(xk), intent(in) :: from(:), to(:) ! Initial and final field orientation
      !
      select case (rotation_mode)
        case default
          write (out,"('spherical_tdse%propagation: Rotation mode ',a,' is not recognized')") trim(rotation_mode)
          stop 'spherical_tdse%propagation - bad rotation_mode'
        case ('none')
        case ('sparse')
          call rt_rotate(wfn_r,from=from,to=to,left=.false.)
          if (.not.skip_left_propagation) then
            call rt_rotate(wfn_l,from=from,to=to,left=.true.)
          end if
        case ('brute force')
          call rt_rotate_bruteforce(wfn_r,from=from,to=to,left=.false.)
          if (.not.skip_left_propagation) then
            call rt_rotate_bruteforce(wfn_l,from=from,to=to,left=.true.)
          end if
      end select
    end subroutine rotate
    !
    subroutine grid_adapt
      integer(ik) :: rmax  ! Current extent of the left/right wavefunctions
      !
      if (.not.sd_adaptive) return
      !
      !  Update radial extent
      !
      call wt_update_lrmax(verbose,wfn_r)
      if (.not.skip_left_propagation) then
        call wt_update_lrmax(verbose,wfn_l)
      end if
      !
      !  If radial grid is already at maximum, there is nothing to do.
      !
      if (sd_nradial==sd_nradial_max) return
      !
      !  If we are still within the safe range, there is nothing to do
      !
      rmax = max(wfn_l%nradial,wfn_r%nradial)
      if (rmax<safe_nr_max()) return
      !
      !  We need to resize now!
      !
      call grid_resize(rmax)
    end subroutine grid_adapt
  end subroutine propagation
  !
  !  Calculate the maximum "safe" extent of the wavefunction which current grid can handle
  !  We need to consider:
  !   a) kinetic-energy artifacts near the edge
  !   b) the absorbing boundary
  !   c) position of the iSURFV sensor sphere
  !   d) how much wavefunction could grow in a single time step
  !
  integer(ik) function safe_nr_max()
                    safe_nr_max = sd_nradial - sd_radial_edge
    if (sd_capped)  safe_nr_max = min(sd_cap_start-1,safe_nr_max)
    if (sts_volkov) safe_nr_max = min(sts_ipt_match-1,safe_nr_max)
    !
    safe_nr_max = safe_nr_max - 2*wt_adaptive_r_buffer
  end function safe_nr_max
  !
  subroutine grid_resize(rmax)
    integer(ik), intent(in) :: rmax ! Current maximum extent of the radial wavefunction
    !
    integer(ik) :: nr
    !
    if (.not.sd_adaptive_r) return
    !
    !  Make sure the new radial grid extent is sensible [see comments above safe_nr_max()]
    !
    nr = int(sd_nradial_scl * rmax,kind=ik) + 1_ik     ! Grow by a factor
    nr = nr + sd_radial_edge + 2*wt_adaptive_r_buffer  ! Add safety buffers
    if (sts_volkov .and. nr>=sts_ipt_match) then       ! We'll reach the iSURF sensor surface, switch to the full grid
      nr = sd_nradial_max
    end if
    !
    !  Make sure we don't try to increment grid size by 1 - this adds too much overhead
    !
    if (sd_nradial<sd_nradial_max) then
      nr = max(nr,int(sd_nradial_scl * sd_nradial,kind=ik))
    end if
    !
    nr = min(nr,sd_nradial_max)                        ! Make sure we do not exceed the maximum
    if (nr==sd_nradial) return  ! Grid is already of the right size, nothing to do here
    !
    if (verbose>=0) then
      write (out,"('Resizing radial wavefunction to ',i8,' grid points')") nr
      call flush_wrapper(out)
    end if
    !
    call wt_resize(wfn_l,nr)
    call wt_resize(wfn_r,nr)
    sd_nradial = nr
    call sd_initialize(repeat=.true.)
    call pt_reset_caches
  end subroutine grid_resize
  !
  subroutine choose_rotation_code
    select case (rotation_mode)
      case default
        write (out,"('spherical_tdse%choose_rotation_mode: rotation mode ',a,' is not recognized')") trim(rotation_mode)
        stop 'spherical_tdse%choose_rotation_mode - bad rotation_mode'
      case ('none','brute force','sparse')
      case ('auto')
        ! If we use a range of M values, presumably we want rotation
        if (sd_mmin/=sd_mmax) rotation_mode = 'sparse'
        ! If we know that VP is along Z, we do not need rotation
        ! This includes 'z Gaussian', 'z Sin2', etc
        if (vp_shape(1:2)=='z ') rotation_mode = 'none'
        if (vp_shape(1:3)=='zz ') rotation_mode = 'none'
        if (vp_shape(1:4)=='zzz ') rotation_mode = 'none'
        if (vp_shape(1:4)=='zero') rotation_mode = 'none'
        ! If we still do not know, we choose the general case (rotation)
        if (rotation_mode=='auto') rotation_mode = 'sparse'
        write (out,"(/'Chosen rotation_mode = ''',a,''''/)") trim(rotation_mode)
    end select 
    !
    !  Issue warnings based on the final choice of rotation_mode
    !
    select case (rotation_mode)
      case default
        stop 'spherical_tdse%choose_rotation_mode - logic error'
      case ('none')
        if (sd_mmin/=sd_mmax) then
          write (out,"(/'WARNING: Multiple, uncoupled M values are propagated.'/)") 
        end if
        if (vp_shape(1:2)/='z ' .and. vp_shape(1:3)/='zz ' .and. vp_shape(1:4)/='zzz ' .and. vp_shape(1:4)/='zero') then
          write (out,"(/'WARNING: Vector-potential may contain components along X/Y, but the propagator')")
          write (out,"( 'WARNING: assumes it is along Z. The results are likely incorrect.'/)")
        end if
      case ('brute force','sparse')
        if (sd_mmin>-sd_lmax .or. sd_mmax<sd_lmax) then
          write (out,"(/'WARNING: Not all possible projections of angular momentum are included in propagation.')")
          write (out,"( 'WARNING: The results are likely incorrect.'/)")
        end if
    end select 
    call flush_wrapper(out)
  end subroutine choose_rotation_code
  !
  subroutine math_init
    real(rk) :: dummy
    !
    if (log10(huge(1._rk))>=4930._rk) then
      ! Quad precision; don't forget "factorial_slack"
      dummy = MathFactorial(1750_ik-5_ik)
      dummy = MathLogFactorial(20000_ik)
    else if (log10(huge(1._rk))>=308._rk) then
      ! Double precision
      dummy = MathFactorial(170_ik-5_ik)
      dummy = MathLogFactorial(10000_ik)
    else
      ! Single precision
      dummy = MathFactorial(34_ik-5_ik)
      dummy = MathLogFactorial(500_ik)
    end if
  end subroutine math_init
  !
! subroutine memory_trace
!   use ISO_C_BINDING
!   interface
!     subroutine mtrace() bind(C,name="mtrace")
!     end subroutine mtrace
!   end interface
!   !
!   !  Environment variable MALLOC_TRACE must also be set, pointing to a writeable filename
!   !
!   call mtrace
! end subroutine memory_trace
  !
  subroutine start
    !$ use OMP_LIB
    logical :: have_openmp
    external :: versions
    !
!   call memory_trace
    !
    write (out,"('Version: ',a/)") __BUILD_ID__
    !
    call TimerStart('start')
    !
    call TimerStart('Initialization')
    !
    !  Some MPI libraries (OpenMPI - I am looking at you!) do not support redirection of
    !  stdin to all ranks. We will therefore try to get the name of the stdin file from
    !  the command-line arguments.
    !
    call nt_force_stdin
    !
    !  The sequence below is a little awkward. We can't complete initialization of a multi-node
    !  run until we've red in the input parameters.
    !
    read (input,nml=sph_tdse)
    !
    call nt_initialize
    !
    if (verbose>=0) then
      write (out,"(t5,a)") trim(rcsid_spherical_tdse)
      call versions
      write (out,"()")
      write (out,"('    Integer kind = ',i0,' (',i0,' decimals)')") kind(1_ik), int(log10(huge(1_ik)+1._rk))
      write (out,"('       Real kind = ',i0,' (',i0,' decimals)')") kind(1._rk), precision(1._rk)
      write (out,"('  Aux. real kind = ',i0,' (',i0,' decimals)')") kind(1._xk), precision(1._xk)
      write (out,"('Max. LAPACK kind = ',i0,' (',i0,' decimals)')") kind(1._lrk), precision(1._lrk)
      write (out,"()")
    end if
    !
    !
    write (out,"(' ===== begin simulation parameters ===== ')")
    write (out,nml=sph_tdse)
    write (out,"(' ====== end simulation parameters ====== ')")
    !
    write (out,"(/a/)") trim(comment)
    !
    have_openmp = .false.
    !$ have_openmp = .true.
    !$ if (omp_num_threads/=0) then
    !$   write (out,"('Forcing number of OpenMP threads to ',i0)") omp_num_threads
    !$   call omp_set_num_threads(omp_num_threads)
    !$ end if
    !$ write (out,"('Maximum number of OpenMP threads for this run is ',i0)") omp_get_max_threads()
    if (omp_num_threads>0 .and. .not.have_openmp) then
      write (out,"(/'WARNING: omp_num_threads input parameter has been specified, but the code was built')")
      write (out,"( 'WARNING: without OpenMP support. Execution will continue on a single CPU core.'/)")
    end if
    if (nts%n_nodes>1) then
      write (out,"('Executing on ',i0,' nodes. This is node ',i0)") nts%n_nodes, nts%this_node
    end if
    !
    if (skip_left_propagation) then
      write (out,"(/'WARNING: Left wavefunction will not be propagated.')")
      write (out,"( 'WARNING: MOST OBSERVABLES WILL BE INCORRECT.')")
      write (out,"( 'WARNING: Don''t do this unless you know exactly what you are doing!'/)")
    end if
    !
    call math_init
    !
    !  Set up potential can CAP evaluation
    !
    call pt_initialize
    call cap_initialize
    !
    !  Set up grids and operators. This is the first call, produce more output
    !
    call sd_initialize(repeat=.false.)
    !
    !  Check accuracy of numerical derivatives
    !
    if (.not.skip_tests) call derivatives_test(verbose)
    !
    !  Construct field-free states using direct diagonalization
    !
    if (.not.skip_tests) call fieldfree_test(verbose)
    !
    ! Caches must be initialized before prepare_initial_wavefunction(), 
    ! which may refer to them.
    !
    call wt_init_memory_caches(verbose)
    !
    ! Calculate transition matrix elements between bound states if asked for
    ! Only makes sense if memory caches for atomic solutions is enabled and cap disabled
    !
    call wt_transition_matrix_elements(cap_name)
    !
    !  We need a wavefunction to start from
    !
    call prepare_initial_wavefunction
    !
    !  Prepare for a t-SURF calculation
    !
    call sts_initialize_global(dt,sqrt(vp_scale**2+vp_scale_x**2))
    !  Warning: visualize_wavefunctions() depends on sts_initialize_instance() calculating
    !           the correct value of tsurf%wfn_scale, even though it is not a tSURF routine.
    call sts_initialize_instance(tsurf,wfn_l,wfn_r)
    !
    if (initial_wf_dump_prefix/=' ') then
      write (out,"(/'Dumping initial wavefunction to disk, prefix = ',a/)") trim(initial_wf_dump_prefix)
      call flush_wrapper(out)
      call nt_merge_all(wfn_l)
      call nt_merge_all(wfn_r)
      if (nts%this_node==1) then
        call dump_wavefunctions(initial_wf_dump_prefix)
      end if
    end if
    !
    call TimerStop('Initialization')
    write (out,"(/'Done with initialization'/)")
    call TimerReport
    !
    select case (task)
      case default
        write (out,"('Task ',a,' is not recognized')") trim(task)
        stop 'spherical_tdse - bad task'
      case ('real time')
        write (out,"(/'Real-time propagation'/)")
        call choose_rotation_code
        call fill_vpot_table
        call unwrap_vpot_table
        if (vp_as_is) then
          call fill_efield_table_single_sided
        else
          call fill_efield_table
        end if
        if (nts%this_node==1) then
          call preview_laser_field 
        end if
        call propagation
      case ('imaginary time')
        write (out,"(/'Imaginary-time propagation, A = ',g24.13/)") vp_scale
        call imaginary_propagation_test(verbose,wfn_l,wfn_r,apot=vp_scale,dt=dt,nstep=timesteps)
    end select
    !
    if (pt_mix_solver=='bi-CG') then
      write (out,"(/'Total number of failures in bi-CG solver = ',i0/)") bicg_failure_count
    end if
    if (dt_subdivision/='off') then
      write (out,"(/'Largest number of microsteps per time step was ',i0)") max_subdivision
      write (out,"( '    Number of time steps requiring subdivision ',i0)") timesteps_subdivided
      write (out,"( '                    Total number of microsteps ',i0)") timesteps_micro
      write (out,"( '    Average number of microsteps per time step ',f0.3/)") real(timesteps_micro,kind=rk)/timesteps
    end if
    if (sd_adaptive) then
      write (out,"()")
      if (sd_adaptive_l) write (out,"(' Largest left/right Lmax ',2i4/)") wfn_l%lmax_top, wfn_r%lmax_top
      if (sd_adaptive_l) write (out,"('   Final left/right Lmax ',2i4/)") wfn_l%lmax, wfn_r%lmax
      if (sd_adaptive_r) write (out,"('Final left/right nradial ',2i8/)") wfn_l%nradial, wfn_r%nradial
      write (out,"()")
    end if
    !
    call flush_wrapper(out)
    !
    call sts_report(tsurf,'In TDSE')
    !
    if (final_wf_dump_prefix/=' ') then
      write (out,"(/'Dumping final wavefunction to disk, prefix = ',a/)") trim(final_wf_dump_prefix)
      call flush_wrapper(out)
      call nt_merge_all(wfn_l)
      call nt_merge_all(wfn_r)
      call flush_wrapper(out)
      if (nts%this_node==1) then
        call dump_wavefunctions(final_wf_dump_prefix)
      end if
    end if
    !
    if (reconstruct_left) then
      write (out,"(/'Reconstructing left wavefunction'/)")
      call wt_reconstruct_left(verbose,wfn_l,wfn_r)
    end if
    write (out,"(/'Analyzing wavefunction composition'/)")
    call flush_wrapper(out)
    call ca_analyze(verbose,composition_threshold,composition_max_energy,wfn_l,wfn_r,tsurf)
    call sts_atend_direct(tsurf,wfn_r)
    !
    call sts_report(tsurf,'At end')
    !
    call TimerStop('start')
    call TimerReport
    !
    call nt_finalize
  end subroutine start
end module spherical_tdse
!
program main
  use spherical_tdse
  call start
end program main
