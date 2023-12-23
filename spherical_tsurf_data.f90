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
!  Public data structures for "spherical_tsurf.f90"
!
module spherical_tsurf_data
  use accuracy
  use constants
  implicit none
  private
  ! Exposed global data
  public sts_active, sts_active_atend, sts_active_direct
  public sts_volkov, sts_volkov_atend, sts_coulomb_atend
  public sts_step_frequency
  public sts_verbose
  public sts_rmatch, sts_ipt_match, sts_kgrid, sts_dgrid
  public sts_kgrid_min, sts_kgrid_max, sts_kgrid_count
  public sts_kgrid_file
  public sts_dgrid_ntheta, sts_dgrid_nphi, sts_dgrid_count
  public sts_dgrid_file
  public sts_asymptotic_q
  public sts_volkov_table, sts_volkov_opendx
  public sts_coulomb_waves, sts_coulomb_table, sts_coulomb_opendx
  public sts_atend_block, sts_r2r_scale, sts_atend_mode
  public sts_data
  public rcsid_spherical_tsurf_data
  !
  character(len=clen), save :: rcsid_spherical_tsurf_data = "$Id: spherical_tsurf_data.f90,v 1.2 2021/04/26 15:44:44 ps Exp $"
  !
  !  Common parameters for all t-SURF instances
  !
  logical, save               :: sts_active        = .true.       ! Set to .true. to activate during the TDSE run
  logical, save               :: sts_active_atend  = .true.       ! Set to .true. to activate during composition analysis 
                                                                  ! The two sts_active* flags are derived from the sts_volkov,
                                                                  ! and sts_volkov_atend, and sts_coulomb_atend input parameters
  logical, save               :: sts_active_direct = .true.       ! Set to .true. to activiate direct algorithm. Mutually
                                                                  ! exclusive with sts_active_atend.
  logical, save               :: sts_volkov        = .false.      ! Accumulate Volkov amplitudes during the TDSE simulation
  logical, save               :: sts_volkov_atend  = .false.      ! Accumulate residual Volkov amplitudes during composition analysis
  logical, save               :: sts_coulomb_atend = .false.      ! Calculate amplitudes of the Coulomb sherical waves during composition
                                                                  ! analysis.
                                                                  ! sts_volkov, sts_volkov_atend, and sts_coulomb_atend can be set 
                                                                  ! independently; all possible combinations of these flags may be 
                                                                  ! meaningful depending on the intent and parameters of the simulation.
  character(len=clen), save   :: sts_atend_mode    = 'direct'     ! Algorithm used for the infinite-time correction. Can be one of:
                                                                  !   'spectral' - Use spectral representation of the atomic Hamiltonian
                                                                  !   'direct'   - Use the direct expression
                                                                  ! 'direct' has much lower resource requirements, and is the default.
  integer(ik), save           :: sts_verbose       = 1            ! Verbosity level.
  integer(ik), save           :: sts_step_frequency= 8            ! Perform a t-SURF timestep at a each sts_step_frequency-th propagation
                                                                  ! time step. This helps to reduce the overhead of the t-SURF calculation,
                                                                  ! at the cost of reducing the accuracy of the spectrum. High-frequency
                                                                  ! componentns are mostly affected by the step frequency reduction.
                                                                  ! sts_step_frequency only applies if sts_volkov=.true.
  real(rk), save              :: sts_rmatch        = -10._rk      ! Matching sphere radius. Negative values will be interpreted as
                                                                  ! the distance before the start of the absorbing boundary. We'll try to
                                                                  ! find the closest matching point in the radial R grid (sd_rtab)
  real(rk), save              :: sts_asymptotic_q  = 1._rk        ! Asymptotic nuclear charge for the coulomb-wave projection; only
                                                                  ! relevant if sts_coulomb_atend is true.
  integer(ik), save           :: sts_ipt_match                    ! Index of the matching point in sd_rtab.
  character(len=clen), save   :: sts_kgrid         =  'uniform'   ! Radial grid for k. Allowed values are:
                                                                  !  'uniform' = (sts_kgrid_min:sts_kgrid_max] grid with sts_kgrid_count points
                                                                  !  'read'    = Read the desired K magnitudes from an external file
  character(len=clen), save   :: sts_kgrid_file    =  'km.table'  ! Name of the file containing K magnitudes. sts_kgrid_count must be specified.
  character(len=clen), save   :: sts_dgrid         =  'product'   ! Angular grid for k. Allowed values are:
                                                                  !  'product' = direct product of uniform grids along theta and phi
                                                                  !  'read'    = Read the desired K directions fron an external file
  character(len=clen), save   :: sts_dgrid_file    = 'kd.table'   ! Name of the file containing K directions. sts_dgrid_count must be specified.
  real(rk), save              :: sts_kgrid_min     =  0.0_rk      ! Minimum k value; sts_kgrid_min is NOT a part of the grid
  real(rk), save              :: sts_kgrid_max     =  3.0_rk      ! Maximum k value; sts_kgrid_max is part of the grid
  integer(ik), save           :: sts_kgrid_count   =  100_ik      ! Number of k values in the radial k grid
  integer(ik), save           :: sts_dgrid_ntheta  =   15_ik      ! Number of points in the theta grid
  integer(ik), save           :: sts_dgrid_nphi    =   30_ik      ! Number of points in the phi grid
  integer(ik), save           :: sts_dgrid_count                  ! Total number of points in the angular grid
  character(len=clen), save   :: sts_volkov_opendx = 'tsurf.dx'   ! Output file for OpenDX visualization; Use blank (' ') to disable.
  character(len=clen), save   :: sts_volkov_table  = ' '          ! Output file for Volkov state amplitudes; Use blank (' ') for standard out
  character(len=clen), save   :: sts_coulomb_waves = ' '          ! Output file for Coulomb partial-wave amplitudes; Use blank (' ') for standard out
  character(len=clen), save   :: sts_coulomb_table = ' '          ! Output file for Coulomb photoelectron spectrum; Use blank (' ') for standard out
  character(len=clen), save   :: sts_coulomb_opendx= 'coulomb.dx' ! Output file for OpenDX visualization; Use blank (' ') to disable.
  integer(ik), save           :: sts_atend_block   = 32_ik        ! Maximum number of solutions to process as a single block in 
                                                                  ! fill_stationary_solutions_at_matching_point()
  real(rk), save              :: sts_r2r_scale     = -1._rk       ! Scale factor to go from the right wavefunction to the real-space wavefunction
                                                                  ! If negative, the scale will be determined from the initial wavefunction at the
                                                                  ! start of the simulation. Please do not adjust unless you know EXACTLY what
                                                                  ! you are doing, and why.
  !
  !  t-SURF state
  !
  type sts_data
    real(rk)                 :: wfn_scale      ! Scale factor for converting the right wavefunction into the
                                               ! real-space wavefunction in the vicinity of the dividing surface.
                                               ! The factor is wavefunction-specific; it will differ between 
                                               ! different time steps and for different stationary states
    !
    !  Stuff needed for dealing with the Volkov solutions
    !
    real(rk), allocatable    :: vphase   (:,:) ! Volkov phase for each k at the current time, modulo 2*pi
                                               ! First index:  1:sts_kgrid_count (radial part of the k grid)
                                               ! Second index: 1:sts_dgrid_count (angular part of the k grid)
                                               ! For distributed-memory parallelization, vphase is independently
                                               ! maintained on each node
    complex(rk), allocatable :: pref     (:,:) ! Volkov pre-factor; see update_volkov_phases() and
                                               ! sts_atend_prepare(). The two routines define the %pref
                                               ! in a different way - so beware!
                                               ! Same indices as vphase.
                                               ! For distributed-memory parallelization, pref is independently
                                               ! maintained on each node
    complex(rk), allocatable :: amplitude(:,:) ! Amplitude of the Volkov state for each k at the current time
                                               ! Same indices as vphase
                                               ! For distributed-memory parallelization, amplitude is accumulated
                                               ! independently on each node, and must be added up together to
                                               ! get the final result.
    !
    !  Data for the Coulomb function evaluation
    !
    real(rk), allocatable    :: fcgr(:,:,:)    ! Coulomb functions and their gradients
                                               ! First index: 1 = function; 2 = gradient. 
                                               ! Second index: k magnitude
                                               ! Third index: L value
    complex(rk), allocatable :: fcamp(:,:,:)   ! Amplitudes of the Coulomb functions
                                               ! First index: k magnitude
                                               ! Second index: angular momentum L
                                               ! Last index: angular momentum projection M
    complex(rk), allocatable :: coulamp(:,:)   ! Amplitude of Coulomb planewave for each k at the current time
                                               ! Same indices as vphase
  end type sts_data
end module spherical_tsurf_data
