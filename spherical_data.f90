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
!  Data structures for the H.G.M.-style spherical TDSE solver
!  Note that the actual operators are sometimes different from HGM's
!  This data may be shared by a number of modules/routines
!
!  2018 Aug 21: The initialization routines were split off to avoid circular
!               dependency on node_tools.f90
!
module spherical_data
  use accuracy
  implicit none
  private
  public sd_nradial, sd_nspin, sd_lmax, sd_mmin, sd_mmax
  public sd_adaptive, sd_adaptive_l, sd_tolerance_l, sd_adaptive_r, sd_tolerance_r
  public sd_nradial_max, sd_nradial_scl, sd_radial_edge
  public sd_rgrid, sd_rgrid_rmin, sd_rgrid_dr, sd_rgrid_r0, sd_rgrid_scale, sd_rgrid_file
  public sd_rgrid_zeta, sd_rgrid_npad
  public sd_rgrid_grad_delta
  public sd_rgrid_report
  public sd_operators
  public sd_rtab, sd_rminus, sd_drtab, sd_pottab, sd_captab, sd_dvdrtab
  public sd_rtab_full, sd_drtab_full
  public sd_pot_nonlocal
  public sd_capped, sd_cap_start
  public sd_d2n_l0, sd_m2n_l0, sd_m2nf_l0, sd_m2np_l0, sd_d2n_lx, sd_m2n_lx, sd_m2nf_lx, sd_m2np_lx
  public sd_d1n_l0, sd_m1n_l0, sd_m1nf_l0, sd_m1np_l0, sd_d1n_lx, sd_m1n_lx, sd_m1nf_lx, sd_m1np_lx
  public sd_d2t_l0, sd_m2t_l0, sd_m2tf_l0, sd_m2tp_l0, sd_d2t_lx, sd_m2t_lx, sd_m2tf_lx, sd_m2tp_lx
  public sd_d1t_l0, sd_m1t_l0, sd_m1tf_l0, sd_m1tp_l0, sd_d1t_lx, sd_m1t_lx, sd_m1tf_lx, sd_m1tp_lx
  public sd_wfn
  public rcsid_spherical_data
  !
  character(len=clen), save :: rcsid_spherical_data = "$Id: spherical_data.f90,v 1.43 2021/04/26 15:44:44 ps Exp ps $"
  !
  integer, parameter          :: iu_temp        = 23         ! An arbitrary unit number, which can be used here
  !
  !  Common data, shared by all wavefunctions which are handled by our propagator.
  !
  !  Parameters directly controlled by the user
  ! 
  integer(ik), save           :: sd_nradial     = 220_ik     ! Current number of radial points in the grid. For historical 
                                                             ! reasons, on input sd_nradial specifies the MAXIMUM number
                                                             ! or points on the radial grid. When sd_adaptive_r is .true.,
                                                             ! sd_nradial may be smaller than the maximum.
  integer(ik), save           :: sd_nradial_max = -1         ! Maximum number of radial points in the grid. It is initially
                                                             ! set from sd_nradial on the first call to sd_initialize.
  real(rk), save              :: sd_nradial_scl = 1.1_rk     ! When sd_nradial needs to be expanded, try to increase the
                                                             ! the grid to (sd_nradial_scl*sd_nradial)
  integer(ik), save           :: sd_radial_edge = -1         ! Number of points close to the outer grid edge where derivative
                                                             ! operators may be affected by resizing the grid
  integer(ik), save           :: sd_nspin       = 1_ik       ! Number of spin components in the wavefunctions; either 1 or 2
  integer(ik), save           :: sd_lmax        = 10_ik      ! Largest angular momentum in the grid (lowest is always 0)
  integer(ik), save           :: sd_mmin        =  0_ik      ! Smallest angular momentum projection 
  integer(ik), save           :: sd_mmax        =  0_ik      ! Largest angular momentum projection 
                                                             ! There are only two allowed usage cases:
                                                             ! a) sd_mmin = sd_mmax = constant (linear polarization)
                                                             ! b) sd_mmin = -sd_lmax; sd_mmax = sd_lmax (arbitrary field)
  logical, save               :: sd_adaptive                 ! Active if any adaptive-grid flags are on (sd_adaptive_l or
                                                             ! sd_adaptive_r).
  logical, save               :: sd_adaptive_l  = .true.     ! Dynamically adjust maximum angular momentum, up to sd_lmax
  real(rk), save              :: sd_tolerance_l(2) = (/ -1, -1 /)  
                                                             ! When the largest element of the angular momentum channel
                                                             ! exceeds sd_tolerance_l, increment the dynamic angular momentum.
                                                             ! Negative values request tolerance of spacing(0.1_rk)
                                                             ! See wt_update_lrmax and wt_reset_lrmax in wavefunction_tools.f90
                                                             ! The first entry applies to the right wavefunction; the second
                                                             ! to the left. If only rhe first (right) entry is set explicitly, 
                                                             ! the corresponding left entry will be chosen to maintain
                                                             ! approximately the same accuracy in the left wavefunction.
  logical, save               :: sd_adaptive_r  = .true.     ! Dynamically adjust maximum extent of the grid, up to sd_nradial
                                                             ! Only some of the terms in the propagator can benefit from the
                                                             ! finite wavefunction extent, so the benefit is quite modest.
  real(rk), save              :: sd_tolerance_r(2) = (/ -1, -1 /)         
                                                             ! When the wavefunction exceeds sd_adaptive_r, increment dynamic
                                                             ! extent. Negative values request tolerance of spacing(0.1_rk)
                                                             ! See wt_update_lrmax and wt_reset_lrmax in wavefunction_tools.f90
                                                             ! The first entry applies to the right wavefunction; the second to
                                                             ! the left. If only rhe first (right) entry is set explicitly, 
                                                             ! the corresponding left entry will be chosen to maintain
                                                             ! approximately the same accuracy in the left wavefunction.
  character(len=clen), save   :: sd_rgrid       = 'uniform'  ! Choice of the radial grid. Can be one of the following:
                                                             ! 'uniform'     - Uniformly distributed points; uses sd_rgrid_dr
                                                             ! 'log'         - Logarithmically distributed points; 
                                                             !                 uses sd_rgrid_r0 and sd_rgrid_scale
                                                             ! 'log-uniform' - Logarithmic grid until step reaches sd_rgrid_dr;
                                                             !                 uniform grid afterwards
                                                             ! 'read'        - Arbitrary user-defined points; uses rd_rgrid_file
                                                             ! Note that both 'log' and 'log-uniform' add a small uniform segment
                                                             ! between the origin and the first 'log' point; otherwise we get
                                                             ! highly oscillatory terms in the derivative operators.
                                                             ! WARNING: Not all grids produced by sd_rgrid/='uniform' are
                                                             ! WARNING: of a good quality. Watch out for rapid oscillations
                                                             ! WARNING: in the left wavefunction: this can (and does!) cause
                                                             ! WARNING: problems with the laser coupling operator
  real(rk), save              :: sd_rgrid_zeta  = 1.0_rk     ! Effective nuclear charge at the origin; this value is only
                                                             ! used to establish the boundary condition on the wavefunction
  real(rk), save              :: sd_rgrid_rmin  = 0.0_rk     ! Starting point of the radial grid. The region of space for
                                                             ! for r<sd_rgrid_rmin is (implicitly) taken to be covered by an
                                                             ! infinitely-high potential barrier. Choosing a non-zero value here
                                                             ! will change the boundary conditions: instead of matching the 
                                                             ! point-charge asymptotics requested by sd_rgrid_zeta, the zero
                                                             ! boundary will apply. This is an expert-level parameter; please
                                                             ! make sure you know what's going on before changing it!
  real(rk), save              :: sd_rgrid_dr    = 0.2_rk     ! Radial grid spacing (sd_rgrid=='uniform')
  real(rk), save              :: sd_rgrid_r0    = 1e-2_rk    ! First point of the logarithmic grid (sd_rgrid=='log')
  real(rk), save              :: sd_rgrid_scale = 1.05_rk    ! Scaling factor of the logarithmic grid (sd_rgrid=='log')
                                                             ! r(i) = sd_rgrid_r0 * sd_rgrid_scale**(i-1)
  integer(ik), save           :: sd_rgrid_npad  = -1_ik      ! Number of extra "padding" uniformly-spaced grid points at the
                                                             ! origin; this only applies to 'log' and 'log-uniform' grids
                                                             ! (-1) is a special value; it will chose number of uniform padding
                                                             ! points such that grid spacing increases monotonically at the inner
                                                             ! boundary.
  character(len=clen), save   :: sd_rgrid_file  = 'grid.dat' ! File containing a list of radial grid point coordinates (in Bohr), 
                                                             ! one point per line. The initial zero point is assumed implicitly,
                                                             ! and should not be included.
                                                             ! (sd_nradial+1) points are expected (see sd_rtab(:))
                                                             ! Coordinates must be in the ascending order (sd_rgrid=='read')
  real(rk), save              :: sd_rgrid_grad_delta &
                                                = 0.01_rk    ! Step-size control for numerical differential of the potential.
                                                             ! The actual step size for point (i) is computed as:
                                                             ! rd_drtab(i)*sd_rgrid_grad_delta
  character(len=clen), save   :: sd_rgrid_report= ' '        ! Output file, containing radial grid and multiplicative potential.
                                                             ! Use blank to suppress output.
  character(len=clen), save   :: sd_operators   = 'channel'  ! Choice of differential operators (laplacian and gradient).
                                                             ! This is a debugging parameter. Do not change the default unless you know
                                                             ! what you are doing.
                                                             ! 'channel' uses channel-specific operators, differerent for L=0 and L>0
                                                             ! 'all common' uses the same operators for all channels
                                                             ! 'grad common' uses the same operators for the radial gradient, but remains
                                                             !               per-channel for the laplacian
                                                             ! 'lap common' uses common operators for the laplacian, but remains per-channel
                                                             !              for the gradient.
                                                             ! According to HGM, 'common' improves numerical stability for high fields.
                                                             ! However, it does not impose the correct boundary conditions at the origin
                                                             ! in the L=0 channel, and decreses the accuracy. 
  !
  !  Derived parameters, computed from the values above
  !
  real(rk), allocatable, save :: sd_rtab (:)                 ! Radial coordinates, (sd_nradial+1) points, in Bohr.
                                                             ! The grid starts at sd_rgrid_rmin (which is usually, but not always, zero).
                                                             ! The very first point is implicit, and is not kept in sd_rtab(). The grud
                                                             ! ends at sd_rtab(sd_nradial+1). 
  real(rk), allocatable, save :: sd_rtab_full(:)             ! Full radial grid, (sd_nradial_max+1) points, in Bohr. This grid is 
                                                             ! prepared at the first initialization call, and is then truncated as
                                                             ! needed.
  real(rk), allocatable, save :: sd_rminus(:)                ! 1/r. Goes from 1 to sd_nradial
  real(rk), allocatable, save :: sd_drtab(:)                 ! Distance between the current and preceeding grid points, in Bohr.
                                                             ! (1:sd_nradial+1)
                                                             ! sd_drtab(i) = sd_rtab(i)-sd_drtab(i-1), where sd_rtab(0) 
                                                             ! is (implicitly) taken as sd_rgrid_rmin.
  real(rk), allocatable, save :: sd_drtab_full(:)            ! Same as sd_drtab, but on the full radial grid.
  real(rk), allocatable, save :: sd_pottab(:,:,:)            ! Per-channel multiplicative potential at grid points
                                                             ! 1: (1:sd_nradial) - Radial grid indices
                                                             ! 2: (1)            - Spin-orbit sub-channel; currently unused
                                                             ! 3: (0:sd_lmax)    - Angular momentum of the channel
  logical, save               :: sd_capped                   ! .True. if sd_captab(:) is not identically zero
  integer(ik), save           :: sd_cap_start                ! First point of the capping potential
  complex(rk), allocatable, save &
                              :: sd_captab(:)                ! Complex absorbing potential at the edge of the grid.
  real(rk), allocatable, save :: sd_dvdrtab(:)               ! Radial gradient of the potential, (sd_nradial) points, in
                                                             ! Hartree per Bohr.
  logical, save               :: sd_pot_nonlocal = .false.   ! Will be set to .true. if the potential appears to be non-local
  !
  !  Differential operators, radial coordinate.
  !  All operators here are in the implicit form, namely:
  !   Op(X) = M^{-1} . D . X
  !  where M and D are three-diagonal real matrices (not necessarily symmetric!).
  !  R=0 boundary conditions are different in L=0 and all the remaining channels,
  !  so that we have to use different implicit operators.
  !
  !  The tri-diagonal matrices are stored as follows:
  !    v(:,1) - the diagonal
  !    v(:,2) - the sub-diagonal; the last element is not used
  !    v(:,3) - the super-diagonal; the last element is not used
  !
  real(rk), allocatable, save :: sd_d2n_l0 (:,:)              ! Laplacian, L=0 channel
  real(rk), allocatable, save :: sd_m2n_l0 (:,:)              ! 
  real(rk), allocatable, save :: sd_m2nf_l0(:,:)              ! sd_m2_l0, factored by m3d_decompose,
  logical,  allocatable, save :: sd_m2np_l0(  :)              ! pivot list for sd_m2nf_l0
  real(rk), allocatable, save :: sd_d2n_lx (:,:)              ! Laplacian, L>0 channel
  real(rk), allocatable, save :: sd_m2n_lx (:,:)              ! 
  real(rk), allocatable, save :: sd_m2nf_lx(:,:)              ! sd_m2_lx, factored by m3d_decompose,
  logical,  allocatable, save :: sd_m2np_lx(  :)              ! pivot list for sd_m2nf_lx
  real(rk), allocatable, save :: sd_d1n_l0 (:,:)              ! Radial gradient, L=0 channel
  real(rk), allocatable, save :: sd_m1n_l0 (:,:)              ! 
  real(rk), allocatable, save :: sd_m1nf_l0(:,:)              ! sd_m1_l0, factored by m3d_decompose,
  logical,  allocatable, save :: sd_m1np_l0(  :)              ! pivot list for sd_m1nf_l0
  real(rk), allocatable, save :: sd_d1n_lx (:,:)              ! Radial gradient, L>0 channel
  real(rk), allocatable, save :: sd_m1n_lx (:,:)              ! 
  real(rk), allocatable, save :: sd_m1nf_lx(:,:)              ! sd_m1_lx, factored by m3d_decompose,
  logical,  allocatable, save :: sd_m1np_lx(  :)              ! pivot list for sd_m1nf_lx
  !
  !  Transposes
  !
  real(rk), allocatable, save :: sd_d2t_l0 (:,:)              ! Laplacian, L=0 channel
  real(rk), allocatable, save :: sd_m2t_l0 (:,:)              ! 
  real(rk), allocatable, save :: sd_m2tf_l0(:,:)              ! sd_m2_l0, factored by m3d_decompose,
  logical,  allocatable, save :: sd_m2tp_l0(  :)              ! pivot list for sd_m2tf_l0
  real(rk), allocatable, save :: sd_d2t_lx (:,:)              ! Laplacian, L>0 channel
  real(rk), allocatable, save :: sd_m2t_lx (:,:)              ! 
  real(rk), allocatable, save :: sd_m2tf_lx(:,:)              ! sd_m2_lx, factored by m3d_decompose,
  logical,  allocatable, save :: sd_m2tp_lx(  :)              ! pivot list for sd_m2tf_lx
  real(rk), allocatable, save :: sd_d1t_l0 (:,:)              ! Radial gradient, L=0 channel
  real(rk), allocatable, save :: sd_m1t_l0 (:,:)              ! 
  real(rk), allocatable, save :: sd_m1tf_l0(:,:)              ! sd_m1_l0, factored by m3d_decompose,
  logical,  allocatable, save :: sd_m1tp_l0(  :)              ! pivot list for sd_m1tf_l0
  real(rk), allocatable, save :: sd_d1t_lx (:,:)              ! Radial gradient, L>0 channel
  real(rk), allocatable, save :: sd_m1t_lx (:,:)              ! 
  real(rk), allocatable, save :: sd_m1tf_lx(:,:)              ! sd_m1_lx, factored by m3d_decompose,
  logical,  allocatable, save :: sd_m1tp_lx(  :)              ! pivot list for sd_m1tf_lx
  !
  !  1-electron wavefunction; there may be more than one floating around;
  !
  type sd_wfn          
    !
    !  IMPORTANT: If you make any changes to this data structure, please do not
    !  IMPORTANT: forget to change checkpoint_tools and node_tools, in addition
    !  IMPORTANT: to the obvious places!
    !
    complex(rk), allocatable :: wfn(:,:,:,:)     ! The indices are:
                                                 ! 1 - [1:sd_nradial]    - Radial grid points
                                                 ! 2 - [1:sd_nspin]      - Spin components
                                                 ! 3 - [0:sd_lmax]       - Angular momentum components
                                                 !                         Note that only the range [abs(m):sd_lmax] is used
                                                 ! 4 - [sd_mmin:sd_mmax] - Projection of the angular momentum on the local Z axis
    integer(ik)              :: lmax             ! Largest non-zero angular momentum channel. Must not exceed sd_lmax
    integer(ik)              :: lmax_top         ! The largest value lmax reached during the simulation.
    integer(ik)              :: nradial          ! Maximum extent of the wavefunction. Must not exceed sd_nradial
                                                 ! The next two entries are managed by the input routines; they are left alone
                                                 ! by the checkpoint tools.
    real(rk)                 :: sd_tol_r         ! The appropriate entry of sd_tolerance_r
    real(rk)                 :: sd_tol_l         ! The appropriate entry of sd_tolerance_l
  end type sd_wfn
  !
end module spherical_data
