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
!  Library of simple atomic potentials
!  These routines should not be used within a critical path; if their output is
!  needed repeatedly, then cache it.
!
module potential_tools
  use accuracy
  use constants
  use timer
  implicit none
  private
  public pot_name, pot_param
  public pot_lmax, pot_tong05, pot_input
  public pot_mask, pot_mask_r0, pot_mask_rx, pot_shift
  public pt_initialize
  public pt_evaluate_potential
  public rcsid_potential_tools
  !
  character(len=clen), save :: rcsid_potential_tools = "$Id: potential_tools.f90,v 1.19 2021/04/26 15:44:44 ps Exp ps $"
  !
  !  We recognize
  !
  character(len=20), save   :: pot_name     = 'hydrogenic'   ! Name of the potential; can be one of:
                                                             ! 'hydrogenic'    - Z/r
                                                             ! 'yukawa'        - Z*exp(-A*r)/r
                                                             ! 'Gauss-Coulomb' - Z*exp(-A*r**2)/r
                                                             ! 'harmonic'      - (1/2)*k*r**2
                                                             ! 'argon'         - Effective potential from HGM PRA 60, 1341 (1999)
                                                             ! 'argon 2P'      - Valence-only effective potential. This is NOT the
                                                             !                   same potential given in HGM; I can't reproduce
                                                             !                   any results claimed for that potential.
                                                             ! 'argon 3S'      - Valence-only effective potential with [3S] core.
                                                             !                   Non-local, handle with care!
                                                             ! Effective potentials from Tong and Lin, J Phys B 38, 2593 (2005).
                                                             ! These are a sum of Coulombic, Yukawa, and exponential terms.
                                                             ! '[Tong05] He'   - "Helium"
                                                             ! '[Tong05] Ne'   - "Neon"
                                                             ! '[Tong05] Ar'   - "Argon"
                                                             ! Effective potentials from Garvey, Jackman, and Green, PRA 12, 1144 (1975).
                                                             ! '[GJC75]'       - This potential requires two input parameters: the charge
                                                             !                   of the nucleus [in pot_param(1)] and the number of electrons
                                                             !                   [in pot_param(2)]. Not all possible combinations are
                                                             !                   meaningful; see the paper for details. It is a good idea
                                                             !                   to set sd_rgrid_zeta to pot_param(1) as well!
                                                             ! Miscellaneous potentials
                                                             ! 'xenon 3/2'     - Valence-only effective potential, J=3/2 core
                                                             ! 'xenon 1/2'     - Valence-only effective potential, J=1/2 core
                                                             ! 'neon 3/2'      - All-electron neon, J=3/2 2p5 core
                                                             ! 'neon 1/2'      - All-electron neon, J=3/2 2p5 core
                                                             ! 'input Tong05'  - Channel-depended effective potential in the
                                                             !                   Tong05 form is supplied in pot_lmax and pot_tong05
                                                             !                   below.
                                                             ! 'input'         - Channel-dependent effective potential in the 
                                                             !                   general-exponential form is supplied in pot_lmax
                                                             !                   and pot_input below.
                                                             ! 'input GJG75'   - Channel-dependent effective potential in the
                                                             !                   GJG75 form is supplied in pot_lmax and pot_input
  real(rk), save            :: pot_param(5) = (/1._rk,0._rk,0._rk,0._rk,0._rk/)
                                                             ! Parameters of the potential; meaning depends on the value of pot_name
  character(len=20), save   :: pot_mask     = 'none'         ! Masking function to apply to the potential; useful for t-SURF,
                                                             ! especially with boundaries close to the origin. Can be:
                                                             ! 'none'       - Use the potential as is
                                                             ! 'becke-3'    - Third-order Becke switching function
  real(rk), save            :: pot_mask_r0  = 1000._rk       ! Starting point for the switch. For r<=pot_mask_r0, the potential
                                                             ! is not modified.
  real(rk), save            :: pot_mask_rx  = 2000._rk       ! Final point for the switch. For r>=pot_mask_rx, the potential is
                                                             ! forced to zero.
  real(rk), save            :: pot_shift    = 0._rk          ! Overall shift in the potetential BEFORE masking function is applied.0
  integer(ik), save         :: pot_lmax     = 0              ! Highest channel L in pot_tong05; all higher L's should use pot_lmax.
                                                             ! Setting pot_lmax=0 requests multiplicative effective potential.
  real(rk), save            :: pot_input (0:15,0:4) = 0._rk  ! Parameters of the general-exponential or general-GJG form. 
                                                             ! The second index is the channel L (0 for a multiplicative potential)
                                                             ! In the general-exponential form (pot_name=='input'), the potential is:
                                                             !    The first index is:
                                                             !    0     = Z
                                                             !    3*n+1 = N; 3*n+2 = C; 3*n+3 = A
                                                             !    u = -Z/r + Sum_{n=0} C * r**N * exp(-A*r)
                                                             ! In the general-GJG form (pot_name=='input GJG'), the potential is:
                                                             !    The first index is:
                                                             !    0     = Z
                                                             !    3*n+1 = C; 3*n+2 = xi; 3*n+3 = eta
                                                             !    u = (1/r) * [ -Z + Sum_{n=0} C/( 1+(eta/xi)*(exp(xi*r)-1) ) ]
  real(rk), save            :: pot_tong05(0:8,0:4)  = 0._rk  ! Parameters of the extended Tong-Lin form.
  !
  !  Parameter tables from Tong and Lin, J Phys B 38, 2593 (2005). 
  !
  !  The potential is in the form: -(1/r)*(Z + a1*exp(-a2*r) + a3*r*exp(-a4*r) + a5*exp(-a6*r) + a7*r**2*exp(-a8*r))
  !
  !  The last term is not actually in Tong&Lin.
  !
  !                                             Zc     a1      a2        a3        a4        a5      a6       a7        a8
  real(rk), parameter   :: tong05_he(0:8) = (/ 1.0,   1.231,  0.662,   -1.325,    1.236,   -0.231,  0.480,   0.000,    1.000   /)
  real(rk), parameter   :: tong05_ne(0:8) = (/ 1.0,   8.069,  2.148,   -3.570,    1.986,    0.931,  0.602,   0.000,    1.000   /)
  real(rk), parameter   :: tong05_ar(0:8) = (/ 1.0,  16.039,  2.007,  -25.543,    4.525,    0.961,  0.443,   0.000,    1.000   /)
  !
  !  My fit to tong05_ne, eliminating 1S core.                                                                                     
  !                                                                                                                                
  real(rk), parameter   :: neon_1S  (0:8) = (/ 1.0, -46.1787, 2.46572, -0.523021, 1.27345, 29.0013, 1.51537, 0.010304, 0.411271/)
  !
  !  My reoptimization of tong05_ne, targeting specific spin-orbit cores.
  !
  real(rk), parameter   :: neon_3_2 (0:8) = (/ 1.0, 8.0731, 2.0833, -3.5921, 1.98051, 0.89827, 0.71919, -0.026825, 0.67424 /)
  real(rk), parameter   :: neon_1_2 (0:8) = (/ 1.0, 8.0854, 2.0374, -3.5867, 1.92889, 0.87084, 0.79550, -0.023105, 0.73210 /)
  !
  !  My potential for the two Xenon spin-orbit cores
  !
  real(rk), parameter   :: xenon_3_2(0:8) = (/ 1.0,  5.88494, 0.875503, 0.0, 1.0, -108.66428, 4.643858, 0.0, 1.0 /)
  real(rk), parameter   :: xenon_1_2(0:8) = (/ 1.0,  7.67004, 0.954979, 0.0, 1.0, -108.67858, 4.284618, 0.0, 1.0 /)
  !
  !  Parameter table defining the [GJG] potentials, Table I of the paper.
  !
  type gjg75_pars
    integer(ik)      :: n                    ! Number of electrons
    character(len=5) :: conf                 ! Electronic configuration of the valence shell
    real(rk)         :: xi0, xi1, eta0, eta1 ! Numerical parameter, see fill_gjg() below
  end type gjg75_pars
  !
  type(gjg75_pars), parameter :: gjg_tab(45) = (/ &
    gjg75_pars(  2, "1s2",  2.625, 12.996, 1.770, 11.402), &
    gjg75_pars(  3, "2s1",  2.164,  9.764, 1.750,  6.821), &
    gjg75_pars(  4, "2s2",  1.300,  6.465, 1.880,  5.547), &
    gjg75_pars(  5, "2p1",  1.031,  4.924, 2.000,  4.939), &
    gjg75_pars(  6, "2p2",  1.065,  4.800, 2.130,  4.434), &
    gjg75_pars(  7, "2p3",  1.179,  4.677, 2.270,  4.143), &
    gjg75_pars(  8, "2p4",  1.360,  4.613, 2.410,  3.925), &
    gjg75_pars(  9, "2p5",  1.508,  4.602, 2.590,  3.755), &
    gjg75_pars( 10, "2p6",  1.792,  4.515, 2.710,  3.671), &
    gjg75_pars( 11, "3s1",  1.712,  3.923, 2.850,  3.469), &
    gjg75_pars( 12, "3s2",  1.492,  3.452, 3.010,  3.269), &
    gjg75_pars( 13, "3p1",  1.170,  3.191, 3.170,  3.087), &
    gjg75_pars( 14, "3p2",  1.012,  2.933, 3.260,  2.958), &
    gjg75_pars( 15, "3p3",  0.954,  2.659, 3.330,  2.857), &
    gjg75_pars( 16, "3p4",  0.926,  2.478, 3.392,  2.739), &
    gjg75_pars( 17, "3p5",  0.933,  2.368, 3.447,  2.633), &
    gjg75_pars( 18, "3p6",  0.957,  2.165, 3.500,  2.560), &
    gjg75_pars( 19, "3d1",  0.964,  2.151, 3.516,  2.509), &
    gjg75_pars( 20, "3d2",  0.941,  2.248, 3.570,  2.404), &
    gjg75_pars( 21, "3d3",  0.950,  2.324, 3.627,  2.328), &
    gjg75_pars( 22, "3d4",  0.998,  2.345, 3.667,  2.238), &
    gjg75_pars( 23, "3d5",  1.061,  2.243, 3.709,  2.171), &
    gjg75_pars( 24, "3d6",  1.138,  2.291, 3.745,  2.187), &
    gjg75_pars( 25, "3d7",  1.207,  2.408, 3.803,  2.090), &
    gjg75_pars( 26, "3d8",  1.308,  2.391, 3.840,  2.088), &
    gjg75_pars( 27, "3d9",  1.397,  2.462, 3.891,  2.048), &
    gjg75_pars( 28, "3d10", 1.455,  2.397, 3.973,  1.925), &
    gjg75_pars( 29, "4s1",  1.520,  2.246, 4.000,  1.985), &
    gjg75_pars( 30, "4s2",  1.538,  2.106, 4.050,  1.378), &
    gjg75_pars( 31, "4p1",  1.541,  1.988, 4.110,  2.001), &
    gjg75_pars( 32, "4p2",  1.512,  1.914, 4.182,  1.897), &
    gjg75_pars( 33, "4p3",  1.492,  1.990, 4.230,  1.782), &
    gjg75_pars( 34, "4p4",  1.460,  1.857, 4.290,  1.772), &
    gjg75_pars( 35, "4p5",  1.407,  1.897, 4.369,  1.686), &
    gjg75_pars( 36, "4p6",  1.351,  1.872, 4.418,  1.611), &
    gjg75_pars( 37, "4d1",  1.286,  1.686, 4.494,  1.619), &
    gjg75_pars( 39, "4d3",  1.129,  1.784, 4.618,  1.509), &
    gjg75_pars( 41, "4d5",  1.139,  1.702, 4.680,  1.485), &
    gjg75_pars( 42, "4d6",  1.136,  1.694, 4.749,  1.412), &
    gjg75_pars( 44, "4d8",  1.197,  1.601, 4.769,  1.435), &
    gjg75_pars( 46, "4d10", 1.246,  1.587, 4.829,  1.397), &
    gjg75_pars( 48, "5s2",  1.205,  1.358, 4.904,  1.414), &
    gjg75_pars( 50, "5p2",  1.130,  1.395, 4.990,  1.324), &
    gjg75_pars( 52, "5p4",  1.050,  1.354, 5.050,  1.314), &
    gjg75_pars( 54, "5p6",  1.044,  1.107, 5.101,  1.316)  &
    /)
  !
  contains
  !
  function tong05(r,a) result(v)
    real(rk), intent(in) :: r      ! Distance from the origin
    real(rk), intent(in) :: a(0:8) ! Parameter table, see tong05_* above
    real(rk)             :: v
    !
    v = -(1._rk/r)*(a(0) + a(1)*exp(-a(2)*r) + a(3)*r*exp(-a(4)*r) + a(5)*exp(-a(6)*r) + a(7)*r**2*exp(-a(8)*r))
  end function tong05
  !
  function genexp(r,a) result(v)
    real(rk), intent(in) :: r       ! Distance from the origin
    real(rk), intent(in) :: a(0:15) ! Parameter table, see pot_input() above [general-exponential form]
    real(rk)             :: v
    !
    integer(ik) :: i
    !
    v = -a(0)/r
    terms: do i=1,ubound(a,dim=1),3
      v = v + a(i+1) * r**a(i+0) * exp(-a(i+2)*r)
    end do terms
  end function genexp
  !
  function gengjg(r,a) result(v)
    real(rk), intent(in) :: r       ! Distance from the origin
    real(rk), intent(in) :: a(0:15) ! Parameter table, see pot_input() above [GJG form]
    real(rk)             :: v
    !
    integer(ik) :: i
    real(rk)    :: c, xi, eta
    !
    v = -a(0)
    terms: do i=1,ubound(a,dim=1),3
      c = a(i+0) ; xi = a(i+1) ; eta = a(i+2)
      if (c==0) cycle terms
      v = v + c / ( 1 + (eta/xi) * (exp(xi*r)-1) )
    end do terms
    v = v / r
  end function gengjg
  !
  subroutine select_gjg(z,n)
    integer(ik), intent(in) :: z, n ! Nuclear charge and the number of electrons
    !
    integer(ik) :: ipr
    real(rk)    :: xi, eta
    !
    pot_lmax  = 0
    pot_input = 0
    pot_input(0,0) = z - n + 1
    find_pars: do ipr=1,size(gjg_tab,dim=1)
      if (gjg_tab(ipr)%n/=n) cycle find_pars
      !
      write (out,"('Using GJG75 N=',i0,' C=',a,' xi0=',f0.3,' 10*xi1=',f0.3,' eta0=',f0.3,' 10*eta1=',f0.3)") &
            gjg_tab(ipr)%n, trim(gjg_tab(ipr)%conf), gjg_tab(ipr)%xi0, gjg_tab(ipr)%xi1, gjg_tab(ipr)%eta0, gjg_tab(ipr)%eta1
      !
      xi  = gjg_tab(ipr)%xi0  + 0.1_rk * gjg_tab(ipr)%xi1  * (z-n)
      eta = gjg_tab(ipr)%eta0 + 0.1_rk * gjg_tab(ipr)%eta1 * (z-n)
      pot_input(1,0) = -(n-1)
      pot_input(2,0) = xi
      pot_input(3,0) = eta
      return
    end do find_pars
    !
    write (out,"('No parameters are available for Z=',i0,' N=',i0,'. Sorry.')") z, n
    call flush_wrapper(out)
    stop 'potential_tools.f90%select_gjg - no parameters'
  end subroutine select_gjg
  !
  subroutine report_gjg
    integer(ik) :: lv

    write (out,"(/' v(r) = (1/r) * (-Z + Sum_{n=0} C_n (1+(eta_n/xi_n)*(exp(xi_n*r)-1))**-1'/)") 
    write (out,"((1x,a3,1x,4(1x,a12)))") &
          ' L ', '  Z  ', ' C_n ', ' xi_n ', ' eta_n ', &
          '---', '-----', '-----', '------', '-------'
    print_gjg: do lv=0,pot_lmax
      write (out,"(1x,i3,1x,f12.6,(t18,1x,f12.6,1x,f12.6,1x,f12.6))") lv, pot_input(:,lv)
    end do print_gjg
    write (out,"()")
  end subroutine report_gjg
  !
  subroutine pt_initialize
    integer(ik) :: lv
    !
    select case (pot_name)
      case default
        write (out,"('potential_tools%pt_setup_potential: potential ',a,' is not recognized')") trim(pot_name)
        stop 'potential_tools%pt_setup_potential - unknown potential'
      case ('hydrogenic')
        write (out,"('Coulombic nuclear potential. Nuclear charge = ',g24.13)") pot_param(1)
      case ('yukawa')
        write (out,"('Yukawa nuclear potential. Nuclear charge = ',g24.13,'. Screening range = ',g24.13)") pot_param(1:2)
      case ('Gauss-Coulomb')
        write (out,"('Gaussian-screened Coulomb. Nuclear charge = ',g24.13,'. Screening parameter = ',g24.13)") pot_param(1:2)
      case ('harmonic')
        write (out,"('Harmonic nuclear potential. Force constant = ',g24.13)") pot_param(1)
      case ('argon')
        write (out,"('""Argon"" potential from H.G. Muller, PRA 60, 1341 (1999)')")
      case ('argon 2P')
        write (out,"('""Argon"" valence-only potential, inspired by H.G. Muller, PRA 60, 1341 (1999).')")
        write (out,"('WARNING: This is NOT the same potential as HGM Eq. 2!')")
      case ('argon 3S')
        write (out,"('""Argon"" valence-only potential with [3S] core, inspired by H.G. Muller, PRA 60, 1341 (1999).')")
        write (out,"('WARNING: This is a non-local potential. Dipole acceleration and possibly other observables are wrong.')")
      case ('[Tong05] He')
        write (out,"('""Helium"" potential from X M Tong and C D Lin, J Phys B 38, 2593 (2005)')")
      case ('[Tong05] Ne')
        write (out,"('""Neon"" potential from X M Tong and C D Lin, J Phys B 38, 2593 (2005)')")
      case ('[Tong05] Ar')
        write (out,"('""Argon"" potential from X M Tong and C D Lin, J Phys B 38, 2593 (2005)')")
      case ('[GJG75]')
        write (out,"('Effective potential for Z=',i0,' N=',i0,' from Garvey, Jackman, and Green, PRA 12, 1144 (1975)')") &
               nint(pot_param(1:2))
        write (out,"('WARNING: GJG75 parameters are not designed for describing atomic excited states or ionization!'/)")
        call select_gjg(nint(pot_param(1)),nint(pot_param(2)))
        call report_gjg
      case ('neon 1S')
        write (out,"('""Neon"" valence-only potential, inspired by X M Tong and C D Lin, J Phys B 38, 2593 (2005)')")
      case ('neon 1/2')
        write (out,"('""Neon"" potential for 2p5(J=1/2) core (unpublished).')")
      case ('neon 3/2')
        write (out,"('""Neon"" potential for 2p5(J=3/2) core (unpublished).')")
      case ('xenon 3/2')
        write (out,"('""Xenon"" valence-only potential, J=3/2 ion core (unpublished)')")
      case ('xenon 1/2')
        write (out,"('""Xenon"" valence-only potential, J=1/2 ion core (unpublished)')")
      case ('input Tong05')
        write (out,"('Tong-Lin-style effective potential provided on input. lmax = ',i0)") pot_lmax
        write (out,"(/' v(r) = -Z/r - v1*exp(-a1*r)/r - v2*exp(-a2*r) - v3*exp(-a3*r)/r - v4*r**exp(-a4*r)'/)")
        if (pot_lmax<0 .or. pot_lmax>ubound(pot_tong05,dim=2)) then
          call flush_wrapper(out)
          stop 'potential_tools%pt_setup_potential - pot_lmax out of range'
        end if
        write (out,"((1x,a3,1x,9(1x,a12)))") &
              ' L ', '  Z  ', ' V1 ', ' A1 ', ' V2 ', ' A2 ', ' V3 ', ' A3 ', ' V4 ', ' A4 ', &
              '---', '-----', '----', '----', '----', '----', '----', '----', '----', '----'
        print_tl05: do lv=0,pot_lmax
          write (out,"(1x,i3,1x,9(1x,f12.6))") lv, pot_tong05(:,lv)
        end do print_tl05
      case ('input')
        write (out,"('General-exponential effective potential provided on input. lmax = ',i0)") pot_lmax
        write (out,"(/' v(r) = -Z/r + Sum_{i=0} C_i * r**N_i * exp(-A_i*r)'/)")
        if (pot_lmax<0 .or. pot_lmax>ubound(pot_input,dim=2)) then
          call flush_wrapper(out)
          stop 'potential_tools%pt_setup_potential - pot_lmax out of range'
        end if
        write (out,"((1x,a3,1x,4(1x,a12)))") &
              ' L ', '  Z  ', ' Cn ', ' Nn ', ' An ', &
              '---', '-----', '----', '----', '----'
        print_tinp: do lv=0,pot_lmax
          write (out,"(1x,i3,1x,f12.6,(t18,1x,f12.6,1x,f12.6,1x,f12.6))") lv, pot_input(:,lv)
        end do print_tinp
      case ('input GJG75')
        write (out,"('Garvey-Jackman-Green type potential provided on input')")
        call report_gjg
    end select 
    !
    select case (pot_mask)
      case default
        write (out,"('potential_tools%pt_setup_potential: mask ',a,' is not recognized')") trim(pot_mask)
        stop 'potential_tools%pt_setup_potential - unknown mask'
      case ('none')
      case ('becke-3')
        write (out,"('Potential is forced to zero at infinity using Becke order-3 step function')")
        write (out,"('       switching starts at ',g24.13,' Bohr')") pot_mask_r0
        write (out,"('  switching is complete at ',g24.13,' Bohr')") pot_mask_rx
    end select
  end subroutine pt_initialize
  !
  function pt_evaluate_potential(l,j,r,centrifugal) result(v)
    integer(ik), intent(in) :: l           ! Orbital angular momentum channel
    integer(ik), intent(in) :: j           ! Total angular momentum sub-channel (currently unused)
                                           ! -1: l-1/2 sub-channel
                                           !  0: weighed-average of the sub-channels
                                           ! +1: l+1/2 sub-channel
    real(rk), intent(in)    :: r           ! Radial position where potential is needed
    logical, intent(in)     :: centrifugal ! Add centrifugal term
    real(rk)                :: v           ! Potential at a grid point, including the centrifugal term
    !
    real(rk) :: zeff, arg
    !
    !  Sanity check: L must be non-negative; J must be in the +1,-1 range, and R must be positive
    !
    if (l<0 .or. abs(j)>1 .or. r<=0._rk) then
      stop 'potential_tools%pt_evaluate_potential - domain error'
    end if
    !
    !  First, the system-specific part of the potential
    !
    select case (pot_name)
      case default
        write (out,"('potential_tools%pt_evaluate_potential: potential ',a,' is not recognized')") trim(pot_name)
        stop 'potential_tools%pt_evaluate_potential - unknown potential'
      case ('hydrogenic')
        v = electron_charge*pot_param(1)/r
      case ('yukawa')
        v = electron_charge*pot_param(1)*exp(-pot_param(2)*r)/r
      case ('Gauss-Coulomb')
        v = electron_charge*pot_param(1)*exp(-pot_param(2)*r**2)/r
      case ('harmonic')
        v = 0.5_rk*pot_param(1)*r**2
      case ('argon')
        zeff = 1._rk + 5.4_rk*exp(-r) + (17._rk-5.4_rk)*exp(-3.682_rk*r)
        v    = electron_charge*zeff/r
      case ('argon 2P')
        zeff = 1._rk + 7.625195_rk*exp(-1.02557_rk*r)
        arg  = 10.00_rk*(r-0.37110_rk)
        if (arg<0.9_rk*log(huge(1._rk))) then
          zeff = zeff - 124.55_rk/(1._rk + exp(arg))
        end if
        v = electron_charge*zeff/r
      case ('argon 3S')
        if (l==0) then 
          zeff = 1._rk + 3.5385_rk*exp(-0.7859_rk*r) - 20._rk*exp(-1.3014_rk*r)
        else
          zeff = 1._rk + 6.5852_rk*exp(-1.0025_rk*r) - 20._rk*exp(-3.4834_rk*r)
        end if
        v = electron_charge*zeff/r
      case ('[Tong05] He')
        v = tong05(r,tong05_he)
      case ('[Tong05] Ne')
        v = tong05(r,tong05_ne)
      case ('[Tong05] Ar')
        v = tong05(r,tong05_ar)
      case ('neon 1S')
        v = tong05(r,neon_1S)
      case ('neon 3/2')
        v = tong05(r,neon_3_2)
      case ('neon 1/2')
        v = tong05(r,neon_1_2)
      case ('xenon 3/2')
        v = tong05(r,xenon_3_2)
      case ('xenon 1/2')
        v = tong05(r,xenon_1_2)
      case ('input Tong05')
        if (l>=0 .and. l<pot_lmax) then
          v = tong05(r,pot_tong05(:,l))
        else
          v = tong05(r,pot_tong05(:,pot_lmax))
        end if
      case ('input')
        if (l>=0 .and. l<pot_lmax) then
          v = genexp(r,pot_input(:,l))
        else
          v = genexp(r,pot_input(:,pot_lmax))
        end if
      case ('input GJG75','[GJG75]')
        if (l>=0 .and. l<pot_lmax) then
          v = gengjg(r,pot_input(:,l))
        else
          v = gengjg(r,pot_input(:,pot_lmax))
        end if
    end select 
    !
    !  Add overall shift
    !
    v = v + pot_shift
    !
    !  Now apply mask to the potential
    !
    select case (pot_mask)
      case default
        write (out,"('potential_tools%pt_evaluate_potential: mask ',a,' is not recognized')") trim(pot_mask)
        stop 'potential_tools%pt_evaluate_potential - unknown mask'
      case ('none')
      case ('becke-3')
        v = v * becke3(r,pot_mask_r0,pot_mask_rx)
    end select
    ! 
    !  Add the universal centrifugal term. Nominally, that's a kinetic energy contribution, but ...
    !
    if (centrifugal) then
      v = v + (1._rk/(2._rk*electron_mass)) * l*(l+1) / r**2
    end if
    !
  end function pt_evaluate_potential
  !
  real(rk) function becke3(r,r0,rx)
    real(rk), intent(in) :: r, r0, rx
    real(rk)             :: x
    !
    if (r<=r0) then
      becke3 = 1._rk
    else if (r>=rx) then
      becke3 = 0._rk
    else
      x      = 2._rk*(r-r0)/(rx-r0) - 1._rk
      becke3 = 0.5_rk - 0.5_rk * step(step(step(x)))
    end if
    !
    contains
    real(rk) function step(x)
      real(rk), intent(in) :: x
      !
      step = 1.5_rk * x - 0.5_rk * x**3
    end function step
  end function becke3
end module potential_tools
