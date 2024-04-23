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
  character(len=clen), save :: rcsid_potential_tools = "$Id: potential_tools.f90,v 1.21 2023/06/09 14:10:24 ps Exp $"
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
  !                                             Zc      a1        a2          a3        a4         a5        a6        a7     a8
  real(rk), parameter   :: tong05_he(0:8) = (/ 1._rk,  1.231_rk, 0.662_rk, -1.325_rk, 1.236_rk,-0.231_rk, 0.480_rk, 0._rk, 1._rk /)
  real(rk), parameter   :: tong05_ne(0:8) = (/ 1._rk,  8.069_rk, 2.148_rk, -3.570_rk, 1.986_rk, 0.931_rk, 0.602_rk, 0._rk, 1._rk /)
  real(rk), parameter   :: tong05_ar(0:8) = (/ 1._rk, 16.039_rk, 2.007_rk,-25.543_rk, 4.525_rk, 0.961_rk, 0.443_rk, 0._rk, 1._rk /)
  !
  !  Parameter tables from Schweizer, Fassbinder, and Gonzalez-Ferez, At. Data Nucl. Data Tables 72, 33-55 (1999).
  !  This potential is a special case of Tong&Lin functional form, with a5 through a8 being zero.
  !
  !   Parameter(Tong)    Expression (SFG99)
  !  -----------------  -------------------
  !        Zc            tilde-Z (or 1 if not specified)
  !        a1            (Z - tilde-Z)
  !        a2            a1
  !        a3            a2
  !        a4            a3
  !                                             Zc      a1     a2         a3        a4        a5-a8
  real(rk), parameter   :: sfg99_li(0:8)  = (/ 1._rk,  2._rk, 3.395_rk,  3.212_rk, 3.207_rk, 0._rk, 0._rk, 0._rk, 0._rk /)
  real(rk), parameter   :: sfg99_na(0:8)  = (/ 1._rk, 10._rk, 7.902_rk, 23.510_rk, 2.688_rk, 0._rk, 0._rk, 0._rk, 0._rk /)
  real(rk), parameter   :: sfg99_k (0:8)  = (/ 1._rk, 18._rk, 3.491_rk, 10.591_rk, 1.730_rk, 0._rk, 0._rk, 0._rk, 0._rk /)
  real(rk), parameter   :: sfg99_rb(0:8)  = (/ 1._rk, 36._rk, 3.431_rk, 10.098_rk, 1.611_rk, 0._rk, 0._rk, 0._rk, 0._rk /)
  real(rk), parameter   :: sfg99_cs(0:8)  = (/ 1._rk, 54._rk, 3.294_rk, 11.005_rk, 1.509_rk, 0._rk, 0._rk, 0._rk, 0._rk /)
  !
  !  My fit to tong05_ne, eliminating 1S core.                                                                                     
  !                                                                                                                                
  real(rk), parameter   :: neon_1S  (0:8) = (/ 1._rk,-46.1787_rk,2.46572_rk,-0.523021_rk,1.27345_rk,29.0013_rk,1.51537_rk, &
                                               0.010304_rk,0.411271_rk/)
  !
  !  My reoptimization of tong05_ne, targeting specific spin-orbit cores.
  !
  real(rk), parameter   :: neon_3_2 (0:8) = (/ 1._rk,8.0731_rk,2.0833_rk,-3.5921_rk,1.98051_rk,0.89827_rk,0.71919_rk, &
                                              -0.026825_rk,0.67424_rk /)
  real(rk), parameter   :: neon_1_2 (0:8) = (/ 1._rk,8.0854_rk,2.0374_rk,-3.5867_rk,1.92889_rk,0.87084_rk,0.79550_rk, &
                                              -0.023105_rk,0.73210_rk /)
  !
  !  My potential for the two Xenon spin-orbit cores
  !
  real(rk), parameter   :: xenon_3_2(0:8) = (/ 1._rk,5.88494_rk,0.875503_rk,0._rk,1._rk,-108.66428_rk,4.643858_rk,0._rk,1._rk /)
  real(rk), parameter   :: xenon_1_2(0:8) = (/ 1._rk,7.67004_rk,0.954979_rk,0._rk,1._rk,-108.67858_rk,4.284618_rk,0._rk,1._rk /)
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
    gjg75_pars(  2_ik, "1s2",  2.625_rk, 12.996_rk, 1.770_rk, 11.402_rk), &
    gjg75_pars(  3_ik, "2s1",  2.164_rk,  9.764_rk, 1.750_rk,  6.821_rk), &
    gjg75_pars(  4_ik, "2s2",  1.300_rk,  6.465_rk, 1.880_rk,  5.547_rk), &
    gjg75_pars(  5_ik, "2p1",  1.031_rk,  4.924_rk, 2.000_rk,  4.939_rk), &
    gjg75_pars(  6_ik, "2p2",  1.065_rk,  4.800_rk, 2.130_rk,  4.434_rk), &
    gjg75_pars(  7_ik, "2p3",  1.179_rk,  4.677_rk, 2.270_rk,  4.143_rk), &
    gjg75_pars(  8_ik, "2p4",  1.360_rk,  4.613_rk, 2.410_rk,  3.925_rk), &
    gjg75_pars(  9_ik, "2p5",  1.508_rk,  4.602_rk, 2.590_rk,  3.755_rk), &
    gjg75_pars( 10_ik, "2p6",  1.792_rk,  4.515_rk, 2.710_rk,  3.671_rk), &
    gjg75_pars( 11_ik, "3s1",  1.712_rk,  3.923_rk, 2.850_rk,  3.469_rk), &
    gjg75_pars( 12_ik, "3s2",  1.492_rk,  3.452_rk, 3.010_rk,  3.269_rk), &
    gjg75_pars( 13_ik, "3p1",  1.170_rk,  3.191_rk, 3.170_rk,  3.087_rk), &
    gjg75_pars( 14_ik, "3p2",  1.012_rk,  2.933_rk, 3.260_rk,  2.958_rk), &
    gjg75_pars( 15_ik, "3p3",  0.954_rk,  2.659_rk, 3.330_rk,  2.857_rk), &
    gjg75_pars( 16_ik, "3p4",  0.926_rk,  2.478_rk, 3.392_rk,  2.739_rk), &
    gjg75_pars( 17_ik, "3p5",  0.933_rk,  2.368_rk, 3.447_rk,  2.633_rk), &
    gjg75_pars( 18_ik, "3p6",  0.957_rk,  2.165_rk, 3.500_rk,  2.560_rk), &
    gjg75_pars( 19_ik, "3d1",  0.964_rk,  2.151_rk, 3.516_rk,  2.509_rk), &
    gjg75_pars( 20_ik, "3d2",  0.941_rk,  2.248_rk, 3.570_rk,  2.404_rk), &
    gjg75_pars( 21_ik, "3d3",  0.950_rk,  2.324_rk, 3.627_rk,  2.328_rk), &
    gjg75_pars( 22_ik, "3d4",  0.998_rk,  2.345_rk, 3.667_rk,  2.238_rk), &
    gjg75_pars( 23_ik, "3d5",  1.061_rk,  2.243_rk, 3.709_rk,  2.171_rk), &
    gjg75_pars( 24_ik, "3d6",  1.138_rk,  2.291_rk, 3.745_rk,  2.187_rk), &
    gjg75_pars( 25_ik, "3d7",  1.207_rk,  2.408_rk, 3.803_rk,  2.090_rk), &
    gjg75_pars( 26_ik, "3d8",  1.308_rk,  2.391_rk, 3.840_rk,  2.088_rk), &
    gjg75_pars( 27_ik, "3d9",  1.397_rk,  2.462_rk, 3.891_rk,  2.048_rk), &
    gjg75_pars( 28_ik, "3d10", 1.455_rk,  2.397_rk, 3.973_rk,  1.925_rk), &
    gjg75_pars( 29_ik, "4s1",  1.520_rk,  2.246_rk, 4.000_rk,  1.985_rk), &
    gjg75_pars( 30_ik, "4s2",  1.538_rk,  2.106_rk, 4.050_rk,  1.378_rk), &
    gjg75_pars( 31_ik, "4p1",  1.541_rk,  1.988_rk, 4.110_rk,  2.001_rk), &
    gjg75_pars( 32_ik, "4p2",  1.512_rk,  1.914_rk, 4.182_rk,  1.897_rk), &
    gjg75_pars( 33_ik, "4p3",  1.492_rk,  1.990_rk, 4.230_rk,  1.782_rk), &
    gjg75_pars( 34_ik, "4p4",  1.460_rk,  1.857_rk, 4.290_rk,  1.772_rk), &
    gjg75_pars( 35_ik, "4p5",  1.407_rk,  1.897_rk, 4.369_rk,  1.686_rk), &
    gjg75_pars( 36_ik, "4p6",  1.351_rk,  1.872_rk, 4.418_rk,  1.611_rk), &
    gjg75_pars( 37_ik, "4d1",  1.286_rk,  1.686_rk, 4.494_rk,  1.619_rk), &
    gjg75_pars( 39_ik, "4d3",  1.129_rk,  1.784_rk, 4.618_rk,  1.509_rk), &
    gjg75_pars( 41_ik, "4d5",  1.139_rk,  1.702_rk, 4.680_rk,  1.485_rk), &
    gjg75_pars( 42_ik, "4d6",  1.136_rk,  1.694_rk, 4.749_rk,  1.412_rk), &
    gjg75_pars( 44_ik, "4d8",  1.197_rk,  1.601_rk, 4.769_rk,  1.435_rk), &
    gjg75_pars( 46_ik, "4d10", 1.246_rk,  1.587_rk, 4.829_rk,  1.397_rk), &
    gjg75_pars( 48_ik, "5s2",  1.205_rk,  1.358_rk, 4.904_rk,  1.414_rk), &
    gjg75_pars( 50_ik, "5p2",  1.130_rk,  1.395_rk, 4.990_rk,  1.324_rk), &
    gjg75_pars( 52_ik, "5p4",  1.050_rk,  1.354_rk, 5.050_rk,  1.314_rk), &
    gjg75_pars( 54_ik, "5p6",  1.044_rk,  1.107_rk, 5.101_rk,  1.316_rk)  &
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
      if (c==0._rk) cycle terms
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
        call select_gjg(nint(pot_param(1),kind=ik),nint(pot_param(2),kind=ik))
        call report_gjg
      case ('[SFG99] Li')
        write (out,"('""Lithium"" from W Schweizer, P. Fassbinder, and R. Gonzalez-Ferez, " // &
                   "At. Data Nucl. Data Tables, 72, 33-55 (1999)')")
      case ('[SFG99] Na')
        write (out,"('""Sodium"" from W Schweizer, P. Fassbinder, and R. Gonzalez-Ferez, " // &
                   "At. Data Nucl. Data Tables, 72, 33-55 (1999)')")
      case ('[SFG99] K')
        write (out,"('""Potassium"" from W Schweizer, P. Fassbinder, and R. Gonzalez-Ferez, " // &
                   "At. Data Nucl. Data Tables, 72, 33-55 (1999)')")
      case ('[SFG99] Rb')
        write (out,"('""Rubidium"" from W Schweizer, P. Fassbinder, and R. Gonzalez-Ferez, " // &
                   "At. Data Nucl. Data Tables, 72, 33-55 (1999)')")
      case ('[SFG99] Cs')
        write (out,"('""Caesium"" from W Schweizer, P. Fassbinder, and R. Gonzalez-Ferez, " // &
                   "At. Data Nucl. Data Tables, 72, 33-55 (1999)')")
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
      case ('[SFG99] Li')
        v = tong05(r,sfg99_li)
      case ('[SFG99] Na')
        v = tong05(r,sfg99_na)
      case ('[SFG99] K')
        v = tong05(r,sfg99_k)
      case ('[SFG99] Rb')
        v = tong05(r,sfg99_rb)
      case ('[SFG99] Cs')
        v = tong05(r,sfg99_cs)
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
      v = v + (1._rk/(2._rk*electron_mass)) * real(l,kind=rk)*real(l+1,kind=rk) / r**2
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
