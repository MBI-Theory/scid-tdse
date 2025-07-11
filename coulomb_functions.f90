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
!  Calculate Coulomb functions using the algorithm of Barnett:
!
!   A.R. Barnett, Comp Phys Comm 21, 297-314 (1981); 27, 147-166 (1982).
!
!  Code is this module is not intended to be efficient, but must 
!  be numerically accurate for all inputs and real kinds.
!
!  Two functions are provided:
!
!   coulombFG calculates regular and irregular solutions in scaled coordinates, 
!             following Barnett's paper. The coordinates are:
!              x = k r; eta = -Z/k
!             Many other special functions can be evaluated in terms of these
!             solutions; see Barnett's paper for details.
!
!             WARNING: Numerical overflow may occur for some input parameters;
!             WARNING: Irregular solutions for large L values are especially
!             WARNING: to overflow
!
!   coulombF  calculates regular hydrogenic solutions. A number of different
!             forms of this function is present in the literature. Specifically:
!             a) The version in Burgess-1965 paper differs by an extra factor of Z**0.5
!             b) The version in Flugge differs by an extra factor of 0.5
!             c) The version in L&L v3 differs by an extra factor of (1/r)
!
!             Our version is normalize to delta((k-k')/(2 Pi))
!
!   coulombF may lose significance if:
!             k*r is very large; in this case, we will lose about log10(k*r) digits
!             k*r is small, Z/k is large, and lambdamin is large; in this case, we
!                 may lose as many digits as are in the ratio of G/F for lambdamin
!                 A warning will appear of the expected accuracy loss exceeds 3 digits
!
!   The analytical form of coulombF is:
!
!     Exp(pi Z/(2 k)) Abs(Gamma(L+1-I Z/k)) ((2 k r)^(L+1)/(2L+1)!) 
!             Exp(-I k r) 1F1(I Z/k + L +1, 2L+2, 2 I k r)
!
!   with the asymptotic form at large R being:
!
!     2 Sin( k r + (Z/k) Log(2 k r) - L pi/2 + Arg(Gamma(L+1 - I Z/k)))
!
!  Additionally, radial part of the bound hydrogenic functions are evaluated 
!  using a simple recursive expression for the generalized Laguerre polynomials. 
!  In order to obtail the functions in R^3, normalized to 1, the radial part
!  should be multiplied by the appropriate spherical harmonics.
!
 module coulomb_functions
   use accuracy
   use math
   implicit none
   private
   public ik, rk
   public coulombFG, coulombF, coulombBound
   public rcsid_coulomb_functions
   !
   character(len=clen), save :: rcsid_coulomb_functions = "$Id: coulomb_functions.f90,v 1.16 2024/04/23 14:30:15 ps Exp $"
   !
 ! integer, parameter      :: out         = 6
 ! integer, parameter      :: ik          = selected_int_kind(15)
 ! integer, parameter      :: rk          = selected_real_kind(28,50)
 ! integer, parameter      :: ik          = selected_int_kind(6)
 ! integer, parameter      :: rk          = selected_real_kind(14,50)
   integer, parameter      :: cf1_maxiter = 10000000_ik
   integer, parameter      :: cf2_maxiter = 10000000_ik
   complex(rk), parameter  :: ii          = (0._rk,1._rk)
   !
   contains
   !
   subroutine coulombF(lmin,lmax,Z,k,r,f,fp)
     integer(ik), intent(in) :: lmin          ! Minimum value of L
     integer(ik), intent(in) :: lmax          ! Maximum value of L
     real(rk), intent(in)    :: Z             ! Nuclear charge; Z>0 means attraction
     real(rk), intent(in)    :: k             ! Asymptotic wavevector (k>0).
     real(rk), intent(in)    :: r             ! Distance where the solution is needed (r>0).
     real(rk), intent(out)   :: f (lmin:lmax) ! Coulomb function for the specified range of L
     real(rk), intent(out)   :: fp(lmin:lmax) ! Radial derivative of f
     !
     real(rk)    :: lambdamin, x, eta
     real(rk)    :: fg(lmin:lmax,4)
     integer(ik) :: nlambda
     !
     !  Sanity checking ...
     !
     if (lmin<0)    stop 'coulomb_functions%coulombF - lmin must be non-negative'
     if (lmax<lmin) stop 'coulomb_functions%coulombF - lmax must be >= lmin'
     if (k<=0)      stop 'coulomb_functions%coulombF - k must be positive'
     if (r<=0)      stop 'coulomb_functions%coulombF - r must be positive'
     !
     !  Calculate Barnett's Coulomb functions 
     !
     lambdamin = real(lmin,kind=rk)
     nlambda   = lmax-lmin+1
     x         = k*r
     eta       = -Z/k
     call coulombFG(lambdamin,nlambda,x,eta,fg,skip_g=.true.)
     !
     !  Convert to the hydrogenic Coulomb functions
     !
     f (:) = 2.0_rk*  fg(:,1)
     fp(:) = 2.0_rk*k*fg(:,2)
   end subroutine coulombF
   !
   subroutine coulombFG(lambdamin,nlambda,x,eta,fg,skip_g)
     real(rk), intent(in)          :: lambdamin ! Minimum value of L
     integer(ik), intent(in)       :: nlambda   ! Number of values of L to calculate
     real(rk), intent(in)          :: x         ! Coordinate x = k*r
     real(rk), intent(in)          :: eta       ! Charge parameter, eta = -Z/k
     real(rk), intent(out)         :: fg(:,:)   ! Table of coulomb functions and their derivatives with respect to x
                                                ! First index:  i [1 .. nlambda]: L = lambdamin + i-1
                                                ! Second index: 1 = F; 2 = F'; 3 = G; 4 = G'
                                                ! F is the regular solution (finite at zero and sine asymptotics)
                                                ! G is the irregular solution (divergent at zero, cosine asymptotics)
     logical, intent(in), optional :: skip_g    ! Set to .true. to skip calculation of the irregular solution
                                                ! Note that space in fg() reserved for G may still be modified!
     !
     integer(ik) :: indl                           ! lambda=lambdamin+indl-1
     real(rk)    :: lambda, lambdamax              ! lambdamax=lambdamin+nlambda-1
     real(rk)    :: rfmax, sfmax                   ! F'/F and sign(F) for lambdamax
     complex(rk) :: pq                             ! (F+I*G)'/(F+I*G)
     real(rk)    :: p, q                           ! Re(pq) and Im(pq)
     real(rk)    :: wrs, om1                       ! Wronskian and normalization factor
     real(rk)    :: rt(2:nlambda), st(2:nlambda)
     !
     if (x<=0) stop 'coulomb_functions%coulombFG - bad x'
     if (nlambda<=0) stop 'coulomb_functions%coulombFG - bad nlambda'
     if (ubound(fg,1)/=nlambda .or. ubound(fg,2)/=4_ik) stop 'coulomb_functions%coulombFG - bad array sizes'
     !
     !  Start by filling coefficient tables: we'll need these repeatedly.
     !
     fill_rs: do indl=2,nlambda
       lambda   = lambdamin + indl - 1
       rt(indl) = brnR(lambda,eta)
       st(indl) = brnS(lambda,x,eta)
     end do fill_rs
     !
     !  Step 1: figure out F'/F and the sign of F for lambda = lambdamax
     !
     lambdamax = lambdamin + nlambda - 1
     rfmax     = cf1(lambdamax,x,eta,sfmax)
     !
     !  Step 2: initialize the downward recursion for F and F'
     !
     fg(nlambda,1) = sfmax
     fg(nlambda,2) = sfmax * rfmax
     !
     !  Step 3: downward recursion for F and F'
     !
     downward_recursion: do indl=nlambda,2,-1
       fg(indl-1,1) = (st(indl)*fg(indl,1) + fg(indl,2))/rt(indl)  ! F (i-1) = (S(i)*F(i) + F'(i))/R(i)
       fg(indl-1,2) =  st(indl)*fg(indl-1,1) - rt(indl)*fg(indl,1) ! F'(i-1) =  S(i)*F(i-1) - R(i)*F(i)
       ! Try rescaling if F or F' grows too large
       if (any(abs(fg(indl-1,1:2))>=sqrt(huge(1._rk)))) then
         fg(indl-1:,1:2) = fg(indl-1:,1:2) / sqrt(huge(1._rk))
       end if
     end do downward_recursion
     !
     !  Step 4: figure out (F+i*G)'/(F+i*G) for lambdamin
     !
     pq = cf2(lambdamin,x,eta)
     p  = real(pq,kind=rk)
     q  = aimag(pq)
     !
     !  Step 5: calculate G and G' for lambdamin
     !
     fg(1,3) = (fg(1,2) - p*fg(1,1))/q   ! (F' - p*F)/q
     fg(1,4) = p*fg(1,3) - q*fg(1,1)     ! p*G - q*F
     !
     !  If F and G are very different, expect problems with the accuracy!
     !
     if (abs(fg(1,3))>=1e3_rk*abs(fg(1,1)) .or. abs(fg(1,1))>=1e3_rk*abs(fg(1,3))) then
       write (out,"(/'WARNING: Accuracy loss of more than 3 digits expected in coulombFG')")
       write (out,"(' lambdamin = ',  1x,g32.24e3)") lambdamin
       write (out,"(' x         = ',  1x,g32.24e3)") x
       write (out,"(' eta       = ',  1x,g32.24e3)") eta
       write (out,"(' f,f''      = ',2(1x,g32.24e3))") fg(1,1:2)
       write (out,"(' g,g''      = ',2(1x,g32.24e3))") fg(1,3:4)
       write (out,"()")
     end if
     !
     !  Step 6: calculate the Wronskian. For normalized F and G it should be unity.
     !
     wrs = fg(1,2)*fg(1,3) - fg(1,1)*fg(1,4) ! F'*G - F*G'
     om1 = 1._rk/sqrt(abs(wrs))
     !
     !  Step 7: normalize available values of F and G
     !
     fg(:,1:2) = om1 * fg(:,1:2)  ! F and F' for all lambda
     !
     if (present(skip_g)) then
       if (skip_g) return ! Early exit: G is not needed.
     end if
     !
     fg(1,3:4) = om1 * fg(1,3:4)  ! G and G' for lambdamin
     !
     !  Step 8: upward recurrence for G and G'
     !
     upward_recursion: do indl=1,nlambda-1
       fg(indl+1,3) = (st(indl+1)*fg(indl,3) - fg(indl,4))/rt(indl+1)  ! G (i+1) = (S(i+1)*G(i) - G'(i))/R(i+1)
       fg(indl+1,4) =  rt(indl+1)*fg(indl,3) - st(indl+1)*fg(indl+1,3) ! G'(i+1) =  R(i+1)*G(i) - S(i+1)*G(i+1)
     end do upward_recursion
     !
     !  We are done
     !
   end subroutine coulombFG
   !
   !  Auxiliary routines below
   !
   real(rk) function brnR(lambda,eta)
     real(rk), intent(in) :: lambda, eta
     !
     brnR = sqrt(1._rk+(eta/lambda)**2)
   end function brnR
   !
   real(rk) function brnS(lambda,x,eta)
     real(rk), intent(in) :: lambda, x, eta
     !
     brnS = lambda/x + eta/lambda
   end function brnS
   !
   !  Continuous fraction 1 in eq. 27, evaluated using Steed's algorithm.
   !  Use Kahan summation algorith for h, just in case the sum converges slowly!
   !
   function cf1(lambda,x,eta,plusminus)
     real(rk), intent(in)  :: lambda, x, eta
     real(rk), intent(out) :: plusminus 
     real(rk)              :: cf1
     !
     real(rk)    :: h, hc, D, deltah
     integer(ik) :: n
     !
     h         = a(0_ik)
     hc        = 0
     D         = 1/b(1_ik)
     deltah    = a(1_ik)*D
     call kahan_add(h,hc,deltah)
     plusminus = sign(1._rk,D)
     iterate: do n=2,cf1_maxiter
       D         = 1/(D*a(n)+b(n))
       deltah    = (b(n)*D-1)*deltah
       call kahan_add(h,hc,deltah)
       plusminus = plusminus*sign(1._rk,D)
       if (abs(deltah)<=1e-6_rk*spacing(abs(h))) then  ! We effectively use double-double representation for the sum!
         cf1       = h
         return
       end if
     end do iterate
     stop 'coulomb_functions%cf1 - Continuous fraction failed to converge'
     !
     contains
     ! Kahan addition algorithm
     subroutine kahan_add(h,hc,deltah)
       real(rk), intent(inout) :: h      ! Sum
       real(rk), intent(inout) :: hc     ! Running correction to the sum
       real(rk), intent(in)    :: deltah ! Increment
       !
       real(rk) :: y, t, dt
       !
       y  = deltah - hc
       t  = h + y
       dt = t - h
       hc = dt - y
       h  = t
     end subroutine kahan_add
     ! Numerator
     real(rk) function a(n)
       integer(ik), intent(in) :: n ! Index, starting at 0
       real(rk)                :: p ! lambda + n
       !
       p = lambda + n
       if (n==0) then
         a = brnS(lambda+1,x,eta)
       else if (n==1) then
         a = -x*(p+1)*(p**2+eta**2)/p
       else if (n>1) then
         a = -x**2*(p**2-1)*(p**2+eta**2)
       else
         stop 'coulomb_functions%cf1%a - invalid order n'
       end if
     end function a
     ! Denominator; b(0) must be one!
     real(rk) function b(n)
       integer(ik), intent(in) :: n ! Index, starting at 0
       real(rk)                :: p ! lambda + n
       !
       p = lambda + n
       if (n==0) then
         b = 1
       else if (n>=1) then
         b = (2*p+1)*(p*(p+1)+eta*x)
       else
         stop 'coulomb_functions%cf1%b - invalid order n'
       end if
     end function b
   end function cf1
   !
   !  Continuous fraction 2 in eq. 29, evaluated using Steed's algorithm.
   !
   function cf2(lambda,x,eta)
     real(rk), intent(in) :: lambda, x, eta
     complex(rk)          :: cf2
     !
     complex(rk) :: h, hc, D, deltah
     integer(ik) :: n
     !
     h      = a(0_ik)
     hc     = 0
     D      = 1/b(1_ik)
     deltah = a(1_ik)*D
     call kahan_add(h,hc,deltah)
     iterate: do n=2,cf2_maxiter
       D      = 1/(D*a(n)+b(n))
       deltah = (b(n)*D-1)*deltah
       call kahan_add(h,hc,deltah)
       if (abs(deltah)<=1e-6_rk*spacing(abs(h))) then
         cf2 = h
         return
       end if
     end do iterate
     stop 'coulomb_functions%cf2 - Continuous fraction failed to converge'
     !
     contains
     ! Kahan addition algorithm
     subroutine kahan_add(h,hc,deltah)
       complex(rk), intent(inout) :: h      ! Sum
       complex(rk), intent(inout) :: hc     ! Running correction to the sum
       complex(rk), intent(in)    :: deltah ! Increment
       !
       complex(rk) :: y, t, dt
       !
       y  = deltah - hc
       t  = h + y
       dt = t - h
       hc = dt - y
       h  = t
     end subroutine kahan_add
     ! Numerator
     complex(rk) function a(n)
       integer(ik), intent(in) :: n ! Index, starting at 0
       !
       if (n==0) then
         a = ii*(1-eta/x)
       else if (n>=1) then
         a = (ii*eta-lambda+n-1)*(ii*eta+lambda+n)
         if (n==1) a = a*ii/x
       else
         stop 'coulomb_functions%cf2%a - invalid order n'
       end if
     end function a
     ! Denominator; b(0) must be one!
     complex(rk) function b(n)
       integer(ik), intent(in) :: n ! Index, starting at 0
       !
       if (n==0) then
         b = 1
       else if (n>=1) then
         b = 2*(x-eta+ii*n)
       else
         stop 'coulomb_functions%cf2%b - invalid order n'
       end if
     end function b
   end function cf2
   !
   function coulombBound(n,l,r,Z) result(rad)
     integer(ik), intent(in) :: n, l   ! Major and angular quantum numbers
     real(rk), intent(in)    :: r      ! Distance from the origin 
     real(rk), intent(in)    :: Z      ! The nuclear charge (must be positive!)
     real(rk)                :: rad
     !
     real(rk)  :: x
     !
     if (n<l+1 .or. l<0 .or. r<0 .or. Z<=0) then
       write (out,"('coulombBound(',i0,',',i0,',',g0.8,',',g0.8,') is non-physical')") n, l, r, Z
       call flush_wrapper(out)
       stop 'coulomb_functions%coulombBound - unphysical'
     end if
     x      = 2*Z*r/n
     ! This naive expression may have problems for large n and l values !
     rad =   sqrt((2*Z/n)**3 * MathFactorial(n-l-1_ik) / (2_ik*n*MathFactorial(n+l))) &
          * exp(-0.5_rk * x) * x**l * l2g(n-l-1_ik,2_ik*l+1_ik,x)
     contains
     !
     !  Assocuated Laguerre polynomial, special case for hydrogenic wavefunctions
     !
     real(rk) function l2g(a,b,z) 
       integer(ik), intent(in) :: a, b
       real(rk), intent(in)    :: z
       !
       integer(ik) :: n
       real(rk)    :: m1, m2
       !
       l2g = 1_ik ; m1 = 0_ik
       recurse: do n=1,a
         m2  = m1 ; m1 = l2g
         l2g = (b+2_ik*n-1_ik-z)*m1/n - (n+b-1_ik)*m2/n
       end do recurse
     end function l2g
   end function coulombBound

 end module coulomb_functions
!
!  tsurf-coulomb-waves.nb can generate compatible input for testing.
!  See "coulomb-test.inp" for some pre-packaged tests.
!
!program test
!  use coulomb_functions
!  !
!  integer(ik), parameter :: ltop = 5000_ik
!  real(rk)    :: Z, k, r
!  integer(ik) :: lmin, lmax, l, laps
!  real(rk)    :: f    (0:ltop), fp    (0:ltop)
!  real(rk)    :: ref_f(0:ltop), ref_fp(0:ltop)
!  real(rk)    :: err_f, err_fp
!  real(rk)    :: rel_f, rel_fp
!  !
!  tests: do
!    read(*,*,end=1000) Z, k, r, lmin, lmax
!    if (lmin<0 .or. lmax>ltop .or. lmax<lmin) stop 'coulomb_functions%test - bad input'
!    read(*,*,end=1000) ref_f (lmin:lmax)
!    read(*,*,end=1000) ref_fp(lmin:lmax)
!    !
!    run_laps: do laps=1,1
!    call coulombF(lmin,lmax,Z,k,r,f(lmin:lmax),fp(lmin:lmax))
!    end do run_laps
!    !
!    report: do l=lmin,lmax
!      call error(f (l),ref_f (l),err_f, rel_f)
!      call error(fp(l),ref_fp(l),err_fp,rel_fp)
!      write (6,"(1x,i4,1x,f12.5,1x,f12.5,1x,f12.5,4(1x,g13.5e3),2(1x,g45.35e3))") &
!             l, Z, k, r, rel_f, rel_fp, err_f, err_fp, f(l), fp(l)
!    end do report
!  end do tests
!  1000 continue
!  !
!  contains
!  subroutine error(v,ref,err,rel)
!    real(rk), intent(in)  :: v, ref   ! Value and reference
!    real(rk), intent(out) :: err, rel ! Error and relative error
!    !
!    err = v - ref
!    rel = abs(err)/max(abs(ref),spacing(ref),spacing(err))
!  end subroutine error
!end program test
!program test2
!  use coulomb_functions
!  !
!  real(rk)              :: Z, k, dr, r
!  integer(ik)           :: lmax, npts, ipt, lv
!  real(rk), allocatable :: c_f(:), c_g(:)
!  !
!  read(*,*) Z, k, dr, lmax, npts
!  allocate (c_f(0:lmax),c_g(0:lmax))
!  !
!  radial_points: do ipt=1,npts
!    r = dr*ipt
!    call coulombF(0_ik,lmax,Z,k,r,c_f,c_g)
!    write (6,"(1x,f14.5,300(1x,i2,1x,g26.15,1x,g26.15))") r, (lv, c_f(lv), c_g(lv), lv=0,lmax)
!  end do radial_points
!end program test2
!program test3
!  use coulomb_functions
!  !
!  integer(ik) :: n, l
!  real(rk)    :: Z, r, rad, r0, re, dr
!  !
!  read(*,*) n, l, Z, r0, re, dr
!  !
!  r = r0
!  radial_points: do while (r<=re)
!    rad = coulombBound(n,l,r,Z)
!    write (6,"(2(1x,g24.14e3))") r, rad
!    r = r + dr
!  end do radial_points
!end program test3
