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
!  Sundry mathematical definitions.
!  This is a subset of the routines exported by multigrid's math.f90
!
module math
  use accuracy
  use constants
  !$ use OMP_LIB
  implicit none
  private
  public MathFactorial, MathLogFactorial
  public MathJacobiPn
  public MathDoubleFactorial, MathLogDoubleFactorial
  public MathBinomial, MathPochhammer
  public Math3J, MathLogSum
  public MathRotationMatrix, MathYJMRotationMatrix
  public MathAllYLM2
  public MathLogGamma
  public MathInterpolate
  public rcsid_math
!
  character(len=clen), save :: rcsid_math = "$Id: math.f90,v 1.15 2023/06/09 14:10:24 ps Exp $"
!
!
  integer(ik), parameter       :: factorial_slack = 5         ! Extra factorials to produce while filling the cache
  integer(ik), save            :: factorial_max = -1          ! Largest value of the factorials cached in the table
  real(rk), allocatable, save  :: factorial_table(:)
  integer(ik), save            :: log_factorial_max = -1      ! Largest value of the factorials cached in the table
  real(rk), allocatable, save  :: log_factorial_table(:)
  integer(ik), save            :: dfactorial_max = -2         ! Largest value of the double factorials cached in the table
  real(rk), allocatable, save  :: dfactorial_table(:)
  integer(ik), save            :: log_dfactorial_max = -2     ! Largest value of the double factorials cached in the table
  real(rk), allocatable, save  :: log_dfactorial_table(:)
!
  contains
  !
  !  External interfaces
  !
  function MathFactorial(n) result(v)
    integer(ik), intent(in)        :: n
    real(rk)                       :: v
    !
    if (n<0) stop 'math%MathFactorialReal - domain error'
    if (n>factorial_max) call fill_factorial_table(n+factorial_slack)
    v = factorial_table(n)
  end function MathFactorial
  !
  function MathLogFactorial(n) result(v)
    integer(ik), intent(in)        :: n
    real(rk)                       :: v
    !
    if (n<0) stop 'math%MathLogFactorial - domain error'
    if (n>log_factorial_max) call fill_log_factorial_table(n+factorial_slack)
    v = log_factorial_table(n)
  end function MathLogFactorial
  !
  function MathDoubleFactorial(n) result(v)
    integer(ik), intent(in)        :: n
    real(rk)                       :: v
    !
    if (n<-1) stop 'math%MathDoubleFactorial - domain error'
    if (n>dfactorial_max) call fill_dfactorial_table(n+factorial_slack)
    v = dfactorial_table(n)
  end function MathDoubleFactorial
  !
  function MathLogDoubleFactorial(n) result(v)
    integer(ik), intent(in) :: n
    real(rk)                :: v
    !
    if (n<-1) stop 'math%MathLogDoubleFactorial - domain error'
    if (n>log_dfactorial_max) call fill_log_dfactorial_table(n+factorial_slack)
    v = log_dfactorial_table(n)
  end function MathLogDoubleFactorial
  !
  !  Pochhammer function: product of (n) integers starting at (a)
  !
  function MathPochhammer(a,n) result(v)
    integer(ik), intent(in) :: a, n
    real(rk)                :: v
    !
    integer(ik) :: aa     ! Starting integer of the equivalent positive sequence
    logical     :: minus
    !
    if (n<0) stop 'math%MathPochhammer - domain error'
    !
    if (n==0) then
      v = 1._rk
      return
    end if
    !
    if (a<=0) then
      !
      !  Catch sequences containing zero factors: these are always zero
      !
      if (a+n-1>=0) then
        v = 0._rk
        return
      end if
      aa    = -(a+n-1)
      minus = mod(n,2)==1
    else
      aa    = a
      minus = .false.
    end if
    !
    v = MathLogFactorial(aa+n-1)-MathLogFactorial(aa-1)
    v = exp(v)
    if (minus) v = -v
  end function MathPochhammer
  !
  !  Binomial coefficients
  !
  function MathBinomial(n,m) result(cnm)
    integer(ik), intent(in)        :: n, m
    real(rk)                       :: cnm
    !
    if (n<0 .or. m<0 .or. m>n) stop 'MathBinomialIntegerReal - domain error'
    cnm = exp(MathLogFactorial(n)-MathLogFactorial(n-m)-MathLogFactorial(m))
  end function MathBinomial
  !
  !  Evaluate log(exp(a)+exp(b)), where a and b may be too large to fit 
  !  in the exponent range.
  !
  function MathLogSum(a,b) result(c)
    real(rk), intent(in) :: a, b
    real(rk)             :: c
    !
    if (a>=b) then
      c = a + log(1._rk + exp(b-a))
    else
      c = b + log(1._rk + exp(a-b))
    end if
  end function MathLogSum
  !
  !  Use recurrence with respect to degree, see Abramowitz & Stegun 22.7.1
  !
  function MathJacobiPn(n,alp,bet,x) result(pn)
    integer(ik), intent(in) :: n        ! Order of the polynomial
    real(rk), intent(in)    :: alp, bet ! Powers of the weight function, must be > -1
    real(rk), intent(in)    :: x        ! Coordinate, abs(x)<=1
    real(rk)                :: pn
    !
    real(rk) :: tab(0:n)
    !
    call jacobiPn_table(n,alp,bet,x,tab)
    pn = tab(n)
  end function MathJacobiPn
  !
  !  Computes Wigner 3J symbols. The code below is a direct implementation
  !  of L&L 3J formulae. The accuracy of this routine is reduced relative to
  !  that is theoretically possible, due to the use of logarithms. The routine
  !  I had in MNDO99 is more accurate and can handle broader range of J values.
  !
  function Math3J(j1,j2,j3,m1,m2,m3) result(v)
    integer(ik), intent(in) :: j1, j2, j3  ! / J1 J2 J3 \ 3-j
    integer(ik), intent(in) :: m1, m2, m3  ! \ M1 M2 M3 / 
    real(rk)                :: v
    !
    integer(ik) :: ij0, ij1, ij2, ij3, im1a, im1b, im2a, im2b, im3a, im3b
    integer(ik) :: t1, t2, t3
    integer(ik) :: z, minz, maxz
    real(rk)    :: logscale, logterm
    !
    !  Before we do anything, check whether this 3J symbol satisfies the
    !  vector addition constraints
    !
    ij0  =   j1 + j2 + j3 + 1
    ij1  =   j1 + j2 - j3
    ij2  =   j1 - j2 + j3
    ij3  = - j1 + j2 + j3
    im1a =   j1 - m1 ; im1b = j1 + m1
    im2a =   j2 - m2 ; im2b = j2 + m2
    im3a =   j3 - m3 ; im3b = j3 + m3
    if (ij1<0 .or. ij2<0 .or. ij3<0 .or. im1a<0 .or. im1b<0 .or. im2a<0 .or. im2b<0 .or. im3a<0 .or. im3b<0 .or. m1+m2+m3/=0) then
      v = 0
      return
    end if
    !
    logscale = MathLogFactorial(ij1)  + MathLogFactorial(ij2)  + MathLogFactorial(ij3)  &
             + MathLogFactorial(im1a) + MathLogFactorial(im1b) + MathLogFactorial(im2a) &
             + MathLogFactorial(im2b) + MathLogFactorial(im3a) + MathLogFactorial(im3b) &
             - MathLogFactorial(ij0)
    logscale = 0.5_rk * logscale
    !
    t1   = j2 - j3 - m1
    t2   = j1 + m2 - j3
    t3   = j1 - j2 - m3
    minz = max(0_ik,t1,t2)
    maxz = min(ij1,im1a,im2b)
    v = 0
    sum_terms: do z=minz,maxz,1
      logterm = logscale - MathLogFactorial(z)      - MathLogFactorial(ij1-z)  - MathLogFactorial(im1a-z) &
                         - MathLogFactorial(im2b-z) - MathLogFactorial(z-t1)   - MathLogFactorial(z-t2)
      if (abs(logterm)>=0.9_rk*log(huge(1._rk))) then
        write (out,"('Math3J: Intermediate logarithm ',g12.5,' exceeds the real(rk) dynamic range.')") logterm
        write (out,"('Math3J: The 3J arguments were: ',6i10)") j1, j2, j3, m1, m2, m3
        stop 'math%Math3J - exceeded dynamic range'
      end if
      if (mod(z+t3,2)==0) then
        v = v + exp(logterm)
      else
        v = v - exp(logterm)
      end if
    end do sum_terms
    !
  end function Math3J
  !
  !  Given Euler angles, construct rotation matrix for coordinate axes 
  !  in 3D space.  The Euler angles are defined as follows:
  !   1. Rotate coordinate axes by alpha around the Z axis
  !   2. Rotate axes by beta around the new Y axis
  !   3. Rotate axes by gamma around the new Z axis
  !  This definition of the Euler angles matches the definition from
  !  section 58 of L&L III.
  !
  !  Note that prior to June 10th, 2010 the definition of all three Euler 
  !  angles in this routine used a wrong sign, corresponding to anti-screw
  !  rotation sense. Thanks, Mike.
  !
  !  If you would like to rotate an object instead of the axes, take the
  !  transpose or (equivalently) replace (alpha,beta,gamma) with
  !  (-gamma,-beta,-alpha).
  !
  !  Some useful special cases are:
  !
  !    a     b     c  
  !   ---   ---   --- 
  !    x     0     0   Rotate coordinate system by x around the Z axis
  !    0     0     x   Rotate coordinate system by x around the Z axis
  !    0     x     0   Rotate coordinate system by x around the Y axis
  !  -pi/2   x   pi/2  Rotate coordinate system by x around the X axis
  !      
  subroutine MathRotationMatrix(euler_angles,mat)
    real(rk), intent(in)  :: euler_angles(3) ! Euler rotation angles: alpha, beta, and gamma
    real(rk), intent(out) :: mat(3,3)        ! Rotation matrix
    !
    real(rk) :: a, b, g, rma(3,3), rmb(3,3), rmg(3,3), rmba(3,3)
    real(rk) :: sin_a, cos_a, sin_b, cos_b, sin_g, cos_g
    !
    a = euler_angles(1) ; sin_a = sin(a) ; cos_a = cos(a)
    b = euler_angles(2) ; sin_b = sin(b) ; cos_b = cos(b)
    g = euler_angles(3) ; sin_g = sin(g) ; cos_g = cos(g)
    !
    !  We can't use array constructors here: gfortran ends up creating array temporaries
    !  in this case.
    !
    rma(1,1) =  cos_a ; rma(1,2) = sin_a ; rma(1,3) =  0._rk
    rma(2,1) = -sin_a ; rma(2,2) = cos_a ; rma(2,3) =  0._rk
    rma(3,1) =  0._rk ; rma(3,2) = 0._rk ; rma(3,3) =  1._rk
    !
    rmb(1,1) =  cos_b ; rmb(1,2) = 0._rk ; rmb(1,3) = -sin_b
    rmb(2,1) =  0._rk ; rmb(2,2) = 1._rk ; rmb(2,3) =  0._rk
    rmb(3,1) =  sin_b ; rmb(3,2) = 0._rk ; rmb(3,3) =  cos_b
    !
    rmg(1,1) =  cos_g ; rmg(1,2) = sin_g ; rmg(1,3) =  0._rk
    rmg(2,1) = -sin_g ; rmg(2,2) = cos_g ; rmg(2,3) =  0._rk
    rmg(3,1) =  0._rk ; rmg(3,2) = 0._rk ; rmg(3,3) =  1._rk
    !
    rmba = matmul(rmb,rma)
    mat  = matmul(rmg,rmba)
  end subroutine MathRotationMatrix
  !
  !  Rotation matrix for angular momentum eigenfunctions, following L&L III Eq. 58.10
  !  Both integer and half-integer J values are OK.
  !
  !  The resulting rotation matrix is accurate to 2ulp for multiplicities up to 6,
  !  with error increasing to 4ulp for multiplicity 20. It loses about 11 decimal places
  !  of accuracy for multiplicity 81, and overflows IEEE double at higher multiplicities.
  !
  !  Note that the rotation matrix uses somewhat weird conventions: it rotates transposed
  !  harmonics from the primed coordinate system defined by the Euler angles back into
  !  the lab system:
  !  
  !    Y(L,M) = Sum Y(L,M') D(M',M)
  !
  !  Furthermore, it looks like the expression for the Wigner matrix in the 5th Russian
  !  edition of L&L is actually incorrect. To get the correct expression, it is necessary
  !  to change the sign of the Euler beta angle. The code below is a literal implementation
  !  of L&L 58.10, so don't forget to flip the sign of beta when calling it!
  !
  subroutine MathYJMRotationMatrix(euler_angles,mult,mat)
    real(rk), intent(in)     :: euler_angles(3) ! Euler rotation angles: alpha, beta, and gamma
                                                 ! See comments in MathRotationMatrix
    integer(ik), intent(in)   :: mult            ! Multipliplicity of the angular-momentum state,
                                                 ! mult = 2*j+1
    complex(rk), intent(out) :: mat(:,:)        ! Rotation matrix
    !
    real(rk)    :: a, b, g
    real(rk)    :: cosb2, sinb2
    complex(rk) :: expa2, expg2
    integer(ik)  :: j2, m2, mp2
    integer(ik)  :: im, imp
    !
    if (mult<1) then
      stop 'math%MathYJMRotationMatrix - multiplicity: domain error'
    end if
    if (size(mat,dim=1)/=mult .or. size(mat,dim=2)/=mult) then
      stop 'math%MathYJMRotationMatrix - rotation matrix dimensions do not match multiplicity'
    end if
    !
    a = euler_angles(1)
    b = euler_angles(2)
    g = euler_angles(3)
    !
    !  We need to take special care when angle beta approaches n*pi. For these angles,
    !  
    !
    sinb2 = sin(0.5_rk*b)
    cosb2 = cos(0.5_rk*b)
    !
    expa2 = exp(cmplx(0._rk,0.5_rk*a,kind=rk))
    expg2 = exp(cmplx(0._rk,0.5_rk*g,kind=rk))
    !
    j2  = mult - 1
    mp2 = -j2
    row_mp: do imp=1,mult
      m2 = -j2
      column_m: do im=1,mult
        mat(imp,im) = sqrt(hf(j2+mp2)*hf(j2-mp2)/(hf(j2+m2)*hf(j2-m2))) &
                    * modJacobiPn((j2-mp2)/2,(mp2-m2)/2,(mp2+m2)/2,sinb2,cosb2) &
                    * expa2**m2 * expg2**mp2
        m2 = m2 + 2
      end do column_m
      mp2 = mp2 + 2
    end do row_mp
    !
    contains
    !
    !  (n/2)! where n must be even
    !
    function hf(n) result(f)
      integer(ik), intent(in) :: n ! Must be even
      real(rk)               :: f
      !
      if (mod(n,2)/=0) stop 'math%MathYJMRotationMatrix%hf - domain error'
      f = MathFactorial(n/2)
    end function hf
    !
    !  Specialized derivative of Jacobi P polynomial:
    !
    !    y^a x^b JacobiP(n,a,b,x^2-y^2)
    !
    !  where x^2+y^2 is equal to 1. Care should be taken in evaluating this function
    !  when either a or b are negative: standard recursion with respect to degree 
    !  becomes undefined in this case.
    !
    !  As the result, we have to use series expansion around +1/-1 argument to
    !  evaluate this function.
    !
    function modJacobiPn(n,a,b,y,x) result(jp)
      integer(ik), intent(in) :: n, a, b ! Parameters of the Jacobi P polynomial
      real(rk), intent(in)   :: y, x    ! Coordinate
      real(rk)               :: jp
      !
      integer(ik) :: k
      !
      if (abs(y)<abs(x)) then
        !
        !  Small y, expand JacobiP around z=+1
        !
        jp = 0
        expand_plus1: do k=max(0,-a),n
          jp = jp + MathPochhammer(a+k+1,n-k) * MathPochhammer(-n,k) * MathPochhammer(a+b+n+1,k) &
                  * y**(2*k+a) / MathFactorial(k)
        end do expand_plus1
        jp = x**b * jp / MathFactorial(n)
      else 
        !
        !  Small x, expand JacobiP around z=-1
        !
        jp = 0
        expand_minus1: do k=max(0,-b),n
          jp = jp + MathPochhammer(b+k+1,n-k) * MathPochhammer(-n,k) * MathPochhammer(a+b+n+1,k) &
                  * x**(2*k+b) / MathFactorial(k)
        end do expand_minus1
        jp = y**a * jp / MathFactorial(n)
        if (mod(n,2)==1) jp = -jp
      endif
!       !
!       !  The general case; do not use
!       !
!       jp = y**a * x**b * MathJacobiPn(n,real(a,kind=rk),real(b,kind=rk),x**2-y**2)
    end function modJacobiPn
  end subroutine MathYJMRotationMatrix
  !
  !  Our old spherical harmonics code has a problem with large values of L, where 
  !  we need to multiply through many large and small factors to give an answer on
  !  the order of 1. The code below integrates evaluation of the prefactor into
  !  the recurrences for the associated Legendre polynomials, so that not large
  !  factors arise. By switching to a different recursion, we can also support 
  !  the case of a sub-set of possible M values being required.
  !
  subroutine MathAllYLM2(l_max,m_min,m_max,dir,ylm,phase)
    integer(ik), intent(in)                :: l_max   ! Desired maximum angular momentum
    integer(ik), intent(in)                :: m_min   ! Desired minimum Z projection of the angular momentum
    integer(ik), intent(in)                :: m_max   ! Desired maximum Z projection of the angular momentum
    real(rk), intent(in)                   :: dir(3)  ! (Unnormalized) direction vector
    complex(rk), intent(out)               :: ylm(m_min:m_max,0:l_max)
    character(len=*), intent(in), optional :: phase   ! Phase convention to use. Recognized values are:
                                                      ! 'L&L'    = Landau and Lifshitz v. 3 phase convention (the default)
                                                      ! 'Arfken' = Arfken and Mathematica phase convention
    !
    real(rk), parameter :: y00 = 1._rk/sqrt(4._rk*pi) ! Y00 does not depend on the spherical angles
    integer(ik)         :: k, l, m, m_low, m_high
    real(rk)            :: sqrtn(0:2*l_max+1)  ! Square roots of small integers
    real(rk)            :: rsqrn(1:2*l_max+1)  ! Inverse square roots of integers
    real(rk)            :: sinth, costh        ! Sine and cosine of the spherical theta angle
    real(rk)            :: r, xymod
    complex(rk)         :: xpy                 ! exp(+(0._rk,1._rk)*phi), where phi is the spherical phi angle
    complex(rk)         :: xmy                 ! exp(-(0._rk,1._rk)*phi)
    complex(rk)         :: vlm                 ! Current YLM in simultaneous L,M recurrence
    !
    !  Sanity testing
    !
    if (l_max<0 .or. abs(m_min)>l_max .or. abs(m_max)>l_max .or. m_min>m_max) then
      write (out,"('MathAllYLM2: l_max= ',i0,' m_min= ',i0,' m_max= ',i0)") l_max, m_min, m_max
      call flush_wrapper(out)
      stop 'math%MathAllYLM2 - called with invalid L,M arguments'
    end if
    !
    !  Spherical parameters
    !
    r = sqrt(sum(dir**2))
    if (r<=0) then
      stop 'math%MathAllYLM2 - direction cannot be zero'
    end if
    xpy   = cmplx(dir(1), dir(2),kind=rk)
    xmy   = cmplx(dir(1),-dir(2),kind=rk)
    xymod = abs(xpy)
    if (xymod>0) then
      xpy = xpy / xymod
      xmy = xmy / xymod
    else
      xpy = 1
      xmy = 1
    end if
    costh = min(1._rk,max(-1._rk,dir(3)/r)) ! Force sine and cosine into the allowed range,
    sinth = min(1._rk,max( 0._rk,xymod /r)) ! to guard against round-off errors
    !
    !  Precompute some of the factors we need
    !
    sqrtn(0) = 0
    fill_square_roots: do k=1,2*l_max+1
      sqrtn(k) = sqrt(real(k,kind=rk))
      rsqrn(k) = 1._rk / sqrtn(k)
    end do fill_square_roots
    !
    if (m_min<=0 .and. m_max>=0) ylm(0,0) = y00
    !
    !  Simultaneously increase L and M, starting at zero
    !
    vlm = y00
    upward_diagonal_lm: do l=1,l_max
      vlm = vlm * (0._rk,1._rk) * sqrtn(2*l+1)*rsqrn(2*l)*xpy*sinth
      if (m_min<=l .and. m_max>=l) ylm(l,l) = vlm
    end do upward_diagonal_lm
    !
    !  Simultaneously increase L and decrease M, starting at zero
    !
    vlm = y00
    downward_diagonal_lm: do l=1,l_max
      vlm = vlm * (0,-1) * sqrtn(2*l+1)*rsqrn(2*l)*xmy*sinth
      if (m_min<=-l .and. m_max>=-l) ylm(-l,l) = vlm
    end do downward_diagonal_lm
    !
    !  For each M increase L, starting at L=M
    !
    upward_l: do l=1,l_max
      m_low  = max(m_min,-(l-1))  ! We also must satisfy |m|<=l-1 for this recursion
      m_high = min(m_max, (l-1))
      upward_l_m_loop: do m=m_low,m_high
        vlm = ylm(m,l-1) * (0._rk,1._rk)*costh*sqrtn(2*l-1)
        if (l>=abs(m)+2) then
          vlm = vlm + ylm(m,l-2) * rsqrn(2*l-3)*sqrtn(l-1-m)*sqrtn(l-1+m) ! Spurious gfortran warning here
        end if
        ylm(m,l) = vlm * sqrtn(2*l+1)*rsqrn(l-m)*rsqrn(l+m)
      end do upward_l_m_loop
    end do upward_l
    !
    if (present(phase)) then
      select case (phase)
        case default; stop 'math%MathAllYLM2 - Invalid phase convention'
        case ('L&L') ! Do nothing
        case ('Arfken')
          adjust_l: do l=0,l_max
            adjust_m: do m=m_min,m_max
              ! For some reason, Intel Fortran has real trouble generating decent code for the 
              ! commented line below. Let's do a hack here ...
              ! ylm(m,l) = ylm(m,l) * (0._rk,1._rk)**(-2*m-l)
              ylm(m,l) = ylm(m,l) * ipow(modulo(-2*m-l,4))
            end do adjust_m
          end do adjust_l
      end select
    end if
  end subroutine MathAllYLM2
  !
  !  Recurrence from A&S 22.7.1. These recurrences are not completely numerically
  !  stable, and lose up to 3 decimal digits for n>10.
  !  Unfortunately, these recurrences do not work for negative integer a and b:
  !  sooner or later, we always hit division by zero. Oops.
  !
  subroutine jacobiPn_table(nmax,a,b,x,pn)
    integer(ik), intent(in) :: nmax  ! Maximum order of the Legendre polynomials desired
    real(rk), intent(in)    :: a, b  ! alpha and beta parameters of the weight function
    real(rk), intent(in)    :: x     ! Coordinate at which P values are needed
    real(rk), intent(out)   :: pn(:) ! Values of LegendreP from n=0 to n=nmax
    !
    integer(ik) :: n
    real(rk)    :: a1n, a2n, a3n, a4n
    !
    if (nmax<0) stop 'math%jacobiPn_table - negative-order polynomial requested'
    pn(1) = 1._rk
    if (nmax<1) return
    pn(2) = 0.5_rk*((a-b) + (2._rk+a+b)*x)
    !
    !  A&S 22.7 recursions are written for a polynomial of degree (n+1).
    !  To avoid confusion, we'll do the same. Note that the a3n coefficient in A&S 22.7.1
    !  is written using non-standard series notation, which is likely to cause confusion ...
    !  Keep in mind that n-th degree polynomial is at pn(n+1)
    !
    n_recursion: do n=1,nmax-1
      a1n = 2 * (n+1)*(n+a+b+1) * (2*n+a+b)
      a2n = (2*n+a+b+1) * (a**2-b**2)
      a3n = (2*n+a+b) * (2*n+a+b+1) * (2*n+a+b+2)
      a4n = 2 * (n+a) * (n+b) * (2*n+a+b+2)
      !
      pn(n+2) = ( (a2n+a3n*x)*pn(n+1) - a4n*pn(n) ) / a1n
    end do n_recursion
  end subroutine jacobiPn_table
  !
  subroutine fill_factorial_table(nmax)
    integer(ik), intent(in) :: nmax
    integer(ik)             :: n, alloc
    real(rk)                :: fac
    !
    !$ if (omp_in_parallel()) then
    !$   stop 'math%fill_factorial_table - unsafe call to MathFactorial'
    !$ end if
    !
    if (factorial_max>=0) then
      deallocate (factorial_table,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i10,' deallocating factorial table')") alloc
        stop 'math%fill_factorial_table - deallocate'
      end if
    end if
    !
    n   = 0
    fac = 1._rk
    !
    allocate (factorial_table(0:nmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('Error ',i10,' allocating ',i10,'-element factorial table')") & 
             alloc, nmax
      stop 'math%fill_factorial_table - allocate'
    end if
    !
    fill_factorials: do while(n<=nmax-1)
      factorial_table(n) = fac
      n = n + 1
      !
      if (huge(fac)/n<=fac) then
        write (out,"(1x,i10,'! would exceed dynamic range of the chosen real kind')") n
        stop 'math%fill_factorial_table - range exceeded'
      end if
      fac = fac * n
    end do fill_factorials
    factorial_table(n) = fac
    !
    factorial_max = nmax
  end subroutine fill_factorial_table
  !
  subroutine fill_log_factorial_table(nmax)
    integer(ik), intent(in) :: nmax
    integer(ik)             :: n, alloc
    real(rk)                :: fac
    !
    !$ if (omp_in_parallel()) then
    !$   stop 'math%fill_factorial_table - unsafe call to MathLogFactorial'
    !$ end if
    !
    if (log_factorial_max>=0) then
      deallocate (log_factorial_table,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i10,' deallocating log-factorial table')") alloc
        stop 'math%fill_factorial_table - deallocate'
      end if
    end if
    !
    n   = 0
    fac = 0._rk
    !
    allocate (log_factorial_table(0:nmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('Error ',i10,' allocating ',i10,'-element log-factorial table')") & 
             alloc, nmax
      stop 'math%fill_log_factorial_table - allocate'
    end if
    !
    fill_factorials: do while(n<=nmax-1)
      log_factorial_table(n) = fac
      n = n + 1
      !
      fac = fac + log(real(n,kind=rk))
    end do fill_factorials
    log_factorial_table(n) = fac
    !
    log_factorial_max = nmax
  end subroutine fill_log_factorial_table
  !
  subroutine fill_dfactorial_table(nmax)
    integer(ik), intent(in) :: nmax
    integer(ik)             :: n, alloc
    !
    !$ if (omp_in_parallel()) then
    !$   stop 'math%fill_dfactorial_table - unsafe call to MathDoubleFactorial'
    !$ end if
    !
    if (dfactorial_max>=0) then
      deallocate (dfactorial_table,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i10,' deallocating double factorial table')") alloc
        stop 'math%fill_dfactorial_table - deallocate'
      end if
    end if
    !
    allocate (dfactorial_table(-1:nmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('Error ',i10,' allocating ',i10,'-element double factorial table')") & 
             alloc, nmax
      stop 'math%fill_dfactorial_table - allocate'
    end if
    !
    dfactorial_table(-1:1) = 1._rk
    n = 2
    fill_factorials: do while(n<=nmax)
      if (huge(1._rk)/n<=dfactorial_table(n-2)) then
        write (out,"(1x,i10,'!! would exceed dynamic range of the chosen real kind')") n
        stop 'math%fill_dfactorial_table - range exceeded'
      end if
      dfactorial_table(n) = dfactorial_table(n-2) * n
      n = n + 1
    end do fill_factorials
    !
    dfactorial_max = nmax
  end subroutine fill_dfactorial_table
  !
  subroutine fill_log_dfactorial_table(nmax)
    integer(ik), intent(in) :: nmax
    integer(ik)             :: n, alloc
    !
    !$ if (omp_in_parallel()) then
    !$   stop 'math%fill_dfactorial_table - unsafe call to MathLogDoubleFactorial'
    !$ end if
    !
    if (log_dfactorial_max>=0) then
      deallocate (log_dfactorial_table,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i10,' deallocating log-double factorial table')") alloc
        stop 'math%fill_dfactorial_table - deallocate'
      end if
    end if
    !
    allocate (log_dfactorial_table(-1:nmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('Error ',i10,' allocating ',i10,'-element log-double factorial table')") & 
             alloc, nmax
      stop 'math%fill_log_dfactorial_table - allocate'
    end if
    !
    log_dfactorial_table(-1:1) = 0._rk
    n = 2
    fill_factorials: do while(n<=nmax)
      log_dfactorial_table(n) = log_dfactorial_table(n-2) + log(real(n,kind=rk))
      n = n + 1
    end do fill_factorials
    !
    log_dfactorial_max = nmax
  end subroutine fill_log_dfactorial_table
  !
  !
  !  Evaluation of the gamma function of complex argument seems to be a bit
  !  complicated. Let's do something very straightforward, and hopefully accurate
  !  enough:
  !    1. Use recurrence relations to _increase_ the real part of the argument
  !       until modulus is large enough to let us use asymprotic expressions
  !    2. Then apply the Stirling formula to get the result.
  !
  !  The resulting routine is not fast, but is reasonably accurate: for
  !  arguments between 0 and 10, the relative error does not exceed 5e-14
  !
  function MathLogGamma(z) result(v)
    complex(rk), intent(in) :: z ! Argument of the gamma function
    complex(rk)             :: v ! Gamma function
    !
    !  Stirling formula coefficients
    !
    real(rk), parameter :: c01 =               1._rk/                12._rk
    real(rk), parameter :: c02 =               1._rk/               288._rk
    real(rk), parameter :: c03 =            -139._rk/             51840._rk
    real(rk), parameter :: c04 =            -571._rk/           2488320._rk
    real(rk), parameter :: c05 =          163879._rk/         209018880._rk
    real(rk), parameter :: c06 =         5246819._rk/       75246796800._rk
    real(rk), parameter :: c07 =      -534703531._rk/      902961561600._rk
    real(rk), parameter :: c08 =     -4483131259._rk/    86684309913600._rk
    real(rk), parameter :: c09 = 432261921612371._rk/514904800886784000._rk
    !
    complex(rk) :: zr    ! Reduced argument
    complex(rk) :: logs  ! Scaling coefficient needed to reduce the result of Stirling formula
    complex(rk) :: logv  ! Logarithm of the result
    real(rk)    :: zcut  ! Minimal safe argument for the Stirling formula
    complex(rk) :: vs    ! Series part of the Stirling formula
    !
    !  To get accurate results from Stirling formula, we must make sure that
    !  the modulus of the argument is large enough to make the last contribution
    !  to the expansion small enough.
    !
    zcut = (c09/spacing(1._rk))**(1._rk/9._rk)
    ! write (out,"(' zcut = ',f25.15)") zcut
    logs = 0._rk
    zr   = z
    inflate_z: do while(abs(zr)<zcut)
      ! logs = logs + log(zr*(zr+1)*(zr+2)) ! This "optimization" does not work, because of the branch cuts of log()
      logs = logs + log(zr) + log(zr+1) + log(zr+2)
      zr   = zr + 3
    end do inflate_z
    ! write (out,"(' zr = ',2(1x,g25.15),' logs = ',2(1x,g25.15))") zr, logs
    !
    !  It is safe to use Stirling formula now
    !
    vs = 1._rk + (c01 + (c02 + (c03 + (c04 + (c05 + (c06 + (c07 + (c08 + c09/zr)/zr)/zr)/zr)/zr)/zr)/zr)/zr)/zr
    ! write (out,"(' vs = ',2(1x,g25.15))") vs
    logv = log(sqrt(2._rk*pi)) - zr + (zr-0.5_rk)*log(zr) + log(vs) - logs
    ! write (out,"(' logv = ',2(1x,g25.15))") logv
    !
    ! v = exp(logv)
    v = logv
  end function MathLogGamma
  !
  !  MathInterpolate is a transcription of the routine 3.1 in Numerical Recipes.
  !
  function MathInterpolate(x,xtab,vtab) result (v)
    real(xk), intent(in) :: x       ! Point where the interpolant is desired
    real(xk), intent(in) :: xtab(:) ! Positions to interpolate
    real(xk), intent(in) :: vtab(:) ! Values to interpolate
    real(xk)             :: v       ! Interpolant
    !
    real(xk)    :: c(size(xtab)), d(size(xtab))
    real(xk)    :: ho, hp, den, scl
    integer(ik) :: ipt, icol, npts, nelem, i
    !
    npts = size(xtab)
    if (size(vtab)/=npts) stop 'math%MathInterpolate - bad input array sizes'
    !
    c    = vtab
    d    = vtab
    ipt  = minloc(abs(x-xtab),dim=1)
    v    = vtab(ipt)
    ipt  = ipt - 1
    tableau_columns: do icol=1,npts-1
      nelem = npts-icol
      update_tableau: do i=1,nelem
        ho   = xtab(i)-x
        hp   = xtab(i+icol)-x
        den  = xtab(i) - xtab(i+icol)
        if (den==0._xk) stop 'math%MathInterpolate - division by zero'
        scl  = (c(i+1) - d(i)) / den
        d(i) = hp*scl
        c(i) = ho*scl
      end do update_tableau
      if (2*ipt<nelem) then
        v   = v + c(ipt+1)
      else
        v   = v + d(ipt)
        ipt = ipt - 1
      end if
    end do tableau_columns
  end function MathInterpolate
end module math
