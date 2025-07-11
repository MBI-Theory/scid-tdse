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
!  Type and precision definitions
!  Theoretically, changing the value of rk and/or ik here should be sufficient to 
!  switch the entire program to the desired accuracy.
!
module accuracy
  implicit none
  private
  public ik, rk, xk
  public input, out
  public clen
  public srk, drk, lrk, lik
!*nq  public qrk
!*qd  public qrk
  public rk_bytes, ik_bytes, xk_bytes
  public flush_wrapper
  public isnan_wrapper
  public rcsid_accuracy
  !
  integer, parameter :: ik     = selected_int_kind(8)       ! "Normal" integers
! integer, parameter :: ik     = selected_int_kind(15)      ! "Long" integers
! integer, parameter :: rk     = selected_real_kind(6,17)   ! Low-precision
  integer, parameter :: rk     = selected_real_kind(14,17)  ! "Normal" reals and complex (complexi? :-)
! integer, parameter :: rk     = selected_real_kind(28,50)  ! quad-precision
  !
  !  ifort 15.0.1 mistakenly issues a warning for the following line, complaining about 
  !  the mismatch in the type/kind of precision(1._rk) and 14. Disregard the warning.
  integer, parameter :: xk     = selected_real_kind(max(precision(1._rk),14),max(range(1._rk),17))
                                                            ! Real kind used for manipulating accuracy-sensitive
                                                            ! quantites. This includes time intervals
                                                            ! and field. Will have at least the range and
                                                            ! accuracy of the (rk) real kind, but not less
                                                            ! than the range and accuracy of double-precision type
                                                            ! This real kind should never appear on time-critical
                                                            ! path.
  integer(ik), parameter :: input  = 5                      ! Standard input I/O channel
  integer(ik), parameter :: out    = 6                      ! Output I/O channel
  integer, parameter     :: clen   = 255                    ! Standard character length; enough for most
                                                            ! keywords and file names
  !
  character(len=clen), save :: rcsid_accuracy = "$Id: accuracy.f90,v 1.47 2024/04/23 14:29:50 ps Exp $"
  !
  !  System kinds; should only be used where we interface to external libraries
  !
  integer, parameter :: srk    = kind(1.)                   ! System default real kind
  integer, parameter :: drk    = kind(1.d0)                 ! System double-precision real kind
!*nq  integer, parameter :: qrk    = kind(1.d0)             ! No quad-precision real kind; use double dummy
!*qd  integer, parameter :: qrk    = kind(1.q0)             ! System quad-precision real kind
  !
  !  (lrk) is the largest real kind where optimized lapack library is available. It is used to
  !        calculate guess for the atomic solutions; these are then refined using iterative methods.
  !
  integer, parameter :: lrk    = selected_real_kind(min(precision(1._rk),precision(1._drk)),min(range(1._rk),range(1._drk)))
! integer, parameter :: lrk    = qrk                        ! Use quad precision for LAPACK
  integer, parameter :: lik    = kind(1)                    ! Use default integer kind for LAPACK
! integer, parameter :: lik    = selected_int_kind(15)      ! Use "long" integer kind for LAPACK
  !
  interface isnan_wrapper
    module procedure isnan_real
    module procedure isnan_complex
  end interface isnan_wrapper
  !
  contains

  integer(ik) function rk_bytes ()
    rk_bytes = real_bytes(radix(1._rk) ,digits(1._rk ),maxexponent(1._rk )-minexponent(1._rk ))
  end function rk_bytes
 
  integer(ik) function xk_bytes ()
    xk_bytes = real_bytes(radix(1._xk) ,digits(1._xk ),maxexponent(1._xk )-minexponent(1._xk ))
  end function xk_bytes
 
  integer(ik) function ik_bytes ()
    ik_bytes = int_bytes(radix(1_ik ),digits(1_ik ))
  end function ik_bytes

  integer(ik) function int_bytes(radix,digits)
    integer, intent(in) :: radix, digits
    !
    real(rk) :: bits
    !
    bits = real(digits,kind=rk)*(log(real(radix,kind=rk))/log(2._rk))
    int_bytes = ceiling((1+bits)/8._rk)
  end function int_bytes

  integer(ik) function real_bytes(radix,digits,exp_range)
    integer, intent(in) :: radix, digits, exp_range
    !
    real(rk) :: exp_bits, mant_bits
    !
    exp_bits   = log(real(exp_range,kind=rk))/log(2._rk)
    mant_bits  = real(digits,kind=rk)*(log(real(radix,kind=rk))/log(2._rk))
    real_bytes = ceiling((exp_bits+mant_bits)/8._rk)
  end function real_bytes

  subroutine flush_wrapper(unit)
    integer(ik), intent(in) :: unit
    !
    flush (unit)  ! Fortran-2003 form
    ! intrinsic flush
    ! call flush(unit) ! Common extension; uncomment if Fortran-2003 statement is not available
  end subroutine flush_wrapper

  elemental function isnan_real(x) result(v)
    ! Depending on the compiler, it may be necessary to uncomment one of the lines below
    ! use IEEE_ARITHMETIC , only: isnan => ieee_is_nan ! This may be very expensive - it keeps switching FPU mode flags.
      intrinsic            :: isnan                    ! This works for relatively recent gfortran and ifort versiobns
    ! logical, intrinsic   :: isnan                    ! ... or somewhat older versions
    ! logical, external    :: isnan
    real(rk), intent(in) :: x
    logical              :: v
    !
    v = isnan(x)
  end function isnan_real

  elemental function isnan_complex(x) result(v)
    complex(rk), intent(in) :: x
    logical                 :: v
    !
    v = isnan_real(real(x,kind=rk)) .or. isnan_real(aimag(x))
  end function isnan_complex

end module accuracy
