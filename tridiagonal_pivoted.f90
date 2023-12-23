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
!  Solution of tridiagonal linear systems with partial pivoting and (optionally) iterative
!  improvement. The code is adapted from LAPACK _GTSVX routines.
!
module tridiagonal_pivoted
  use accuracy
  use constants
  use timer
  implicit none
  private
  public m3dp_decompose, m3dp_solve
  public m3dp_decompose_x, m3dp_solve_x
  public rcsid_tridiagonal_pivoted
  !
  character(len=clen), save :: rcsid_tridiagonal_pivoted = "$Id: tridiagonal_pivoted.f90,v 1.6 2021/04/26 15:44:44 ps Exp $"
  !
  interface m3dp_decompose
    module procedure m3dp_decompose_r
    module procedure m3dp_decompose_c
  end interface m3dp_decompose
  interface m3dp_decompose_x
    module procedure m3dp_decompose_xr
  end interface m3dp_decompose_x
  !
  interface m3dp_solve
    module procedure m3dp_solve_rr
    module procedure m3dp_solve_rc
    module procedure m3dp_solve_cc
    module procedure m3dp_solve_rr_g
    module procedure m3dp_solve_rc_g
    module procedure m3dp_solve_cc_g
  end interface m3dp_solve
  interface m3dp_solve_x
    module procedure m3dp_solve_xrr_g
  end interface m3dp_solve_x
  !
  contains
  !
  !  m3dp_decompose() is a clone of LAPACK _GTTRF(). It computes LU factorization of
  !  a tridiagonal matrix using elimination with partial pivoting and row interchanges.
  !
  subroutine m3dp_decompose_r(m,mf,lpiv,fail)
    real(rk), intent(in)           :: m (:,:) ! Tridiagonal matrix to factorize;
                                              ! m(i,1) is main diagonal (i,i)
                                              ! m(i,2) is sub-diagonal (i+1,i)
                                              ! m(i,3) is super-diagonal (i,i+1)
    real(rk), intent(out)          :: mf(:,:) ! LU factorization;
                                              ! mf(i,1) is the inverse of the diagonal of U [LAPACK's 1/D]
                                              ! mf(i,2) is the first super-diagonal of U [LAPACK's DU]
                                              ! mf(i,3) is the second superdiagonal of U [LAPACK's DU2]
                                              ! mf(i,4) contains L matrix factors [LAPACK's DL]
    logical, intent(out)           :: lpiv(:) ! lpiv(i) is true if row i was interchanged with row i+1
    logical, intent(out), optional :: fail    ! Set to true if decomposition fails
    !
    real(rk)      :: fact, temp         ! Must match the type of mf()
    !
    include "tridiagonal_pivoted_m3dp_decompose_common.f90"
  end subroutine m3dp_decompose_r
  !
  subroutine m3dp_decompose_xr(m,mf,lpiv,fail)
    real(xk), intent(in)           :: m (:,:) ! Tridiagonal matrix to factorize;
                                              ! m(i,1) is main diagonal (i,i)
                                              ! m(i,2) is sub-diagonal (i+1,i)
                                              ! m(i,3) is super-diagonal (i,i+1)
    real(xk), intent(out)          :: mf(:,:) ! LU factorization;
                                              ! mf(i,1) is the inverse of the diagonal of U [LAPACK's 1/D]
                                              ! mf(i,2) is the first super-diagonal of U [LAPACK's DU]
                                              ! mf(i,3) is the second superdiagonal of U [LAPACK's DU2]
                                              ! mf(i,4) contains L matrix factors [LAPACK's DL]
    logical, intent(out)           :: lpiv(:) ! lpiv(i) is true if row i was interchanged with row i+1
    logical, intent(out), optional :: fail    ! Set to true if decomposition fails
    !
    real(xk)      :: fact, temp         ! Must match the type of mf()
    !
    include "tridiagonal_pivoted_m3dp_decompose_common.f90"
  end subroutine m3dp_decompose_xr
  !
  subroutine m3dp_decompose_c(m,mf,lpiv,fail)
    complex(rk), intent(in)        :: m (:,:) ! Tridiagonal matrix to factorize;
                                              ! m(i,1) is main diagonal (i,i)
                                              ! m(i,2) is sub-diagonal (i+1,i)
                                              ! m(i,3) is super-diagonal (i,i+1)
    complex(rk), intent(out)       :: mf(:,:) ! LU factorization;
                                              ! mf(i,1) is the inverse of the diagonal of U [LAPACK's 1/D]
                                              ! mf(i,2) is the first super-diagonal of U [LAPACK's DU]
                                              ! mf(i,3) is the second superdiagonal of U [LAPACK's DU2]
                                              ! mf(i,4) contains L matrix factors [LAPACK's DL]
    logical, intent(out)           :: lpiv(:) ! lpiv(i) is true if row i was interchanged with row i+1
    logical, intent(out), optional :: fail    ! Set to true if decomposition fails
    !
    complex(rk)      :: fact, temp         ! Must match the type of mf()
    !
    include "tridiagonal_pivoted_m3dp_decompose_common.f90"
  end subroutine m3dp_decompose_c
  !
  !  This routine is a subset of LAPACK _DTTS2.
  !
  subroutine m3dp_solve_rr(mf,lpiv,r,x)
    real(rk), intent(in)     :: mf(:,:)    ! LU factorization;
                                           ! mf(i,1) is the inverse of the diagonal of U [LAPACK's 1/D]
                                           ! mf(i,2) is the first super-diagonal of U [LAPACK's DU]
                                           ! mf(i,3) is the second superdiagonal of U [LAPACK's DU2]
                                           ! mf(i,4) contains L matrix factors [LAPACK's DL]
    logical, intent(in)      :: lpiv(:)    ! lpiv(i) is true if row i was interchanged with row i+1
    real(rk), intent(in)     :: r(:)       ! Right-hand-side vector
    real(rk), intent(out)    :: x(:)       ! Solution vector
    !
    include "tridiagonal_pivoted_m3dp_solve_common.f90"
  end subroutine m3dp_solve_rr
  !
  subroutine m3dp_solve_rc(mf,lpiv,r,x)
    real(rk), intent(in)     :: mf(:,:)    ! LU factorization;
                                           ! mf(i,1) is the inverse of the diagonal of U [LAPACK's 1/D]
                                           ! mf(i,2) is the first super-diagonal of U [LAPACK's DU]
                                           ! mf(i,3) is the second superdiagonal of U [LAPACK's DU2]
                                           ! mf(i,4) contains L matrix factors [LAPACK's DL]
    logical, intent(in)      :: lpiv(:)    ! lpiv(i) is true if row i was interchanged with row i+1
    complex(rk), intent(in)  :: r(:)       ! Right-hand-side vector
    complex(rk), intent(out) :: x(:)       ! Solution vector
    !
    include "tridiagonal_pivoted_m3dp_solve_common.f90"
  end subroutine m3dp_solve_rc
  !
  subroutine m3dp_solve_cc(mf,lpiv,r,x)
    complex(rk), intent(in)  :: mf(:,:)    ! LU factorization;
                                           ! mf(i,1) is the inverse of the diagonal of U [LAPACK's 1/D]
                                           ! mf(i,2) is the first super-diagonal of U [LAPACK's DU]
                                           ! mf(i,3) is the second superdiagonal of U [LAPACK's DU2]
                                           ! mf(i,4) contains L matrix factors [LAPACK's DL]
    logical, intent(in)      :: lpiv(:)    ! lpiv(i) is true if row i was interchanged with row i+1
    complex(rk), intent(in)  :: r(:)       ! Right-hand-side vector
    complex(rk), intent(out) :: x(:)       ! Solution vector
    !
    include "tridiagonal_pivoted_m3dp_solve_common.f90"
  end subroutine m3dp_solve_cc
  !
  subroutine m3dp_solve_rr_g(mf,lpiv,r,x)
    real(rk), intent(in)     :: mf(:,:)    ! LU factorization;
                                           ! mf(i,1) is the inverse of the diagonal of U [LAPACK's 1/D]
                                           ! mf(i,2) is the first super-diagonal of U [LAPACK's DU]
                                           ! mf(i,3) is the second superdiagonal of U [LAPACK's DU2]
                                           ! mf(i,4) contains L matrix factors [LAPACK's DL]
    logical, intent(in)      :: lpiv(:)    ! lpiv(i) is true if row i was interchanged with row i+1
    real(rk), intent(in)     :: r(:,:)     ! Right-hand-side vectors
    real(rk), intent(out)    :: x(:,:)     ! Solution vectors
    !
    include "tridiagonal_pivoted_m3dp_solve_g_common.f90"
  end subroutine m3dp_solve_rr_g
  !
  subroutine m3dp_solve_xrr_g(mf,lpiv,r,x)
    real(xk), intent(in)     :: mf(:,:)    ! LU factorization;
                                           ! mf(i,1) is the inverse of the diagonal of U [LAPACK's 1/D]
                                           ! mf(i,2) is the first super-diagonal of U [LAPACK's DU]
                                           ! mf(i,3) is the second superdiagonal of U [LAPACK's DU2]
                                           ! mf(i,4) contains L matrix factors [LAPACK's DL]
    logical, intent(in)      :: lpiv(:)    ! lpiv(i) is true if row i was interchanged with row i+1
    real(xk), intent(in)     :: r(:,:)     ! Right-hand-side vectors
    real(xk), intent(out)    :: x(:,:)     ! Solution vectors
    !
    include "tridiagonal_pivoted_m3dp_solve_g_common.f90"
  end subroutine m3dp_solve_xrr_g
  !
  subroutine m3dp_solve_rc_g(mf,lpiv,r,x)
    real(rk), intent(in)     :: mf(:,:)    ! LU factorization;
                                           ! mf(i,1) is the inverse of the diagonal of U [LAPACK's 1/D]
                                           ! mf(i,2) is the first super-diagonal of U [LAPACK's DU]
                                           ! mf(i,3) is the second superdiagonal of U [LAPACK's DU2]
                                           ! mf(i,4) contains L matrix factors [LAPACK's DL]
    logical, intent(in)      :: lpiv(:)    ! lpiv(i) is true if row i was interchanged with row i+1
    complex(rk), intent(in)  :: r(:,:)     ! Right-hand-side vectors
    complex(rk), intent(out) :: x(:,:)     ! Solution vectors
    !
    include "tridiagonal_pivoted_m3dp_solve_g_common.f90"
  end subroutine m3dp_solve_rc_g
  !
  subroutine m3dp_solve_cc_g(mf,lpiv,r,x)
    complex(rk), intent(in)  :: mf(:,:)    ! LU factorization;
                                           ! mf(i,1) is the inverse of the diagonal of U [LAPACK's 1/D]
                                           ! mf(i,2) is the first super-diagonal of U [LAPACK's DU]
                                           ! mf(i,3) is the second superdiagonal of U [LAPACK's DU2]
                                           ! mf(i,4) contains L matrix factors [LAPACK's DL]
    logical, intent(in)      :: lpiv(:)    ! lpiv(i) is true if row i was interchanged with row i+1
    complex(rk), intent(in)  :: r(:,:)     ! Right-hand-side vectors
    complex(rk), intent(out) :: x(:,:)     ! Solution vectors
    !
    include "tridiagonal_pivoted_m3dp_solve_g_common.f90"
  end subroutine m3dp_solve_cc_g
end module tridiagonal_pivoted
