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
! !
! !  This routine is a subset of LAPACK _DTTS2.
! !
! subroutine m3dp_solve_rr(mf,lpiv,r,x)
!   real(rk), intent(in)     :: mf(:,:)    ! LU factorization;
!                                          ! mf(1,i) is the inverse of the diagonal of U [LAPACK's 1/D]
!                                          ! mf(2,i) is the first super-diagonal of U [LAPACK's DU]
!                                          ! mf(3,i) is the second superdiagonal of U [LAPACK's DU2]
!                                          ! mf(4,i) contains L matrix factors [LAPACK's DL]
!   logical, intent(in)      :: lpiv(:)    ! lpiv(i) is true if row i was interchanged with row i+1
!   real(rk), intent(in)     :: r(:)       ! Right-hand-side vector
!   real(rk), intent(out)    :: x(:)       ! Solution vector
    !
!   character(len=clen), save :: rcsid_tridiagonal_pivoted_m3dp_solve_common = "$Id: tridiagonal_pivoted_m3dp_solve_common.f90,v 1.6 2021/04/26 15:44:44 ps Exp $"
    integer(ik) :: n, i
    !
    n = size(mf,dim=1)
    if (size(mf,dim=2)<4 .or. size(lpiv)/=n .or. size(r,dim=1)/=n .or. size(x,dim=1)/=n) then
      stop 'tridiagonal_pivoted%m3dp_solve - bad array dimensions'
    end if
    !
    !  Solve L*x' = r.
    !
    x(1) = r(1)
    solve_l: do i=1,n-1
      if (.not.lpiv(i)) then
        x(i+1) = r(i+1) - mf(i,4)*x(i)    ! mf(4,i) is the sub-diagonal in column i
      else ! lpiv(i)==.true.
        x(i+1) = x(i) - mf(i,4)*r(i+1)    ! r(i+1) here is actually the solution for the i-th variable, which is isolated
        x(i)   = r(i+1)                   ! 
      end if
    end do solve_l
    !
    !  Solve U*x = x', where U is upper tri-diagonal
    !
    x(n) = x(n) * mf(n,1)     ! mf(1,n) is the (inverse of the) last value on the diagonal of U
    if (n>1) then
      x(n-1) = (x(n-1)-mf(n-1,2)*x(n)) * mf(n-1,1)
    end if
    backsubstitute: do i=n-2,1,-1
      x(i) = (x(i)-mf(i,2)*x(i+1)-mf(i,3)*x(i+2)) * mf(i,1)
    end do backsubstitute
! end subroutine m3dp_solve_rr
