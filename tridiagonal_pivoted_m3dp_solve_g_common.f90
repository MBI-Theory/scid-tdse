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
! subroutine m3dp_solve_rr_g(mf,lpiv,r,x)
!   real(rk), intent(in)     :: mf(:,:)    ! LU factorization;
!                                          ! mf(1,i) is the inverse of the diagonal of U [LAPACK's 1/D]
!                                          ! mf(2,i) is the first super-diagonal of U [LAPACK's DU]
!                                          ! mf(3,i) is the second superdiagonal of U [LAPACK's DU2]
!                                          ! mf(4,i) contains L matrix factors [LAPACK's DL]
!   logical, intent(in)      :: lpiv(:)    ! lpiv(i) is true if row i was interchanged with row i+1
!   real(rk), intent(in)     :: r(:,:)     ! Right-hand-side vectors
!   real(rk), intent(out)    :: x(:,:)     ! Solution vectors
    !
!   character(len=clen), save :: rcsid_tridiagonal_pivoted_m3dp_solve_g_common = "$Id: tridiagonal_pivoted_m3dp_solve_g_common.f90,v 1.7 2021/04/26 15:44:44 ps Exp ps $"
    integer(ik) :: n, i, nrhs, ir
    !
    n    = size(mf,dim=1)
    nrhs = size(r,dim=2)
    if (size(mf,dim=2)<4 .or. size(lpiv)/=n .or. size(r,dim=1)/=n .or. size(x,dim=1)/=n .or. size(x,dim=2)/=nrhs) then
      write (out,"(/'tridiagonal_pivoted%m3dp_solve_g: Bad array dimensions in pivoted solver'/)")
      write (out,"('  mf: ',2i8)") ubound(mf)
      write (out,"('lpiv: ',2i8)") ubound(lpiv)
      write (out,"('   r: ',2i8)") ubound(r)
      write (out,"('   x: ',2i8)") ubound(x)
      call flush_wrapper(out)
      stop 'tridiagonal_pivoted%m3dp_solve_g - bad array dimensions'
    end if
    !
    !  Solve L*x' = r.
    !
    solve_l_rhs: do ir=1,nrhs
      x(1,ir) = r(1,ir)
      solve_l: do i=1,n-1
        if (.not.lpiv(i)) then
          x(i+1,ir) = r(i+1,ir) - mf(i,4)*x(i,ir)    ! mf(4,i) is the sub-diagonal in column i
        else ! lpiv(i)==.true.
          x(i+1,ir) = x(i,ir) - mf(i,4)*r(i+1,ir)    ! r(i+1) here is actually the solution for the i-th variable, which is isolated
          x(i,ir)   = r(i+1,ir)                      ! 
        end if
      end do solve_l
    end do solve_l_rhs
    !
    solve_u_rhs: do ir=1,nrhs
      !
      !  Solve U*x = x', where U is upper tri-diagonal
      !
      x(n,ir) = x(n,ir) * mf(n,1)     ! mf(1,n) is the (inverse of the) last value on the diagonal of U
      if (n>1) then
        x(n-1,ir) = (x(n-1,ir)-mf(n-1,2)*x(n,ir)) * mf(n-1,1)
      end if
      backsubstitute: do i=n-2,1,-1
        x(i,ir) = (x(i,ir)-mf(i,2)*x(i+1,ir)-mf(i,3)*x(i+2,ir)) * mf(i,1)
      end do backsubstitute
    end do solve_u_rhs
! end subroutine m3dp_solve_rr_g
