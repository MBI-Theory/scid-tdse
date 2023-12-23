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
! subroutine m3u_solve_rr_g(mf,r,x)
!   real(rk), intent(in)  :: mf(:,:) ! Factorized tridiagonal matrix
!   real(rk), intent(in)  :: r (:,:) ! Right-hand sizes
!   real(rk), intent(out) :: x (:,:) ! Solution vectors
    !
!   character(len=clen), save :: rcsid_tridiagonal_tools_m3u_solve_g_common = "$Id: tridiagonal_tools_m3u_solve_g_common.f90,v 1.5 2021/04/26 15:44:44 ps Exp $"
    !
    integer(ik) :: i, sz, ir
    !
    sz = size(mf,dim=1)
    if (size(mf,dim=2)<3 .or. size(r,dim=1)/=sz .or. size(x,dim=1)/=sz .or. size(x,dim=2)/=size(r,dim=2)) then
      stop 'tridiagonal_tools%m3u_solve_g_common - bad input sizes'
    end if
    ! A decent compiler will rearrange this loop nest as needed ....
    scan_rhs: do ir=1,size(r,dim=2)
      x(1,ir) = mf(1,1)*r(1,ir)
      transform_rhs: do i=2,sz
        x(i,ir) = mf(i,1)*r(i,ir) + mf(i,2)*x(i-1,ir)
      end do transform_rhs
      !
      backsubstitute: do i=sz-1,1,-1
        x(i,ir) = x(i,ir) - mf(i,3)*x(i+1,ir)
      end do backsubstitute
    end do scan_rhs
! end subroutine m3u_solve_rr_g
