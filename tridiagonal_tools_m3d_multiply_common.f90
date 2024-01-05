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
! subroutine m3d_multiply_rr(m,v,mv)
!   real(rk), intent(in)  :: m(:,:) ! Tri-diagonal matrix
!   real(rk), intent(in)  :: v(:)   ! Vector 
!   real(rk), intent(out) :: mv(:)  ! Vector 
    !
!   character(len=clen), save :: rcsid_tridiagonal_tools_m3d_multiply_common = "$Id: tridiagonal_tools_m3d_multiply_common.f90,v 1.5 2021/04/26 15:44:44 ps Exp $"
    !
    integer(ik) :: sz
    !
    sz = size(v)
    if (size(m,dim=2)<3 .or. size(m,dim=1)/=sz .or. size(mv)/=sz) then
      stop 'tridiagonal_tools%m3d_multiply_common - bad input sizes'
    end if
    mv(1:sz  ) =              m(:    ,1)*v(1:sz  ) 
    mv(2:sz  ) = mv(2:sz  ) + m(:sz-1,2)*v(1:sz-1)
    mv(1:sz-1) = mv(1:sz-1) + m(:sz-1,3)*v(2:sz  )
! end subroutine m3d_multiply_rr
