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
! subroutine m3d_left_scale_cr(s,m,sm)
!   complex(rk), intent(in)  :: s   (:) ! Diagonal matrix
!   real(rk), intent(in)     :: m (:,:) ! Tridiagonal matrix
!   complex(rk), intent(out) :: sm(:,:) ! Tridiagonal s . m
    !
!   character(len=clen), save :: rcsid_tridiagonal_tools_m3d_left_scale_common = "$Id: tridiagonal_tools_m3d_left_scale_common.f90,v 1.5 2021/04/26 15:44:44 ps Exp $"
    integer(ik) :: sz
    !
    sz = size(m,dim=1)
    if (size(m,dim=2)<3 .or. size(sm,dim=2)<3 .or. size(s)/=sz .or. size(m,dim=1)/=sz) then
      stop 'tridiagonal_tools%m3d_left_scale_common - bad input sizes'
    end if
    !
    sm(1:sz  ,1) = s(1:sz  ) * m(1:sz  ,1)
    sm(1:sz-1,2) = s(2:sz  ) * m(1:sz-1,2)
    sm(1:sz-1,3) = s(1:sz-1) * m(1:sz-1,3)
! end subroutine m3d_left_scale_cr
