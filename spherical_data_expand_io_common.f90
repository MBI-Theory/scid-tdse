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
! subroutine sd_expand_implicit_operator_r(op,mop,mopf,mopp,block)
!   real(rk), intent(in)  :: op   (:,:) ! Tridiagonal operator matrix
!   real(rk), intent(in)  :: mop  (:,:) ! Tridiagonal modifier operator
!   real(rk), intent(in)  :: mopf (:,:) ! Factored modifier operator
!   logical, intent(in)   :: mopp (  :) ! Pivot list for the factored modifier
!   real(rk), intent(out) :: block(:,:) ! Explicit operator matrix
    !
!   character(len=clen), save :: rcsid_spherical_data_expand_io_common = "$Id: spherical_data_expand_io_common.f90,v 1.7 2021/04/26 15:44:44 ps Exp $"
    integer(ik) :: sz, icol
!   real(rk)    :: rhs(size(block,dim=1))
!   real(rk)    :: tmp(size(block,dim=1))
!   real(rk)    :: scr(size(block,dim=1),m3d_sc_size)
    !
    sz = size(block,dim=1)
    !
    !  Sanity check
    !
    if (size(op,dim=2)<3 .or. size(mop,dim=2)<3 .or. size(mopf,dim=2)<m3d_dc_size .or. &
        size(op,dim=1)/=sz .or. size(mop,dim=1)/=sz .or. size(mopf,dim=1)/=sz .or. &
        size(mopp)/=sz .or. size(block,dim=1)/=sz) then
      stop 'spherical_data%sd_expand_implicit_operator - dimensions mismatch'
    end if
    !
    rhs(:) = 0
    scan_columns: do icol=1,sz
      rhs(icol) = 1._rk
      call m3d_multiply(op,rhs,tmp)
      call m3d_solve(mop,mopf,mopp,tmp,block(:,icol),scr)
      rhs(icol) = 0._rk
    end do scan_columns
! end subroutine sd_expand_implicit_operator_r
