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
! subroutine m3d_decompose_r(m,mf,lpiv,fail)
!   real(rk), intent(in)           :: m (:,:) ! Tridiagonal matrix, assumed to be diagonally dominant
!   real(rk), intent(out)          :: mf(:,:) ! Factorized tridiagonal matrix
!   logical, intent(out)           :: lpiv(:) ! Pivot list
!   logical, intent(out), optional :: fail    ! Set to .true. if decomposition fails; 
!                                             ! if fail is absent, abort on decomposition failure.
!   !
!   character(len=clen), save :: rcsid_tridiagonal_tools_m3d_decompose_x_common = "$Id: tridiagonal_tools_m3d_decompose_x_common.f90,v 1.5 2021/04/26 15:44:44 ps Exp $"
    !
    if (size(mf,dim=2)/=m3d_dc_size) then
      stop 'tridiagonal_tools%m3d_decompose_x: Bug - decomposion array wrong size'
    end if
    !
    select case (m3d_solver)
      case default
        write (out,"('tridiagonal_tools%m3d_decompose: Unrecognized solver ',a)") trim(m3d_solver)
        stop 'tridiagonal_tools%m3d_decompose: Bad solver choice'
      case ('unpivoted','cyclic')
        ! This routine is never on the critical path; there is no point in calling the CR routine?
        if (present(fail)) then
          call m3u_decompose_x(m,mf,fail)
        else
          call m3u_decompose_x(m,mf)
        end if
      case ('pivoted','refined')
        if (present(fail)) then
          call m3dp_decompose_x(m,mf,lpiv,fail)
        else
          call m3dp_decompose_x(m,mf,lpiv)
        end if
    end select
    !
! end subroutine m3d_decompose_r
