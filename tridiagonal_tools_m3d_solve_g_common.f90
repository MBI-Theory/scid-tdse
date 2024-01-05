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
! subroutine m3d_solve_rr_g(m,mf,lpiv,r,x,scr)
!   real(rk), intent(in)           :: m (:,:)  ! Tridiagonal matrix
!   real(rk), intent(in)           :: mf(:,:)  ! Factorized tridiagonal matrix
!   logical, intent(in)            :: lpiv(:)  ! Pivot list
!   real(rk), intent(in)           :: r(:,:)   ! Right-hand sides
!   real(rk), intent(out)          :: x(:,:)   ! Solutions of the linear systems
!   real(rk), intent(out)          :: scr(:,:) ! Scratch for iterative refinement
    !
!   character(len=clen), save :: rcsid_tridiagonal_tools_m3d_solve_g_common = "$Id: tridiagonal_tools_m3d_solve_g_common.f90,v 1.8 2021/04/26 15:44:44 ps Exp $"
    ! 
    integer(ik) :: ipass
    integer(ik) :: nrhs
    !
    if (size(mf,dim=2)/=m3d_dc_size .or. size(scr,dim=2)<m3d_sc_size*size(r,dim=2)) then
      stop 'tridiagonal_tools%tridiagonal_tools_m3d_solve_g: Bug: bad size of the decomposion or scratch arrays'
    end if
    !
    select case (m3d_solver)
      case default
        write (out,"('tridiagonal_tools%m3d_solve_g: Unrecognized solver ',a)") trim(m3d_solver)
        stop 'tridiagonal_tools%m3d_solve_g: Bad solver choice'
      case ('unpivoted')
        call m3u_solve(mf,r,x)
      case ('cyclic')
        call m3c_solve(mf,r,x,scr)
      case ('pivoted')
        call m3dp_solve(mf,lpiv,r,x)
      case ('refined')
        nrhs = size(r,dim=2)
        if (size(scr,dim=1)/=size(r,dim=1) .or. size(scr,dim=2)<2*nrhs) then
          write (out,"(/'tridiagonal_tools%m3d_solve_g: Bad scratch array for iterative refinement'/)")
          write (out,"('need at least: ',2i8)") size(r,dim=1), 2*nrhs
          write (out,"('          got: ',2i8)") ubound(scr)
          call flush_wrapper(out)
          stop 'tridiagonal_tools%m3d_solve_g: Bad scratch array for iterative refinement'
        end if
        call m3dp_solve(mf,lpiv,r,x)
        refine_iterations: do ipass=1,m3d_iterations
          call m3d_multiply(m,x,scr(:,1:nrhs))
          scr(:,1:nrhs) = scr(:,1:nrhs) - r(:,1:nrhs)
          call m3dp_solve(mf,lpiv,scr(:,1:nrhs),scr(:,nrhs+1:2*nrhs))
          x(:,1:nrhs) = x(:,1:nrhs) - scr(:,nrhs+1:2*nrhs)
        end do refine_iterations
    end select
! end subroutine m3d_solve_rr_g
