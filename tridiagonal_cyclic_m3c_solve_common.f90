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
! subroutine m3c_solve(mf,rhs,x,scr)
!   real(rk), intent(in)  :: mf(:,:)  ! Factorization data
!   real(rk), intent(in)  :: rhs(:)   ! Right-hand side
!   real(rk), intent(out) :: x(:)     ! Solution of the linear system
!   real(rk), intent(out) :: scr(:,:) ! Scratch space
    !
    integer(ik) :: nm        ! Current number of variables
    integer(ik) :: nf        ! Number of variables eliminated in this round
    integer(ik) :: nr        ! Number of variables in the reduced system of equations
    integer(ik) :: ip        ! Current position in the reduction buffer
    integer(ik) :: ip0       ! Previous position in the reduction buffer
    integer(ik) :: is        ! Step counter
    integer(ik) :: nsteps    ! Number of recursive steps; can't be more than max_cr_steps
    integer(ik) :: nm_stack(max_cr_steps)
    integer(ik) :: ip_stack(max_cr_steps)
    !
    nm = size(mf,dim=1)
    if (nm<1) stop 'tridiagonal_cyclic%m3c_solve - must have at least one linear equation!'
    if (size(mf,dim=2)<8 .or. size(scr,dim=2)<2 .or. size(rhs)/=nm .or. size(x)/=nm .or. size(scr,dim=1)/=nm) then
      stop 'tridiagonal_cyclic%m3c_solve - bad array dimensions'
    end if
    ip     = 1
    nsteps = 0
    !
    !  First forward substitution step is special: the RHS is in the input array
    !
    nsteps = nsteps + 1
    nm_stack (nsteps) = nm
    ip_stack (nsteps) = ip
    !
    nr  = nm/2 
    nf  = nm - nr
    call cr_solve_forward_step(mf(ip:ip+nf-1,:),rhs,scr(ip:ip+nr-1,1))
    ip0 = ip
    ip  = ip + nf
    nm  = nr
    !
    !  Continue reduction using our running rhs
    !
    forward_reduction: do while(nm>0)
      nsteps = nsteps + 1
      if (nsteps>max_cr_steps) stop 'tridiagonal_cyclic%m3c_solve - stack logic failed!'
      nm_stack (nsteps) = nm
      ip_stack (nsteps) = ip
      !
      nr  = nm/2 
      nf  = nm - nr
      call cr_solve_forward_step(mf(ip:ip+nf-1,:),scr(ip0:ip0+nm-1,1),scr(ip:ip+nr-1,1))
      ip0 = ip
      nm  = nr
      ip  = ip + nf
    end do forward_reduction
    !
    !  Backsubstitution
    !
    backward_substitution: do is=nsteps,2,-1
      nm  = nm_stack (is)
      ip  = ip_stack (is)
      ip0 = ip_stack (is-1)
      nr  = nm/2 
      nf  = nm - nr
      call cr_solve_backward_step(mf(ip:ip+nf-1,:),scr(ip0:ip0+nm-1,1),scr(ip:ip+nr-1,2),scr(ip0:ip0+nm-1,2))
    end do backward_substitution
    !
    !  The final backsubstitution step is different
    !
    nm  = nm_stack(1)
    ip  = 1
    nr  = nm/2
    nf  = nm - nr
    call cr_solve_backward_step(mf(ip:ip+nf-1,:),rhs,scr(ip:ip+nr-1,2),x)
! end subroutine m3c_solve
