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
  !  Single forward solution step
  !
! subroutine cr_solve_forward_step_g(f,r,rr)
!   real(rk), intent(in)  :: f (:,:)  ! Reduction data from cr_reduce_step
!   real(rk), intent(in)  :: r (:,:)  ! Full RHS set
!   real(rk), intent(out) :: rr(:,:)  ! Reduced RHS set on this step
    !
    integer(ik) :: nm, nf, nr, nrhs, j
    !
    nm = size(r,dim=1)    ! Total size of the linear problem
    nf = size(f,dim=1)    ! Number of eliminated degrees of freedom
    nr = size(rr,dim=1)   ! Reduced number of the degrees of freedom
    nrhs = size(r,dim=2)  ! Number of right-hand sides
    !
    if (size(f,dim=2)<8) stop 'tridiagonal_cyclic%cr_solver_forward_step_g - bad argument dimensions (1)'
    if (nm/=nf+nr .or. (nf/=nr .and. nf/=nr+1)) stop 'tridiagonal_cyclic%cr_solver_forward_step_g - bad argument dimensions (2)'
    if (size(rr,dim=2)/=nrhs) stop 'tridiagonal_cyclic%cr_solver_forward_step_g - bad argument dimensions (3)'
    !
    !  gfortran, as of version 12.2, creates explicit array temporaries for the spread()
    !  intrinsic. rewrite with forall
    !
    ! rr(1:nr,:)   = r(2:2*nr:2,:) + spread(f(1:nr,4),2,nrhs)*r(1:2*nr-1:2,:) 
    ! rr(1:nf-1,:) = rr(1:nf-1,:) + spread(f(1:nf-1,5),2,nrhs)*r(3:2*nf-1:2,:)
    rhs_loop: forall (j=1:nrhs) 
      rr(1:nr,j)   = r(2:2*nr:2,j) + f(1:nr,4)  *r(1:2*nr-1:2,j) 
      rr(1:nf-1,j) = rr(1:nf-1,j)  + f(1:nf-1,5)*r(3:2*nf-1:2,j)
    end forall rhs_loop
! end subroutine cr_solve_forward_step_g
