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
  !  Single backsubstitution step
  !
! subroutine cr_solve_backward_step_g(f,r,xr,x)
!   real(rk), intent(in)  :: f (:,:) ! Reduction data from cr_reduce_step
!   real(rk), intent(in)  :: r (:,:) ! RHS set
!   real(rk), intent(in)  :: xr(:,:) ! Reduced solution vectors
!   real(rk), intent(out) :: x (:,:) ! Full Solution vectors
    !
    integer(ik) :: nm, nf, nr, nrhs, j
    !
    nm = size(r,dim=1) ! Total size of the linear problem
    nf = size(f,dim=1) ! Number of eliminated degrees of freedom
    nr = size(xr,dim=1)! Reduced number of the degrees of freedom
    nrhs = size(r,dim=2) ! Number of right-hand sides
    !
    if (size(f,dim=2)<8) stop 'tridiagonal_cyclic%cr_solver_backward_step_g - bad argument dimensions (1)'
    if (nm/=nf+nr .or. (nf/=nr .and. nf/=nr+1)) stop 'tridiagonal_cyclic%cr_solver_bakward_step_g - bad argument dimensions (2)'
    if (size(x,dim=1)/=nm) stop 'tridiagonal_cyclic%cr_solver_backward_step_g - bad argument dimensions (3)'
    if (size(x,dim=2)/=nrhs .or. size(xr,dim=2)/=nrhs) &
       stop 'tridiagonal_cyclic%cr_solver_backward_step_g - bad argument dimensions (4)'
    !
    !  gfortran (at least up ro 12.2) generates bad code for the vector version,
    !  using spread() intrinsic, which involves creation of completely unnecessary 
    !  array temporaries. We'll use the more verbose but equivalent forall instead.
    !
    ! x(2:2*nr:2,:)   = xr(1:nr,:)
    ! x(1:2*nf-1:2,:) = spread(f(1:nf,1),2,nrhs) * r(1:2*nf-1:2,:)
    ! x(1:2*nr-1:2,:) = x(1:2*nr-1:2,:) + spread(f(1:nr,2),2,nrhs)*xr(1:nr,:)
    ! x(3:2*nf-1:2,:) = x(3:2*nf-1:2,:) + spread(f(1:nf-1,3),2,nrhs)*xr(1:nf-1,:)
    !
    rhs_solve: forall (j=1:nrhs) 
      x(2:2*nr:2,j)   = xr(1:nr,j)
      x(1:2*nf-1:2,j) = f(1:nf,1) * r(1:2*nf-1:2,j)
      x(1:2*nr-1:2,j) = x(1:2*nr-1:2,j) + f(1:nr,2)  *xr(1:nr,j)
      x(3:2*nf-1:2,j) = x(3:2*nf-1:2,j) + f(1:nf-1,3)*xr(1:nf-1,j)
    end forall rhs_solve
! end subroutine cr_solve_backward_step_g
