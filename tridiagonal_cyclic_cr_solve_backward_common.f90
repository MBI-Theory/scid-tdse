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
! subroutine cr_solve_backward_step(f,r,xr,x)
!   real(rk), intent(in)  :: f(:,:) ! Reduction data from cr_reduce_step
!   real(rk), intent(in)  :: r (:)  ! RHS
!   real(rk), intent(in)  :: xr(:)  ! Reduced solution vector
!   real(rk), intent(out) :: x (:)  ! Full Solution vector
    !
    integer(ik) :: nm, nf, nr
    !
    nm = size(r)       ! Total size of the linear problem
    nf = size(f,dim=1) ! Number of eliminated degrees of freedom
    nr = size(xr)      ! Reduced number of the degrees of freedom
    !
    if (size(f,dim=2)<8) stop 'tridiagonal_cyclic%cr_solver_backward_step - bad argument dimensions (1)'
    if (nm/=nf+nr .or. (nf/=nr .and. nf/=nr+1)) stop 'tridiagonal_cyclic%cr_solver_bakward_step - bad argument dimensions (2)'
    if (size(x)/=size(r)) stop 'tridiagonal_cyclic%cr_solver_backward_step - bad argument dimensions (3)'
    !
    x(2:2*nr:2)   = xr(1:nr)
    x(1:2*nf-1:2) = f(1:nf,1) * r(1:2*nf-1:2)
    x(1:2*nr-1:2) = x(1:2*nr-1:2) + f(1:nr,2)*xr(1:nr)
    x(3:2*nf-1:2) = x(3:2*nf-1:2) + f(1:nf-1,3)*xr(1:nf-1)
! end subroutine cr_solve_backward_step
