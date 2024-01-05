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
!  Solution of tridiagonal linear systems using unpivoted cyclic-reduction algorthm.
!
!  See Gander and Golub, "Cyclic Reduction - History and Applications", Proceedings of
!  the Workshop on Scientific Computing 10-12 March 1997 (Springer, NY), available at
!  http://www.inf.ethz.ch/personal/gander/papers/cyclic.pdf
!
!  We choose bot to provide the *_x (high-precision on low-precision build) version
!  of this solver. Numerical properties of the cyclic-reduction solver are the
!  same as for the unpivoted solver, so there is no point in implementing it for
!  calls not on the critical path.
!
module tridiagonal_cyclic
  use accuracy
  use constants
  use timer
  implicit none
  private
  public m3c_decompose, m3c_solve
  public rcsid_tridiagonal_cyclic
  !
  character(len=clen), save :: rcsid_tridiagonal_cyclic = "$Id: tridiagonal_cyclic.f90,v 1.2 2021/04/26 15:44:44 ps Exp $"
  !
  integer(ik), parameter :: max_cr_steps = 64 ! Enough to exhaust address space on a 64-bit system several times over
  !
  !  Public interfaces
  !
  interface m3c_decompose
    module procedure m3c_decompose_r
    module procedure m3c_decompose_c
  end interface m3c_decompose
  !
  interface m3c_solve
    module procedure m3c_solve_rr
    module procedure m3c_solve_rc
    module procedure m3c_solve_cc
    module procedure m3c_solve_rr_g
    module procedure m3c_solve_rc_g
    module procedure m3c_solve_cc_g
  end interface m3c_solve
  !
  !  Internal interfaces
  !
  interface cr_reduce_step
    module procedure cr_reduce_step_r
    module procedure cr_reduce_step_c
  end interface cr_reduce_step
  !
  interface cr_solve_forward_step
    module procedure cr_solve_forward_step_rr
    module procedure cr_solve_forward_step_rc
    module procedure cr_solve_forward_step_cc
    module procedure cr_solve_forward_step_rr_g
    module procedure cr_solve_forward_step_rc_g
    module procedure cr_solve_forward_step_cc_g
  end interface cr_solve_forward_step
  !
  interface cr_solve_backward_step
    module procedure cr_solve_backward_step_rr
    module procedure cr_solve_backward_step_rc
    module procedure cr_solve_backward_step_cc
    module procedure cr_solve_backward_step_rr_g
    module procedure cr_solve_backward_step_rc_g
    module procedure cr_solve_backward_step_cc_g
  end interface cr_solve_backward_step
  !
  contains
  !
  subroutine m3c_decompose_r(m,mf,fail)
    real(rk), intent(in)          :: m (:,:) ! Tridiagonal matrix, assumed to be diagonally dominant
    real(rk), intent(out)         :: mf(:,:) ! Factorized tridiagonal matrix
    logical, intent(out),optional :: fail    ! Set to .true. if decomposition fails; 
                                             ! if fail is absent, abort on decomposition failure.
    include 'tridiagonal_cyclic_m3c_decompose_common.f90'
  end subroutine m3c_decompose_r
  !
  subroutine m3c_decompose_c(m,mf,fail)
    complex(rk), intent(in)       :: m (:,:) ! Tridiagonal matrix, assumed to be diagonally dominant
    complex(rk), intent(out)      :: mf(:,:) ! Factorized tridiagonal matrix
    logical, intent(out),optional :: fail    ! Set to .true. if decomposition fails; 
                                             ! if fail is absent, abort on decomposition failure.
    include 'tridiagonal_cyclic_m3c_decompose_common.f90'
  end subroutine m3c_decompose_c
  !
  subroutine m3c_solve_rr(mf,rhs,x,scr)
    real(rk), intent(in)  :: mf(:,:)  ! Factorization data
    real(rk), intent(in)  :: rhs(:)   ! Right-hand side
    real(rk), intent(out) :: x(:)     ! Solution of the linear system
    real(rk), intent(out) :: scr(:,:) ! Scratch space
    !
    include 'tridiagonal_cyclic_m3c_solve_common.f90'
  end subroutine m3c_solve_rr
  !
  subroutine m3c_solve_rc(mf,rhs,x,scr)
    real(rk), intent(in)     :: mf(:,:)  ! Factorization data
    complex(rk), intent(in)  :: rhs(:)   ! Right-hand side
    complex(rk), intent(out) :: x(:)     ! Solution of the linear system
    complex(rk), intent(out) :: scr(:,:) ! Scratch space
    !
    include 'tridiagonal_cyclic_m3c_solve_common.f90'
  end subroutine m3c_solve_rc
  !
  subroutine m3c_solve_cc(mf,rhs,x,scr)
    complex(rk), intent(in)  :: mf(:,:)  ! Factorization data
    complex(rk), intent(in)  :: rhs(:)   ! Right-hand side
    complex(rk), intent(out) :: x(:)     ! Solution of the linear system
    complex(rk), intent(out) :: scr(:,:) ! Scratch space
    !
    include 'tridiagonal_cyclic_m3c_solve_common.f90'
  end subroutine m3c_solve_cc
  !
  subroutine m3c_solve_rr_g(mf,rhs,x,scr)
    real(rk), intent(in)  :: mf(:,:)  ! Factorization data
    real(rk), intent(in)  :: rhs(:,:) ! Right-hand sides
    real(rk), intent(out) :: x(:,:)   ! Solutions of the linear system
    real(rk), intent(out) :: scr(:,:) ! Scratch space
    !
    include 'tridiagonal_cyclic_m3c_solve_g_common.f90'
  end subroutine m3c_solve_rr_g
  !
  subroutine m3c_solve_rc_g(mf,rhs,x,scr)
    real(rk), intent(in)     :: mf(:,:)  ! Factorization data
    complex(rk), intent(in)  :: rhs(:,:) ! Right-hand sides
    complex(rk), intent(out) :: x(:,:)   ! Solutions of the linear system
    complex(rk), intent(out) :: scr(:,:) ! Scratch space
    !
    include 'tridiagonal_cyclic_m3c_solve_g_common.f90'
  end subroutine m3c_solve_rc_g
  !
  subroutine m3c_solve_cc_g(mf,rhs,x,scr)
    complex(rk), intent(in)  :: mf(:,:)  ! Factorization data
    complex(rk), intent(in)  :: rhs(:,:) ! Right-hand sides
    complex(rk), intent(out) :: x(:,:)   ! Solutions of the linear system
    complex(rk), intent(out) :: scr(:,:) ! Scratch space
    !
    include 'tridiagonal_cyclic_m3c_solve_g_common.f90'
  end subroutine m3c_solve_cc_g
  !
  !  Internal subroutines below this point
  !
  subroutine cr_reduce_step_r(m,f,fail)
    real(rk), intent(in)  :: m (:,:) ! Tridiagonal matrix to be reduced
    real(rk), intent(out) :: f(:,:)  ! Reduction data.
    logical, intent(out)  :: fail
    !
    include 'tridiagonal_cyclic_cr_reduce_step_common.f90'
  end subroutine cr_reduce_step_r
  !
  subroutine cr_reduce_step_c(m,f,fail)
    complex(rk), intent(in)  :: m (:,:) ! Tridiagonal matrix to be reduced
    complex(rk), intent(out) :: f(:,:)  ! Reduction data.
    logical, intent(out)     :: fail
    !
    include 'tridiagonal_cyclic_cr_reduce_step_common.f90'
  end subroutine cr_reduce_step_c
  !
  subroutine cr_solve_forward_step_rr(f,r,rr)
    real(rk), intent(in)  :: f(:,:) ! Reduction data from cr_reduce_step
    real(rk), intent(in)  :: r (:)  ! Full RHS
    real(rk), intent(out) :: rr(:)  ! Reduced RHS on this step
    !
    include 'tridiagonal_cyclic_cr_solve_forward_common.f90'
  end subroutine cr_solve_forward_step_rr
  !
  subroutine cr_solve_forward_step_rc(f,r,rr)
    real(rk), intent(in)     :: f(:,:) ! Reduction data from cr_reduce_step
    complex(rk), intent(in)  :: r (:)  ! Full RHS
    complex(rk), intent(out) :: rr(:)  ! Reduced RHS on this step
    !
    include 'tridiagonal_cyclic_cr_solve_forward_common.f90'
  end subroutine cr_solve_forward_step_rc
  !
  subroutine cr_solve_forward_step_cc(f,r,rr)
    complex(rk), intent(in)  :: f(:,:) ! Reduction data from cr_reduce_step
    complex(rk), intent(in)  :: r (:)  ! Full RHS
    complex(rk), intent(out) :: rr(:)  ! Reduced RHS on this step
    !
    include 'tridiagonal_cyclic_cr_solve_forward_common.f90'
  end subroutine cr_solve_forward_step_cc
  !
  subroutine cr_solve_forward_step_rr_g(f,r,rr)
    real(rk), intent(in)  :: f (:,:) ! Reduction data from cr_reduce_step
    real(rk), intent(in)  :: r (:,:) ! Full RHS set
    real(rk), intent(out) :: rr(:,:) ! Reduced RHS set on this step
    !
    include 'tridiagonal_cyclic_cr_solve_forward_g_common.f90'
  end subroutine cr_solve_forward_step_rr_g
  !
  subroutine cr_solve_forward_step_rc_g(f,r,rr)
    real(rk), intent(in)     :: f (:,:) ! Reduction data from cr_reduce_step
    complex(rk), intent(in)  :: r (:,:) ! Full RHS set
    complex(rk), intent(out) :: rr(:,:) ! Reduced RHS set on this step
    !
    include 'tridiagonal_cyclic_cr_solve_forward_g_common.f90'
  end subroutine cr_solve_forward_step_rc_g
  !
  subroutine cr_solve_forward_step_cc_g(f,r,rr)
    complex(rk), intent(in)  :: f (:,:) ! Reduction data from cr_reduce_step
    complex(rk), intent(in)  :: r (:,:) ! Full RHS set
    complex(rk), intent(out) :: rr(:,:) ! Reduced RHS set on this step
    !
    include 'tridiagonal_cyclic_cr_solve_forward_g_common.f90'
  end subroutine cr_solve_forward_step_cc_g
  !
  subroutine cr_solve_backward_step_rr(f,r,xr,x)
    real(rk), intent(in)  :: f(:,:) ! Reduction data from cr_reduce_step
    real(rk), intent(in)  :: r (:)  ! RHS
    real(rk), intent(in)  :: xr(:)  ! Reduced solution vector
    real(rk), intent(out) :: x (:)  ! Full Solution vector
    !
    include 'tridiagonal_cyclic_cr_solve_backward_common.f90'
  end subroutine cr_solve_backward_step_rr
  !
  subroutine cr_solve_backward_step_rc(f,r,xr,x)
    real(rk), intent(in)     :: f(:,:) ! Reduction data from cr_reduce_step
    complex(rk), intent(in)  :: r (:)  ! RHS
    complex(rk), intent(in)  :: xr(:)  ! Reduced solution vector
    complex(rk), intent(out) :: x (:)  ! Full Solution vector
    !
    include 'tridiagonal_cyclic_cr_solve_backward_common.f90'
  end subroutine cr_solve_backward_step_rc
  !
  subroutine cr_solve_backward_step_cc(f,r,xr,x)
    complex(rk), intent(in)  :: f(:,:) ! Reduction data from cr_reduce_step
    complex(rk), intent(in)  :: r (:)  ! RHS
    complex(rk), intent(in)  :: xr(:)  ! Reduced solution vector
    complex(rk), intent(out) :: x (:)  ! Full Solution vector
    !
    include 'tridiagonal_cyclic_cr_solve_backward_common.f90'
  end subroutine cr_solve_backward_step_cc
  !
  subroutine cr_solve_backward_step_rr_g(f,r,xr,x)
    real(rk), intent(in)  :: f (:,:) ! Reduction data from cr_reduce_step
    real(rk), intent(in)  :: r (:,:) ! RHS set
    real(rk), intent(in)  :: xr(:,:) ! Reduced solution vectors
    real(rk), intent(out) :: x (:,:) ! Full Solution vectors
    !
    include 'tridiagonal_cyclic_cr_solve_backward_g_common.f90'
  end subroutine cr_solve_backward_step_rr_g
  !
  subroutine cr_solve_backward_step_rc_g(f,r,xr,x)
    real(rk), intent(in)     :: f (:,:) ! Reduction data from cr_reduce_step
    complex(rk), intent(in)  :: r (:,:) ! RHS set
    complex(rk), intent(in)  :: xr(:,:) ! Reduced solution vectors
    complex(rk), intent(out) :: x (:,:) ! Full Solution vectors
    !
    include 'tridiagonal_cyclic_cr_solve_backward_g_common.f90'
  end subroutine cr_solve_backward_step_rc_g
  !
  subroutine cr_solve_backward_step_cc_g(f,r,xr,x)
    complex(rk), intent(in)  :: f (:,:) ! Reduction data from cr_reduce_step
    complex(rk), intent(in)  :: r (:,:) ! RHS set
    complex(rk), intent(in)  :: xr(:,:) ! Reduced solution vectors
    complex(rk), intent(out) :: x (:,:) ! Full Solution vectors
    !
    include 'tridiagonal_cyclic_cr_solve_backward_g_common.f90'
  end subroutine cr_solve_backward_step_cc_g
  !
end module tridiagonal_cyclic
