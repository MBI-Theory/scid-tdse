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
!  Simple tools for dealing with three-diagonal matrices and systems of equations.
!
!  All tridiagonal matrices used by this module are in the following format:
!   m(1,i) = Main diagonal; element (i,i)
!   m(2,i) = Sub-diagonal; element (i+1,i)
!   m(3,i) = Super-diagonal; element (i,i+1)
!
module tridiagonal_tools
  use accuracy
  use constants
  use timer
  use tridiagonal_pivoted
  use tridiagonal_cyclic
  implicit none
  private
  public m3d_solver, m3d_iterations
  public m3d_dc_size, m3d_sc_size
  public m3d_multiply, m3d_decompose, m3d_solve
  public m3d_multiply_x, m3d_decompose_x, m3d_solve_x
  public m3u_decompose, m3u_solve
  public m3u_decompose_x, m3u_solve_x
  public m3d_left_scale, m3d_right_scale, m3d_transpose
  public rcsid_tridiagonal_tools
  !
  character(len=clen), save :: rcsid_tridiagonal_tools = "$Id: tridiagonal_tools.f90,v 1.22 2021/04/26 15:44:44 ps Exp $"
  !
! character(len=10), save :: m3d_solver     = 'cyclic'     ! Choice of the linear solver. Can be one of:
  character(len=10), save :: m3d_solver     = 'unpivoted'  ! Choice of the linear solver. Can be one of:
                                                           ! 'unpivoted' = use a simple, unpivoted solver.
                                                           ! 'cyclic'    = use an unpivoted cyclic-reduction sover. Not
                                                           !               recommended on scalar or short-vector
                                                           !               architectures.
                                                           ! 'pivoted'   = use a pivoted solver
                                                           ! 'refined'   = use a pivoted solver with iterative refinement
                                                           ! 'unpivoted' is the fastest and is the default.
                                                           ! 'refined' is the slowest and most accurate.
  integer(ik), save       :: m3d_iterations = 2_ik         ! Number of iterations for m3d_solver=='refined'
  integer, parameter      :: m3d_dc_size    = 8            ! Largest array size needed for decomposition for any solver
  integer, parameter      :: m3d_sc_size    = 2            ! Largest array size needed for scratch for any solver
  !
  interface m3d_multiply
    module procedure m3d_multiply_rr
    module procedure m3d_multiply_rc
    module procedure m3d_multiply_cc
    module procedure m3d_multiply_rr_g
    module procedure m3d_multiply_rc_g
    module procedure m3d_multiply_cc_g
  end interface m3d_multiply
  interface m3d_multiply_x
    module procedure m3d_multiply_xrr_g
  end interface m3d_multiply_x
  !
  interface m3d_decompose
    module procedure m3d_decompose_r
    module procedure m3d_decompose_c
  end interface m3d_decompose
  interface m3d_decompose_x
    module procedure m3d_decompose_xr
  end interface m3d_decompose_x
  !
  interface m3d_solve
    module procedure m3d_solve_rr
    module procedure m3d_solve_rc
    module procedure m3d_solve_cc
    module procedure m3d_solve_rr_g
    module procedure m3d_solve_rc_g
    module procedure m3d_solve_cc_g
  end interface m3d_solve
  interface m3d_solve_x
    module procedure m3d_solve_xrr_g
  end interface m3d_solve_x
  !
  interface m3u_decompose
    module procedure m3u_decompose_r
    module procedure m3u_decompose_c
  end interface m3u_decompose
  interface m3u_decompose_x
    module procedure m3u_decompose_xr
  end interface m3u_decompose_x
  !
  interface m3u_solve
    module procedure m3u_solve_rr
    module procedure m3u_solve_rc
    module procedure m3u_solve_cc
    module procedure m3u_solve_rr_g
    module procedure m3u_solve_rc_g
    module procedure m3u_solve_cc_g
  end interface m3u_solve 
  interface m3u_solve_x
    module procedure m3u_solve_xrr_g
  end interface m3u_solve_x
  !
  interface m3d_left_scale
    module procedure m3d_left_scale_cr
  end interface m3d_left_scale
  !
  interface m3d_right_scale
    module procedure m3d_right_scale_rc
  end interface m3d_right_scale
  !
  interface m3d_transpose
    module procedure m3d_transpose_r
  end interface m3d_transpose
  !
  contains
  !
  subroutine m3d_multiply_rr(m,v,mv)
    real(rk), intent(in)  :: m(:,:) ! Tri-diagonal matrix
    real(rk), intent(in)  :: v(:)   ! Vector 
    real(rk), intent(out) :: mv(:)  ! Vector 
    !
    include "tridiagonal_tools_m3d_multiply_common.f90"
  end subroutine m3d_multiply_rr
  !
  subroutine m3d_multiply_rc(m,v,mv)
    real(rk), intent(in)     :: m(:,:) ! Tri-diagonal matrix
    complex(rk), intent(in)  :: v(:)   ! Vector 
    complex(rk), intent(out) :: mv(:)  ! Vector 
    !
    include "tridiagonal_tools_m3d_multiply_common.f90"
  end subroutine m3d_multiply_rc
  !
  subroutine m3d_multiply_cc(m,v,mv)
    complex(rk), intent(in)  :: m(:,:) ! Tri-diagonal matrix
    complex(rk), intent(in)  :: v(:)   ! Vector 
    complex(rk), intent(out) :: mv(:)  ! Vector 
    !
    include "tridiagonal_tools_m3d_multiply_common.f90"
  end subroutine m3d_multiply_cc
  !
  subroutine m3d_multiply_rr_g(m,v,mv)
    real(rk), intent(in)  :: m(:,:)  ! Tri-diagonal matrix
    real(rk), intent(in)  :: v(:,:)  ! Vectors
    real(rk), intent(out) :: mv(:,:) ! Vectors
    !
    include "tridiagonal_tools_m3d_multiply_g_common.f90"
  end subroutine m3d_multiply_rr_g
  !
  subroutine m3d_multiply_xrr_g(m,v,mv)
    real(xk), intent(in)  :: m(:,:)  ! Tri-diagonal matrix
    real(xk), intent(in)  :: v(:,:)  ! Vectors
    real(xk), intent(out) :: mv(:,:) ! Vectors
    !
    include "tridiagonal_tools_m3d_multiply_g_common.f90"
  end subroutine m3d_multiply_xrr_g
  !
  subroutine m3d_multiply_rc_g(m,v,mv)
    real(rk), intent(in)     :: m(:,:)  ! Tri-diagonal matrix
    complex(rk), intent(in)  :: v(:,:)  ! Vectors
    complex(rk), intent(out) :: mv(:,:) ! Vectors
    !
    include "tridiagonal_tools_m3d_multiply_g_common.f90"
  end subroutine m3d_multiply_rc_g
  !
  subroutine m3d_multiply_cc_g(m,v,mv)
    complex(rk), intent(in)  :: m(:,:)  ! Tri-diagonal matrix
    complex(rk), intent(in)  :: v(:,:)  ! Vectors
    complex(rk), intent(out) :: mv(:,:) ! Vectors
    !
    include "tridiagonal_tools_m3d_multiply_g_common.f90"
  end subroutine m3d_multiply_cc_g
  !
  !  Linear solver drivers
  !
  subroutine m3d_decompose_r(m,mf,lpiv,fail)
    real(rk), intent(in)           :: m (:,:) ! Tridiagonal matrix
    real(rk), intent(out)          :: mf(:,:) ! Factorized tridiagonal matrix
    logical, intent(out)           :: lpiv(:) ! Pivot list
    logical, intent(out), optional :: fail    ! Set to .true. if decomposition fails; 
                                              ! if fail is absent, abort on decomposition failure.
    !
    include "tridiagonal_tools_m3d_decompose_common.f90"
  end subroutine m3d_decompose_r
  !
  subroutine m3d_decompose_c(m,mf,lpiv,fail)
    complex(rk), intent(in)        :: m (:,:) ! Tridiagonal matrix
    complex(rk), intent(out)       :: mf(:,:) ! Factorized tridiagonal matrix
    logical, intent(out)           :: lpiv(:) ! Pivot list
    logical, intent(out), optional :: fail    ! Set to .true. if decomposition fails; 
                                              ! if fail is absent, abort on decomposition failure.
    !
    include "tridiagonal_tools_m3d_decompose_common.f90"
  end subroutine m3d_decompose_c
  !
  subroutine m3d_decompose_xr(m,mf,lpiv,fail)
    real(xk), intent(in)           :: m (:,:) ! Tridiagonal matrix
    real(xk), intent(out)          :: mf(:,:) ! Factorized tridiagonal matrix
    logical, intent(out)           :: lpiv(:) ! Pivot list
    logical, intent(out), optional :: fail    ! Set to .true. if decomposition fails; 
                                              ! if fail is absent, abort on decomposition failure.
    !
    include "tridiagonal_tools_m3d_decompose_x_common.f90"
  end subroutine m3d_decompose_xr
  !
  !  Solver driver. We may use unpivoted, pivoted, or iteratively-refined pivoted
  !  routines, depending on the m3d_solver value.
  !
  subroutine m3d_solve_rr(m,mf,lpiv,r,x,scr)
    real(rk), intent(in)           :: m (:,:)  ! Tridiagonal matrix
    real(rk), intent(in)           :: mf(:,:)  ! Factorized tridiagonal matrix
    logical, intent(in)            :: lpiv(:)  ! Pivot list
    real(rk), intent(in)           :: r(:)     ! Right-hand side
    real(rk), intent(out)          :: x(:)     ! Solution of the linear system
    real(rk), intent(out)          :: scr(:,:) ! Scratch space for iterative refinement
    !
    include "tridiagonal_tools_m3d_solve_common.f90"
  end subroutine m3d_solve_rr
  !
  subroutine m3d_solve_rc(m,mf,lpiv,r,x,scr)
    real(rk), intent(in)           :: m (:,:)  ! Tridiagonal matrix
    real(rk), intent(in)           :: mf(:,:)  ! Factorized tridiagonal matrix
    logical, intent(in)            :: lpiv(:)  ! Pivot list
    complex(rk), intent(in)        :: r(:)     ! Right-hand side
    complex(rk), intent(out)       :: x(:)     ! Solution of the linear system
    complex(rk), intent(out)       :: scr(:,:) ! Scratch space for iterative refinement
    !
    include "tridiagonal_tools_m3d_solve_common.f90"
  end subroutine m3d_solve_rc
  !
  subroutine m3d_solve_cc(m,mf,lpiv,r,x,scr)
    complex(rk), intent(in)        :: m (:,:)  ! Tridiagonal matrix
    complex(rk), intent(in)        :: mf(:,:)  ! Factorized tridiagonal matrix
    logical, intent(in)            :: lpiv(:)  ! Pivot list
    complex(rk), intent(in)        :: r(:)     ! Right-hand side
    complex(rk), intent(out)       :: x(:)     ! Solution of the linear system
    complex(rk), intent(out)       :: scr(:,:) ! Scratch space for iterative refinement
    !
    include "tridiagonal_tools_m3d_solve_common.f90"
  end subroutine m3d_solve_cc
  !
  subroutine m3d_solve_rr_g(m,mf,lpiv,r,x,scr)
    real(rk), intent(in)           :: m (:,:)  ! Tridiagonal matrix
    real(rk), intent(in)           :: mf(:,:)  ! Factorized tridiagonal matrix
    logical, intent(in)            :: lpiv(:)  ! Pivot list
    real(rk), intent(in)           :: r(:,:)   ! Right-hand sides
    real(rk), intent(out)          :: x(:,:)   ! Solutions of the linear system
    real(rk), intent(out)          :: scr(:,:) ! Scratch space for iterative refinement
    !
    include "tridiagonal_tools_m3d_solve_g_common.f90"
  end subroutine m3d_solve_rr_g
  !
  subroutine m3d_solve_xrr_g(m,mf,lpiv,r,x,scr)
    real(xk), intent(in)           :: m (:,:)  ! Tridiagonal matrix
    real(xk), intent(in)           :: mf(:,:)  ! Factorized tridiagonal matrix
    logical, intent(in)            :: lpiv(:)  ! Pivot list
    real(xk), intent(in)           :: r(:,:)   ! Right-hand sides
    real(xk), intent(out)          :: x(:,:)   ! Solutions of the linear system
    real(xk), intent(out)          :: scr(:,:) ! Scratch space for iterative refinement
    !
    include "tridiagonal_tools_m3d_solve_x_g_common.f90"
  end subroutine m3d_solve_xrr_g
  !
  subroutine m3d_solve_rc_g(m,mf,lpiv,r,x,scr)
    real(rk), intent(in)           :: m (:,:)  ! Tridiagonal matrix
    real(rk), intent(in)           :: mf(:,:)  ! Factorized tridiagonal matrix
    logical, intent(in)            :: lpiv(:)  ! Pivot list
    complex(rk), intent(in)        :: r(:,:)   ! Right-hand sides
    complex(rk), intent(out)       :: x(:,:)   ! Solutions of the linear system
    complex(rk), intent(out)       :: scr(:,:) ! Scratch space for iterative refinement
    !
    include "tridiagonal_tools_m3d_solve_g_common.f90"
  end subroutine m3d_solve_rc_g
  !
  subroutine m3d_solve_cc_g(m,mf,lpiv,r,x,scr)
    complex(rk), intent(in)        :: m (:,:)  ! Tridiagonal matrix
    complex(rk), intent(in)        :: mf(:,:)  ! Factorized tridiagonal matrix
    logical, intent(in)            :: lpiv(:)  ! Pivot list
    complex(rk), intent(in)        :: r(:,:)   ! Right-hand sides
    complex(rk), intent(out)       :: x(:,:)   ! Solutions of the linear system
    complex(rk), intent(out)       :: scr(:,:) ! Scratch space for iterative refinement
    !
    include "tridiagonal_tools_m3d_solve_g_common.f90"
  end subroutine m3d_solve_cc_g
  !
  !  Unpivoted LU decomposition
  !
  subroutine m3u_decompose_r(m,mf,fail)
    real(rk), intent(in)           :: m (:,:) ! Tridiagonal matrix, assumed to be diagonally dominant
                                              ! As the result, we do not need to bother with pivoting
    real(rk), intent(out)          :: mf(:,:) ! Factorized tridiagonal matrix
    logical, intent(out), optional :: fail    ! Set to .true. if decomposition fails; 
                                              ! if fail is absent, abort on decomposition failure.
    !
    real(rk)    :: denom
    !
    include "tridiagonal_tools_m3u_decompose_common.f90"
  end subroutine m3u_decompose_r
  !
  subroutine m3u_decompose_xr(m,mf,fail)
    real(xk), intent(in)           :: m (:,:) ! Tridiagonal matrix, assumed to be diagonally dominant
                                              ! As the result, we do not need to bother with pivoting
    real(xk), intent(out)          :: mf(:,:) ! Factorized tridiagonal matrix
    logical, intent(out), optional :: fail    ! Set to .true. if decomposition fails; 
                                              ! if fail is absent, abort on decomposition failure.
    !
    real(xk)    :: denom
    !
    include "tridiagonal_tools_m3u_decompose_common.f90"
  end subroutine m3u_decompose_xr
  !
  subroutine m3u_decompose_c(m,mf,fail)
    complex(rk), intent(in)        :: m (:,:) ! Tridiagonal matrix, assumed to be diagonally dominant
                                              ! As the result, we do not need to bother with pivoting
    complex(rk), intent(out)       :: mf(:,:) ! Factorized tridiagonal matrix
    logical, intent(out), optional :: fail    ! Set to .true. if decomposition fails; 
                                              ! if fail is absent, abort on decomposition failure.
    !
    complex(rk) :: denom
    !
    include "tridiagonal_tools_m3u_decompose_common.f90"
  end subroutine m3u_decompose_c
  !
  subroutine m3u_solve_rr(mf,r,x)
    real(rk), intent(in)  :: mf(:,:) ! Factorized tridiagonal matrix
    real(rk), intent(in)  :: r (:)   ! Right-hand size
    real(rk), intent(out) :: x (:)   ! Solution vector
    !
    include "tridiagonal_tools_m3u_solve_common.f90"
  end subroutine m3u_solve_rr
  !
  subroutine m3u_solve_rc(mf,r,x)
    real(rk), intent(in)     :: mf(:,:) ! Factorized tridiagonal matrix
    complex(rk), intent(in)  :: r (:)   ! Right-hand size
    complex(rk), intent(out) :: x (:)   ! Solution vector
    !
    include "tridiagonal_tools_m3u_solve_common.f90"
  end subroutine m3u_solve_rc
  !
  subroutine m3u_solve_cc(mf,r,x)
    complex(rk), intent(in)  :: mf(:,:) ! Factorized tridiagonal matrix
    complex(rk), intent(in)  :: r (:)   ! Right-hand size
    complex(rk), intent(out) :: x (:)   ! Solution vector
    !
    include "tridiagonal_tools_m3u_solve_common.f90"
  end subroutine m3u_solve_cc
  !
  subroutine m3u_solve_rr_g(mf,r,x)
    real(rk), intent(in)  :: mf(:,:) ! Factorized tridiagonal matrix
    real(rk), intent(in)  :: r (:,:) ! Right-hand sizes
    real(rk), intent(out) :: x (:,:) ! Solution vectors
    !
    include "tridiagonal_tools_m3u_solve_g_common.f90"
  end subroutine m3u_solve_rr_g
  !
  subroutine m3u_solve_xrr_g(mf,r,x)
    real(xk), intent(in)  :: mf(:,:) ! Factorized tridiagonal matrix
    real(xk), intent(in)  :: r (:,:) ! Right-hand sizes
    real(xk), intent(out) :: x (:,:) ! Solution vectors
    !
    include "tridiagonal_tools_m3u_solve_g_common.f90"
  end subroutine m3u_solve_xrr_g
  !
  subroutine m3u_solve_rc_g(mf,r,x)
    real(rk), intent(in)     :: mf(:,:) ! Factorized tridiagonal matrix
    complex(rk), intent(in)  :: r (:,:) ! Right-hand sizes
    complex(rk), intent(out) :: x (:,:) ! Solution vectors
    !
    include "tridiagonal_tools_m3u_solve_g_common.f90"
  end subroutine m3u_solve_rc_g
  !
  subroutine m3u_solve_cc_g(mf,r,x)
    complex(rk), intent(in)  :: mf(:,:) ! Factorized tridiagonal matrix
    complex(rk), intent(in)  :: r (:,:) ! Right-hand sizes
    complex(rk), intent(out) :: x (:,:) ! Solution vectors
    !
    include "tridiagonal_tools_m3u_solve_g_common.f90"
  end subroutine m3u_solve_cc_g
  !
  subroutine m3d_left_scale_cr(s,m,sm)
    complex(rk), intent(in)  :: s   (:) ! Diagonal matrix
    real(rk), intent(in)     :: m (:,:) ! Tridiagonal matrix
    complex(rk), intent(out) :: sm(:,:) ! Tridiagonal s . m
    !
    include "tridiagonal_tools_m3d_left_scale_common.f90"
  end subroutine m3d_left_scale_cr
  !
  subroutine m3d_right_scale_rc(m,s,ms)
    real(rk), intent(in)     :: m (:,:) ! Tridiagonal matrix
    complex(rk), intent(in)  :: s   (:) ! Diagonal matrix
    complex(rk), intent(out) :: ms(:,:) ! Tridiagonal m . s
    !
    include "tridiagonal_tools_m3d_right_scale_common.f90"
  end subroutine m3d_right_scale_rc
  !
  subroutine m3d_transpose_r(m,mt)
    real(rk), intent(in)  :: m (:,:) ! Tridiagonal matrix
    real(rk), intent(out) :: mt(:,:) ! Transpose of m
    !
    if (any(ubound(m)/=ubound(mt))) then
      stop 'tridiagonal_tools%m3d_transpose_r - bad input sizes'
    end if
    !
    mt(:,1) = m(:,1)
    mt(:,2) = m(:,3)
    mt(:,3) = m(:,2)
  end subroutine m3d_transpose_r
end module tridiagonal_tools
