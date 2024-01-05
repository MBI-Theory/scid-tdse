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
!  Cubic-spline interpolation. The implementation follows Numerical Recipes 3.3.
!  We only implement the "natural" cubic splines, with the second derivatives
!  constrained to be zero at the outer edges of the interval.
!
!  By default, we'll use the pivoted tridiagonal solver to construct the interpolant,
!  which is not expected to be on the critical execution path. The global choice of
!  the solver can be applied by setting cs_use_default_solver=.true.
!
!  There does not appear to be any reason to prefer any of the solvers her - the 
!  differences are on the order of a few ulp.
!
module cubic_spline
  use accuracy
  use constants
  use timer
  use tridiagonal_pivoted
  use tridiagonal_tools
  implicit none
  private
  public cs_data, cs_use_default_solver
  public cs_initialize, cs_destroy, cs_evaluate
  public rcsid_cubic_spline
  !
  type cs_data
    integer               :: npts      ! Number of data points defining the fit
    real(xk), allocatable :: fit(:,:)  ! Data points defining the fit: (1:3,1:npts)
                                       ! The first index is:
                                       !   1 = Argument value [x]
                                       !   2 = Function value at this argument [f(x)]
                                       !   3 = 2nd derivative of the function
                                       ! The second index gives points, which are 
                                       ! expected to be in the ascending order
  end type cs_data
  !
  character(len=clen), save :: rcsid_cubic_spline    = "$Id: cubic_spline.f90,v 1.2 2021/04/26 15:44:44 ps Exp $"
  logical, save             :: cs_use_default_solver = .false. ! Use the default solver, not the pivoted.
  !
  contains
  !
  subroutine cs_destroy(csd)
    type(cs_data), intent(inout) :: csd
    !
    if (allocated(csd%fit)) deallocate(csd%fit)
    csd%npts = 0
  end subroutine cs_destroy
  !
  subroutine cs_initialize(csd,x,y)
    type(cs_data), intent(inout) :: csd
    real(xk), intent(in)         :: x(:), y(:)  ! x, y(x) pairs. X must be sorted in the ascending order.
                                                ! Specifying less than 2 (ie 0 or 1) points will make
                                                ! cs_evaluate() return zero for all arguments.
    !
    real(xk), allocatable :: d(:), m(:,:), mf(:,:), scr(:,:), rhs(:)
    logical, allocatable  :: lpiv(:)
    integer(ik)           :: sz, alloc
    logical               :: fail
    !
    sz = size(x)
    if (size(y)/=sz) stop 'cubic_spline%cs_initialize - x and y are not the same size'
    !
    !  Copy the input data to the fit table
    !
    csd%npts = sz
    if (sz<=1) return ! Special case: the fit is zero everywhere
    allocate (csd%fit(3,sz),stat=alloc)
    if (alloc/=0) stop 'cubic_spline%cs_initialize - allocation failed (1)'
    csd%fit(1,:) = x(:)
    csd%fit(2,:) = y(:)
    csd%fit(3,:) = 0
    !
    if (sz<=2) return ! Early return if possible
    !
    !  If we have interior points, we need to determine the value of the second derivative
    !  at these points, which requires solving a tridiagonal linear problem
    !
    allocate (d(sz-1),m(sz-2,3),mf(sz-2,m3d_dc_size),scr(sz-2,m3d_sc_size),lpiv(sz-2),rhs(sz-2),stat=alloc)
    if (alloc/=0) stop 'cubic_spline%cs_initialize - allocation failed (2)' 
    !
    !  Prepare linear-system we use to solve for the second derivatives
    !
    d(1:sz-1)   = x(2:sz) - x(1:sz-1)                    ! Spacing between the grid points
    m(1:sz-2,1) = (1._xk/3._xk) * (d(2:sz-1)+d(1:sz-2))  ! Diagonal
    m(1:sz-3,2) = (1._xk/6._xk) * d(2:sz-2)              ! Sub-diagonal
    m(  sz-2,2) = 0._xk                    
    m( :    ,3) = m(:,2)                                 ! Super-diagonal
    rhs(1:sz-2) = (y(1:sz-2)-y(2:sz-1))/d(1:sz-2) &
                + (y(3:sz  )-y(2:sz-1))/d(2:sz-1)        ! Right-hand side
    !
    !  Factor and solve the linear system
    !
    if (cs_use_default_solver) then
      call m3d_decompose(m,mf,lpiv,fail)
      if (fail) stop 'ubic_spline%cs_initialize - m3d_decompose() failed'
      call m3d_solve(m,mf,lpiv,rhs,csd%fit(3,2:sz-1),scr)
    else
      call m3dp_decompose(m,mf,lpiv,fail)
      if (fail) stop 'ubic_spline%cs_initialize - m3dp_decompose() failed'
      call m3dp_solve(mf,lpiv,rhs,csd%fit(3,2:sz-1))
    end if
    !
    deallocate (d,m,mf,scr,lpiv,rhs)
  end subroutine cs_initialize
  !
  function cs_evaluate(csd,x) result(y)
    type(cs_data), intent(in) :: csd
    real(xk), intent(in)      :: x
    real(xk)                  :: y
    !
    integer(ik) :: il, ir, ic
    real(xk)    :: d, t, u1, u2, u3, u4
    !
    !  Are we outsize of the interpolation range? If yes, the result is zero.
    !
    y = 0._xk
    if (csd%npts<=1) return
    if (x<csd%fit(1,1) .or. x>csd%fit(1,csd%npts)) return
    !
    !  Use bisection to find the interpolation patch
    !  We know that fit(1,il)<=x<=fit(1,ir) at the beginning.
    !  We will maintain the condition on each bisection.
    !
    il = 1 ; ir = csd%npts
    narrow_the_range: do while (ir>il+1)
      ic = (il+ir)/2
      if (csd%fit(1,ic)>x) then
        ir = ic
      else ! csd%fit(1,ic)<=x
        il = ic
      end if
    end do narrow_the_range
    if (ir-il/=1) stop 'cubic_spline%cs_evaluate - logic error'
    !
    !  Evaluate the reduced argument on this interolation patch.
    !  The reduced argument is in the [-0.5:0.5] range.
    !
    d  = csd%fit(1,il+1) - csd%fit(1,il)  ! Grid spacing
    t  = (x - csd%fit(1,il))/d - 0.5_xk   ! Reduced argument, [-0.5:0.5] range
    !
    !  Basis polynomials
    !
    u1 = 0.5_xk - t                                         ! Prop. to the function on the left edge
    u2 = 0.5_xk + t                                         ! Prop. to the function on the right edge
    u3 = (-1._xk/48._xk)*(1._xk-4._xk*t**2)*(3._xk-2._xk*t) ! Prop. to the second derivative on the left
    u4 = (-1._xk/48._xk)*(1._xk-4._xk*t**2)*(3._xk+2._xk*t) ! Prop. to the second derivative on the right
    !
    !  Now assemble the result
    !
    y = csd%fit(2,il)*u1 + csd%fit(2,il+1)*u2 + csd%fit(3,il)*u3*d**2 + csd%fit(3,il+1)*u4*d**2
  end function cs_evaluate
  !
end module cubic_spline
!
! program test
!   use accuracy
!   use cubic_spline
!   use tridiagonal_tools
!   !
!   integer(ik)           :: n, ios, i
!   real(xk)              :: vx, vy
!   real(xk), allocatable :: x(:), y(:)
!   type (cs_data)        :: csd
!   !
!   cs_use_default_solver = .true.
!   read(input,*) n, m3d_solver
!   allocate (x(n),y(n))
!   read(input,*) x, y
!   do i=1,100000
!     call cs_destroy(csd)
!     call cs_initialize(csd,x,y)
!   end do
!   interpol: do 
!     read(input,*,iostat=ios) vx
!     if (ios/=0) exit interpol
!     do i=1,10000
!       vy = cs_evaluate(csd,vx)
!     end do
!     write(out,"(1x,g48.32,1x,g48.32)") vx, vy
!   end do interpol
!   call cs_destroy(csd)
! end program test
