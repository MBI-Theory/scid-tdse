!
!   SCID-TDSE: Simple 1-electron atomic TDSE solver
!   Copyright (C) 2015-2024 Serguei Patchkovskii, Serguei.Patchkovskii@mbi-berlin.de
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
!   Laser field routines
!
module spherical_tdse_field
  use spherical_tdse_data
  use node_tools
  implicit none
  private
  public rcsid_spherical_tdse_field
  public initialize_field
  !
  character(len=clen), save :: rcsid_spherical_tdse_field = "$Id: spherical_tdse_field.f90,v 1.1 2024/02/13 14:22:14 ps Exp $"
  !
  contains
  !
  subroutine initialize_field
    call fill_vpot_table
    call unwrap_vpot_table
    if (vp_as_is) then
      call fill_efield_table_single_sided
    else
      call fill_efield_table
    end if
    if (nts%this_node==1) then
      call preview_laser_field 
    end if
  end subroutine initialize_field
  !
  subroutine fill_vpot_table
    integer(ik) :: its, alloc, iu_temp
    real(xk)    :: time, vp, th, ph
    !
    call TimerStart('Prepare VP table')
    allocate (vpot_table(0:3,-1:2*timesteps+1),stat=alloc)
    if (alloc/=0) stop 'spherical_tdse_field%fill_vpot_table - allocation failed'
    if (vp_shape=='table') then
      open (newunit=iu_temp,form='formatted',recl=256,action='read',position='rewind',status='old',file=trim(vp_table))
      read (iu_temp,*) vpot_table
      close (iu_temp)
    else
      vp = vp_apot(0._xk) ! Initialize internal structures of vp_apot
      !
      !  Start by filling the table 
      !
      time_steps: do its=-1,2*timesteps+1
        time = dt*0.5_xk*its
        vp   = vp_apot(time,theta=th,phi=ph)
        vpot_table(:,its) = (/ time, vp, th, ph /)
      end do time_steps
    end if
    !
    !  Issue a warning if vpot_table at its=-1,0 and 2*timesteps is not zero
    !
    if (any(vpot_table(1:3,(/-1,0,2*timesteps,2*timesteps+1/))/=0._rk)) then
      write (out,"(/'WARNING: Vector potential at time steps <= 0 and/or the end of the simulation is not zero')")
      if (.not. vp_as_is) then
        write (out,"( 'WARNING: Resetting vector-potential to zero, and initial/final orientation to lab.'/)")
      else
        write (out,"( 'WARNING: USING VECTOR-POTENTIAL AS IS. THE RESULTS ARE LIKELY UNPHYSICAL.'/)")
      end if
    end if
    !
    !  Force vector-potential to be zero at its=-1,0 and 2*timesteps
    !  Force local coordinate system at these time steps to be the laboratory system
    !
    if (.not. vp_as_is) then
      vpot_table(1:3,(/-1,0,2*timesteps,2*timesteps+1/)) = 0._rk
    end if
    !
    call TimerStop('Prepare VP table')
  end subroutine fill_vpot_table
  !
  subroutine unwrap_vpot_table
    integer(ik) :: its
    real(xk)    :: vp, th, ph, vp_extrap, th_ref, ph_ref
    real(xk)    :: th_v1, ph_v1, th_v2, ph_v2, r2_v1, r2_v2
    real(xk)    :: vp_max
    !
    if (.not.field_unwrap .or. rotation_mode=='none') then
      write (out,"(/'Skipping vector-potential unwrapping'/)")
      return
    end if
    call TimerStart('Unwrap VP table')
    !
    !  Spherical representation is not necessarily continuous in the normal parameter domain
    !  vp [0:inf), th [0:pi], ph [0:2pi)
    !  Our propagator is much happier when all three parameters to be continuos though the 
    !  first order (althogh it should be able to handle rapid change in theta and phi), 
    !  so that we need to unwrap the v.p.
    !
    !  Step 1: make sure V.P. magnitude's derivative is continuous. Points -1 and 0 are
    !          assumed to be "good" already, and will be used to start the process.
    ! 
    vp_max = maxval(abs(vpot_table(1,:)))
    unwrap_magnitude: do its=1,2*timesteps+1
      vp_extrap = 2*vpot_table(1,its-1) - vpot_table(1,its-2)
      vp        =   vpot_table(1,its)
      ! Try to avoid messing with vector-potential when it is nearly zero
      if (abs(vp_extrap+vp)>abs(vp_extrap-vp) .or. (abs(vp)+abs(vp_extrap))<unwrap_threshold*vp_max) cycle unwrap_magnitude
      ! We get smoother vector-potential by flipping the magnitude; do it!
      vpot_table(1,its) = -vpot_table(1,its)
      vpot_table(2,its) = -pi_xk + vpot_table(2,its)
    end do unwrap_magnitude
    !
    !  Step 2: keeping VP magnitude constant, try to choose the smoothest possible version of (theta,phi)
    !
    !   We have two main possibilities:
    !   1.  th =  th0 + 2pi*n ; ph = ph0 + 2pi*k
    !   2   th = -th0 + 2pi*n ; ph = ph0 + 2pi*k + pi
    !
    unwrap_theta_phi: do its=0,2*timesteps+1
      th     = vpot_table(2,its) 
      ph     = vpot_table(3,its) 
      th_ref = vpot_table(2,its-1)
      ph_ref = vpot_table(3,its-1)
      th_v1  = nearest_2pi( th,      th_ref)
      ph_v1  = nearest_2pi( ph,      ph_ref)
      th_v2  = nearest_2pi(-th,      th_ref)
      ph_v2  = nearest_2pi( ph+pi_xk,ph_ref)
      r2_v1  = (th_v1-th_ref)**2 + (ph_v1-ph_ref)**2
      r2_v2  = (th_v2-th_ref)**2 + (ph_v2-ph_ref)**2
      if (r2_v1<=r2_v2) then
        vpot_table(2,its) = th_v1
        vpot_table(3,its) = ph_v1
      else
        vpot_table(2,its) = th_v2
        vpot_table(3,its) = ph_v2
      end if
    end do unwrap_theta_phi
    call TimerStop('Unwrap VP table')
    contains
    function nearest_2pi(x,ref) result(v)
      real(xk), intent(in) :: x    
      real(xk), intent(in) :: ref 
      real(xk)             :: v
      !
      v = x + 2*pi_xk*nint((ref-x)/(2._xk*pi_xk),kind=ik)
    end function nearest_2pi
  end subroutine unwrap_vpot_table
  !
  !  We calculate the electric field using numerical differentiation of the 
  !  vector-potential. Since we require high numerical accuracy, we use
  !  implicit derivative expressions.
  !
  !  WARNING: This routine implicitly assumes that the vector-potential starts
  !  WARNING: at zero and ends at zero. It will produce garbage at the ends of
  !  WARNING: the interval if this condition is violated
  !
  subroutine fill_efield_table
    integer(ik)            :: its, alloc
    real(xk), allocatable  :: vpot_xyz(:,:)              ! Vector-potential in Cartesian laboratory coordinates
    real(xk), allocatable  :: dt(:)                      ! dt(i) is the timestep leading to the time at point (i)
    real(xk), allocatable  :: d1(:,:), m1(:,:), m1f(:,:) ! tri-diagonal matrices defining the implict gradient
    logical, allocatable   :: m1p(:)                     ! pivot table 
    real(xk), allocatable  :: tmp1(:,:), tmp2(:,:)       ! Temporary arrays
    real(xk), allocatable  :: scr(:,:)                   ! Scratch for the linear solver
    !
    call TimerStart('Prepare efield table')
    !
    !  efield_table will remain until the end of the run. Remaining arrays must be
    !  deallocated before leaving fill_efied_table.
    !
    allocate (efield_table(3,0:2*timesteps), &
              vpot_xyz(0:2*timesteps,3), dt(0:2*timesteps+1), &
              d1(0:2*timesteps,3), m1(0:2*timesteps,3), m1f(0:2*timesteps,m3d_dc_size), m1p(0:2*timesteps), &
              tmp1(0:2*timesteps,3), tmp2(0:2*timesteps,3), scr(0:2*timesteps,3*m3d_sc_size), &
              stat=alloc)
    if (alloc/=0) stop 'spherical_tdse%fill_efield_table - allocation failed'
    !
    !  Calculate Cartesian vector-potential using spherical vector-potential in vp_apot()
    !  Also build the table of the timesteps and the rotation matrices.
    !
    fill_vpot_xyz: do its=0,2*timesteps
      vpot_xyz(its,:) = vpot_sph2xyz(vpot_table(1:3,its))
      dt(its)         = vpot_table(0,its) - vpot_table(0,its-1)
    end do fill_vpot_xyz
    dt(2*timesteps+1) = vpot_table(0,2*timesteps+1)-vpot_table(0,2*timesteps)
    !
    !  We need to prepare the Delta_1 and M_1 matrices for the first derivative.
    !  The code snipped below is lifted from initialize_radial_gradient() in spherical_data.f90
    !  We can't just reuse initialize_radial_gradient, since the boundary conditions are different
    !
    grad_tables: do its=0,2*timesteps
      ! Diagonal
      d1(its,1) = 1/dt(its) - 1/dt(its+1) + (-dt(its) + dt(its+1))/ (dt(its)**2 + dt(its)*dt(its+1) + dt(its+1)**2)
      m1(its,1) = (dt(its) + dt(its+1))**2/ (2*(dt(its)**2 + dt(its)*dt(its+1) + dt(its+1)**2))
      if (its>=2*timesteps) cycle grad_tables
      ! Sub-diagonal
      d1(its,2) = -((dt(its+2)**2*(2*dt(its+1) + dt(its+2)))/ (dt(its+1)*(dt(its+1) + dt(its+2))* &
                         (dt(its+1)**2 + dt(its+1)*dt(its+2) + dt(its+2)**2)))
      m1(its,2) = dt(its+2)**2/(2*(dt(its+1)**2 + dt(its+1)*dt(its+2) + dt(its+2)**2))
      ! Super-diagonal
      d1(its,3) = (dt(its)**2*(dt(its) + 2*dt(its+1)))/(dt(its+1)* (dt(its) + dt(its+1))*(dt(its)**2 + &
                         dt(its)*dt(its+1) + dt(its+1)**2))
      m1(its,3) = dt(its)**2/(2*(dt(its)**2 + dt(its)*dt(its+1) + dt(its+1)**2))
    end do grad_tables
    call m3d_decompose_x(m1,m1f,m1p)
    !
    !  We are done with all the matrices; the derivative is now given by:
    !
    !    M1^-1 Delta_1 A
    !
    call m3d_multiply_x(d1,vpot_xyz,tmp1)
    call m3d_solve_x(m1,m1f,m1p,tmp1,tmp2,scr)
    efield_table = -transpose(tmp2)
    !
    deallocate (vpot_xyz,dt,d1,m1,m1f,m1p,tmp1,tmp2,scr)
    call TimerStop('Prepare efield table')
  end subroutine fill_efield_table
  !
  !  Differentiate vector-potential without making assumptions about its
  !  values outside the time domain. Less accurate than fill_efield_table()
  !
  subroutine fill_efield_table_single_sided
    integer(ik)            :: its, alloc
    integer(ik)            :: il, ih
    real(xk), allocatable  :: vpot_xyz(:,:)              ! Vector-potential in Cartesian laboratory coordinates
    !
    call TimerStart('Prepare efield table (special)')
    !
    !  efield_table will remain until the end of the run. Remaining arrays must be
    !  deallocated before leaving fill_efied_table.
    !
    allocate (efield_table(3,0:2*timesteps), vpot_xyz(3,-1:2*timesteps+1), &
              stat=alloc)
    if (alloc/=0) stop 'spherical_tdse%fill_efield_table_single_sided - allocation failed'
    !
    !  Calculate Cartesian vector-potential using spherical vector-potential in vp_apot()
    !  Also build the table of the timesteps.
    !
    fill_vpot_xyz: do its=-1,2*timesteps+1
      vpot_xyz(:,its) = vpot_sph2xyz(vpot_table(1:3,its))
    end do fill_vpot_xyz
    !
    differentiate_vpot: do its=0,2*timesteps
      il = max(its-2,-1)
      ih = min(its+2,2*timesteps+1)
      efield_table(:,its) = - poly_gradient(vpot_table(0,il:ih),vpot_xyz(:,il:ih),its-il+1)
    end do differentiate_vpot
    !
    deallocate (vpot_xyz)
    call TimerStop('Prepare efield table (special)')
  end subroutine fill_efield_table_single_sided
  !
  !  Evaluate gradient at a grid point from the first derivative of the Lagrange interpolant
  !
  function poly_gradient(x,f,pos) result (g)
    real(xk), intent(in)    :: x(  :)             ! Time grid
    real(xk), intent(in)    :: f(:,:)             ! (Vector) function on a time frid
    integer(ik), intent(in) :: pos                ! Time point where we need the gradient
    real(xk)                :: g(size(f,dim=1))   ! Time derivative of the function
    !
    integer(ik) :: npts
    integer(ik) :: i, j
    real(xk)    :: term(size(f,dim=1))
    !
    npts = size(x)
    if (size(f,dim=1)<=0)    stop 'spherical_tdse%poly_gradient - bad vector length'
    if (size(f,dim=2)/=npts) stop 'spherical_tdse%poly_gradient - inconsistent sizes'
    if (pos<1 .or. pos>npts) stop 'spherical_tdse%poly_gradient - bad pos'
    !
    g(:) = 0
    accumulate_gradient: do i=1,npts
      if (i==pos) cycle accumulate_gradient
      term = f(:,i) - f(:,pos)
      accumulate_weights: do j=1,npts
        if (j==i) cycle accumulate_weights
        term = term / (x(i)-x(j))
        if (j==pos) cycle accumulate_weights
        term = term * (x(pos)-x(j))
      end do accumulate_weights
      g = g + term
    end do accumulate_gradient
  end function poly_gradient
  !
  function vpot_sph2xyz(atp)
    real(xk), intent(in) :: atp(3)
    real(xk)             :: vpot_sph2xyz(3)
    !
    vpot_sph2xyz(1) = atp(1)*sin(atp(2))*cos(atp(3))
    vpot_sph2xyz(2) = atp(1)*sin(atp(2))*sin(atp(3))
    vpot_sph2xyz(3) = atp(1)*cos(atp(2))
  end function vpot_sph2xyz
  !
  subroutine preview_laser_field
    integer(ik) :: its, iu_temp
    real(xk)    :: tatp(0:3) ! (time, a, theta, phi)
    real(xk)    :: efield(3) ! (Ex,Ey,Ez)
    !
    if (field_preview==' ') return
    call TimerStart('Dump VP table')
    !
    open (newunit=iu_temp,form='formatted',action='write',position='rewind',status='replace',file=trim(field_preview))
    !
    write (iu_temp,"(('#',a11,1x,a26,3(1x,a30,4x),1x,3(1x,a20,4x)))") &
           ' step ', ' Time, au[t] ', ' V.P. magnitude ', ' V.P. theta ', ' V.P. phi ', ' Ex ', ' Ey ', ' Ez ', &
           '------', '-------------', '----------------', '------------', '----------', '----', '----', '----'
    time_steps: do its=0,2*timesteps
      tatp   = vpot_table(:,its)
      efield = efield_table(:,its)
      write (iu_temp,"(1x,f11.1,1x,f26.16,3(1x,g30.19e3),1x,3(1x,g24.13e3))") 0.5_rk*its, tatp, efield
    end do time_steps
    !
    close (iu_temp)
    call TimerStop('Dump VP table')
  end subroutine preview_laser_field
end module spherical_tdse_field
