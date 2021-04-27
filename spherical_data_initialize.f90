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
!  Initialization for data structures in spherical_data.f90. These routines had to be
!  split off to avoid circular dependency on node_tools.
!
!  2017 Oct 10: Modify the routines to allow for on-the-fly changes to grid parameters
!               and repeated calls to the initializer. This is intended to allow dynamical
!               resizing of the radial grid (the adaptive L grid is already handled
!               efficiently with no need to explicitly resize it).
!
module spherical_data_initialize
  use accuracy
  use constants
  use timer
  use tridiagonal_tools
  use potential_tools
  use cap_tools
  use spherical_data
  use node_tools
  implicit none
  private
  public sd_initialize
  public sd_expand_implicit_operator, sd_laser_clm
  public rcsid_spherical_data_initialize
  !
  character(len=clen), save :: rcsid_spherical_data_initialize = &
     "$Id: spherical_data_initialize.f90,v 1.5 2021/04/26 15:44:44 ps Exp ps $"
  !
  interface sd_expand_implicit_operator
    module procedure sd_expand_implicit_operator_r
    module procedure sd_expand_implicit_operator_c
  end interface sd_expand_implicit_operator
  !
  integer, parameter          :: iu_temp        = 23         ! An arbitrary unit number, which can be used here
  !
  contains
  !
  subroutine sd_initialize(repeat)
    logical, intent(in) :: repeat ! Indicates grid-reinitialization. Repeated calls are quieter
    real(rk)            :: fac_lr
    !
    call TimerStart('Grid (re-)initialization')
    !
    if (sd_nradial_max<0) sd_nradial_max = sd_nradial
    if (sd_nradial>sd_nradial_max) then
      stop 'spherical_data_initialize%sd_initialize - reinitialization with sd_nradial>sd_nradial_max'
    end if
    if (sd_radial_edge<0) then
      sd_radial_edge = 8*rk_bytes()
      write (out,"(/'sd_radial_edge reset to ',i6)") sd_radial_edge
    end if
    !
    call allocate_arrays
    if (.not.repeat) then
      !
      !  On the first call, we will initialize the full grid, and store it for later use
      !
      call initialize_radial_grid
      sd_rtab_full  = sd_rtab
      sd_drtab_full = sd_drtab
      !
      !  Set wavefunction tolerances. fac_lr is the average grid spacing.
      !
      fac_lr = sum(sd_drtab_full)/size(sd_drtab_full)
      if (sd_tolerance_l(1)<0) sd_tolerance_l(1) = spacing(0.1_rk)
      if (sd_tolerance_l(2)<0) sd_tolerance_l(2) = sd_tolerance_l(1) * fac_lr
      if (sd_tolerance_r(1)<0) sd_tolerance_r(1) = spacing(0.1_rk)
      if (sd_tolerance_r(2)<0) sd_tolerance_r(2) = sd_tolerance_r(1) * fac_lr
      sd_adaptive = sd_adaptive_l .or. sd_adaptive_r
      !
      if (sd_adaptive_l) write (out,"('Right/left adaptive tolerance (L) = ',2e14.6)") sd_tolerance_l
      if (sd_adaptive_r) write (out,"('Right/left adaptive tolerance (R) = ',2e14.6)") sd_tolerance_r
      if (sd_adaptive)   write (out,"()")
    end if
    !
    !  Copy the relevant part of the grid
    !
    sd_rtab  = sd_rtab_full (:sd_nradial+1)
    sd_drtab = sd_drtab_full(:sd_nradial+1)
    !
    !  Calculate the inverse - we need it in a few places.
    !
    sd_rminus(:) = 1._rk / sd_rtab(:sd_nradial)
    !
    call initialize_radial_gradient(repeat)
    call initialize_radial_laplacian(repeat)
    call initialize_channel_potentials(repeat)
    !
    call TimerStop('Grid (re-)initialization')
  end subroutine sd_initialize 
  !
  !  Construct explict operator matrix for an operator given in implicit form:
  !
  !  Conceptually: block = mop . op
  !
  subroutine sd_expand_implicit_operator_r(op,mop,mopf,mopp,block)
    real(rk), intent(in)  :: op   (:,:) ! Tridiagonal operator matrix
    real(rk), intent(in)  :: mop  (:,:) ! Tridiagonal modifier operator
    real(rk), intent(in)  :: mopf (:,:) ! Factored modifier operator
    logical, intent(in)   :: mopp (:  ) ! Pivot list for the factored modifier
    real(rk), intent(out) :: block(:,:) ! Explicit operator matrix
    !
    real(rk)    :: rhs(size(block,dim=1))
    real(rk)    :: tmp(size(block,dim=1))
    real(rk)    :: scr(size(block,dim=1),m3d_sc_size)
    !
    include "spherical_data_expand_io_common.f90"
  end subroutine sd_expand_implicit_operator_r
  !
  subroutine sd_expand_implicit_operator_c(op,mop,mopf,mopp,block)
    complex(rk), intent(in)  :: op   (:,:) ! Tridiagonal operator matrix
    complex(rk), intent(in)  :: mop  (:,:) ! Tridiagonal modifier operator
    complex(rk), intent(in)  :: mopf (:,:) ! Factored modifier operator
    logical, intent(in)      :: mopp (:  ) ! Pivot list for the factored modifier
    complex(rk), intent(out) :: block(:,:) ! Explicit operator matrix
    !
    complex(rk) :: rhs(size(block,dim=1))
    complex(rk) :: tmp(size(block,dim=1))
    complex(rk) :: scr(size(block,dim=1),m3d_sc_size)
    !
    include "spherical_data_expand_io_common.f90"
  end subroutine sd_expand_implicit_operator_c
  !
  !  Angular coupling coefficients; see propagator_tools.f90 and wavefunction_tools.f90
  !
  function sd_laser_clm(l,m) result(clm)
    integer(ik), intent(in) :: l, m ! Quantum numbers
    real(xk)                :: clm  ! Coupling coefficient
    !
    clm = sqrt(real(l**2 - m**2,kind=xk)/real(4*l**2 - 1,kind=xk))
  end function sd_laser_clm
  !
  !  Internal routines below
  !
  subroutine allocate_arrays
    integer(ik) :: alloc
    !
    if (allocated(sd_rtab)   ) deallocate(sd_rtab)
    if (allocated(sd_drtab)  ) deallocate(sd_drtab)
    if (allocated(sd_pottab) ) deallocate(sd_pottab)
    if (allocated(sd_captab) ) deallocate(sd_captab)
    if (allocated(sd_dvdrtab)) deallocate(sd_dvdrtab)
    if (allocated(sd_rminus) ) deallocate(sd_rminus)
    if (allocated(sd_d2n_l0) ) deallocate(sd_d2n_l0)
    if (allocated(sd_m2n_l0) ) deallocate(sd_m2n_l0)
    if (allocated(sd_m2nf_l0)) deallocate(sd_m2nf_l0)
    if (allocated(sd_m2np_l0)) deallocate(sd_m2np_l0)
    if (allocated(sd_d2n_lx) ) deallocate(sd_d2n_lx)
    if (allocated(sd_m2n_lx) ) deallocate(sd_m2n_lx)
    if (allocated(sd_m2nf_lx)) deallocate(sd_m2nf_lx)
    if (allocated(sd_m2np_lx)) deallocate(sd_m2np_lx)
    if (allocated(sd_d1n_l0) ) deallocate(sd_d1n_l0)
    if (allocated(sd_m1n_l0) ) deallocate(sd_m1n_l0)
    if (allocated(sd_m1nf_l0)) deallocate(sd_m1nf_l0)
    if (allocated(sd_m1np_l0)) deallocate(sd_m1np_l0)
    if (allocated(sd_d1n_lx) ) deallocate(sd_d1n_lx)
    if (allocated(sd_m1n_lx) ) deallocate(sd_m1n_lx)
    if (allocated(sd_m1nf_lx)) deallocate(sd_m1nf_lx)
    if (allocated(sd_m1np_lx)) deallocate(sd_m1np_lx)
    if (allocated(sd_d2t_l0) ) deallocate(sd_d2t_l0)
    if (allocated(sd_m2t_l0) ) deallocate(sd_m2t_l0)
    if (allocated(sd_m2tf_l0)) deallocate(sd_m2tf_l0)
    if (allocated(sd_m2tp_l0)) deallocate(sd_m2tp_l0)
    if (allocated(sd_d2t_lx) ) deallocate(sd_d2t_lx)
    if (allocated(sd_m2t_lx) ) deallocate(sd_m2t_lx)
    if (allocated(sd_m2tf_lx)) deallocate(sd_m2tf_lx)
    if (allocated(sd_m2tp_lx)) deallocate(sd_m2tp_lx)
    if (allocated(sd_d1t_l0) ) deallocate(sd_d1t_l0)
    if (allocated(sd_m1t_l0) ) deallocate(sd_m1t_l0)
    if (allocated(sd_m1tf_l0)) deallocate(sd_m1tf_l0)
    if (allocated(sd_m1tp_l0)) deallocate(sd_m1tp_l0)
    if (allocated(sd_d1t_lx) ) deallocate(sd_d1t_lx)
    if (allocated(sd_m1t_lx) ) deallocate(sd_m1t_lx)
    if (allocated(sd_m1tf_lx)) deallocate(sd_m1tf_lx)
    if (allocated(sd_m1tp_lx)) deallocate(sd_m1tp_lx)
    !
    allocate (sd_rtab(sd_nradial+1),sd_drtab(sd_nradial+1), &
              sd_pottab(sd_nradial,1,0:sd_lmax), sd_captab(sd_nradial), &
              sd_dvdrtab(sd_nradial), sd_rminus(sd_nradial), &
              sd_d2n_l0(sd_nradial,3),sd_m2n_l0(sd_nradial,3), sd_m2nf_l0(sd_nradial,m3d_dc_size), sd_m2np_l0(sd_nradial), &
              sd_d2n_lx(sd_nradial,3),sd_m2n_lx(sd_nradial,3), sd_m2nf_lx(sd_nradial,m3d_dc_size), sd_m2np_lx(sd_nradial), &
              sd_d1n_l0(sd_nradial,3),sd_m1n_l0(sd_nradial,3), sd_m1nf_l0(sd_nradial,m3d_dc_size), sd_m1np_l0(sd_nradial), &
              sd_d1n_lx(sd_nradial,3),sd_m1n_lx(sd_nradial,3), sd_m1nf_lx(sd_nradial,m3d_dc_size), sd_m1np_lx(sd_nradial), &
              sd_d2t_l0(sd_nradial,3),sd_m2t_l0(sd_nradial,3), sd_m2tf_l0(sd_nradial,m3d_dc_size), sd_m2tp_l0(sd_nradial), &
              sd_d2t_lx(sd_nradial,3),sd_m2t_lx(sd_nradial,3), sd_m2tf_lx(sd_nradial,m3d_dc_size), sd_m2tp_lx(sd_nradial), &
              sd_d1t_l0(sd_nradial,3),sd_m1t_l0(sd_nradial,3), sd_m1tf_l0(sd_nradial,m3d_dc_size), sd_m1tp_l0(sd_nradial), &
              sd_d1t_lx(sd_nradial,3),sd_m1t_lx(sd_nradial,3), sd_m1tf_lx(sd_nradial,m3d_dc_size), sd_m1tp_lx(sd_nradial), &
              stat=alloc)
    if (alloc/=0) then
      write (out,"('spherical_data_initialize%allocate_array: memory allocation (1) failed with code ',i0)") alloc
      stop 'spherical_data_initialize%allocate_array - memory allocation failed (1)'
    end if
    if (.not.allocated(sd_rtab_full)) then
      allocate (sd_rtab_full(sd_nradial_max+1),sd_drtab_full(sd_nradial_max+1),stat=alloc)
      if (alloc/=0) then
        write (out,"('spherical_data_initialize%allocate_array: memory allocation (2) failed with code ',i0)") alloc
        stop 'spherical_data_initialize%allocate_array - memory allocation failed (2)'
      end if
    end if
  end subroutine allocate_arrays
  !
  subroutine initialize_radial_grid
    integer(ik)         :: ipt, npt_log
    integer             :: ios
    !
    write (out,"()")
    if (sd_rgrid_rmin<0) sd_rgrid_rmin = 0 ! Guard agains stupid user input, which can lead to mysterious errors
    if (sd_rgrid_rmin/=0) then
      write (out,"('Radial grid starts at ',g0.5,' instead of the origin')") sd_rgrid_rmin
      write (out,"('Infinite-barrier boundary conditions in effect.')")
    end if
    select case (sd_rgrid)
      case default
        write (out,"('spherical_data_initialize%initialize_radial_grid: Radial grid ',a,' is not recognized')") trim(sd_rgrid)
        stop 'spherical_data_initialize%initialize_radial_grid - bad grid'
      case ('uniform')
        write (out,"('Creating uniform grid with step ',g24.13,' Bohr')") sd_rgrid_dr
        fill_uniform_grid: do ipt=1,sd_nradial+1
          sd_rtab(ipt) = sd_rgrid_rmin + ipt*sd_rgrid_dr
        end do fill_uniform_grid
      case ('log')
        if (sd_rgrid_r0<=sd_rgrid_rmin) then
          write (out,"('spherical_data_initialize%initialize_radial_grid: r0 (',g0.5,') is less than rmin (',g0.5,')')") &
                 sd_rgrid_r0, sd_rgrid_rmin
          stop 'spherical_data_initialize%initialize_radial_grid - logarithmic grid starts inside infinite barrier!'
        end if
        if (sd_rgrid_npad<0) then
          sd_rgrid_npad = max(0,nint((1._rk-sd_rgrid_rmin/sd_rgrid_r0)/(sd_rgrid_scale-1._rk)))
          !
          !  Choose master node's result for sd_rgrid_npad. We do not want small numerical differences
          !  between the nodes (which -may- have different architectures) to produce different grids
          !  on different nodes!
          !
          call nt_broadcast(sd_rgrid_npad)
        end if
        write (out,"('Creating logarithmic grid starting at ',g24.13,' Bohr; scale factor ',g24.13)") sd_rgrid_r0, sd_rgrid_scale
        write (out,"('Inserting ',i0,'-point uniform padding grid at the origin')") sd_rgrid_npad
        pad_log_grid: do ipt=1,sd_rgrid_npad
          sd_rtab(ipt) = sd_rgrid_rmin + (ipt*(sd_rgrid_r0-sd_rgrid_rmin))/(sd_rgrid_npad+1)
        end do pad_log_grid
        fill_log_grid: do ipt=sd_rgrid_npad+1,sd_nradial+1
          sd_rtab(ipt) = sd_rgrid_r0 * sd_rgrid_scale**(ipt-1-sd_rgrid_npad)
        end do fill_log_grid
      case ('log-uniform')
        if (sd_rgrid_r0<=sd_rgrid_rmin) then
          write (out,"('spherical_data_initialize%initialize_radial_grid: r0 (',g0.5,') is less than rmin (',g0.5,')')") &
                 sd_rgrid_r0, sd_rgrid_rmin
          stop 'spherical_data_initialize%initialize_radial_grid - log-linear grid starts inside infinite barrier!'
        end if
        if (sd_rgrid_npad<0) then
          sd_rgrid_npad = max(0,nint((1._rk-sd_rgrid_rmin/sd_rgrid_r0)/(sd_rgrid_scale-1._rk)))
          !
          !  Choose master node's result for sd_rgrid_npad. We do not want small numerical differences
          !  between the nodes (which -may- have different architectures) to produce different grids
          !  on different nodes!
          !
          call nt_broadcast(sd_rgrid_npad)
        end if
        npt_log = 1 + nint(log(sd_rgrid_dr/(sd_rgrid_r0*(sd_rgrid_scale-1)))/log(sd_rgrid_scale))
        call nt_broadcast(npt_log)
        write (out,"('Creating log-uniform grid')")
        write (out,"('Logarithmic grid starting at ',g24.13,' Bohr; scale factor ',g24.13)") sd_rgrid_r0, sd_rgrid_scale
        write (out,"('Inserting ',i0,'-point uniform padding grid at the origin')") sd_rgrid_npad
        write (out,"('Uniform grid step ',g24.13,' Bohr')") sd_rgrid_dr
        write (out,"('Logarithmic to uniform switch after point ',i0)") npt_log
        !
        if (npt_log<1 .or. npt_log>=sd_nradial-sd_rgrid_npad) then
          write (out,"('Logarithmic to uniform switching point ',i0,' is not between 1 and ',i0)") npt_log, sd_nradial-sd_rgrid_npad
          stop 'spherical_data_initialize%initialize_radial_grid - bad log-uniform parameters'
        end if
        !
        pad_lu_grid: do ipt=1,sd_rgrid_npad
          sd_rtab(ipt) = sd_rgrid_rmin + (ipt*(sd_rgrid_r0-sd_rgrid_rmin))/(sd_rgrid_npad+1)
        end do pad_lu_grid
        fill_lu_grid_log_part: do ipt=1+sd_rgrid_npad,npt_log+sd_rgrid_npad
          sd_rtab(ipt) = sd_rgrid_r0 * sd_rgrid_scale**(ipt-1-sd_rgrid_npad)
        end do fill_lu_grid_log_part
        fill_lu_grid_uniform_part: do ipt=npt_log+1+sd_rgrid_npad,sd_nradial+1
          sd_rtab(ipt) = sd_rtab(ipt-1) + sd_rgrid_dr
        end do fill_lu_grid_uniform_part
        !
        write (out,"('    Uniform grid: ',i6,' points between ',g24.13,' and ',g24.13,' Bohr')") &
               sd_rgrid_npad, sd_rtab(1), sd_rtab(sd_rgrid_npad)
        write (out,"('Logarithmic grid: ',i6,' points between ',g24.13,' and ',g24.13,' Bohr')") &
               npt_log, sd_rtab(sd_rgrid_npad+1), sd_rtab(sd_rgrid_npad+npt_log)
        write (out,"('    Uniform grid: ',i6,' points between ',g24.13,' and ',g24.13,' Bohr')") &
               sd_nradial-npt_log-sd_rgrid_npad, sd_rtab(npt_log+sd_rgrid_npad+1), sd_rtab(sd_nradial)
      case ('read')
        write (out,"('Reading user-defined from ',a)") trim(sd_rgrid_file)
        open (iu_temp,form='formatted',action='read',position='rewind',status='old',file=trim(sd_rgrid_file))
        fill_user_grid: do ipt=1,sd_nradial+1
          read (iu_temp,*,iostat=ios) sd_rtab(ipt)
          if (ios/=0) then
            write (out,"('Error ',i0,' reading point ',i0,' of the radial grid.')") ios, ipt
            stop 'spherical_data_initialize%initialize_radial_grid - read error'
          end if
          if (sd_rtab(ipt)<=sd_rgrid_rmin) then
            write (out,"('Tabulated grid point ',i0,' is inside infinite barrier: ',g0.5)") ipt, sd_rtab(ipt)
            stop 'spherical_data_initialize%initialize_radial_grid - bad tabular grid'
          end if
        end do fill_user_grid
        close (iu_temp)
    end select
    !
    !  Now calculate the distance table from the grid positions
    ! 
    sd_drtab(1) = sd_rtab(1) - sd_rgrid_rmin
    distance_table: do ipt=2,sd_nradial+1
      if (sd_rtab(ipt)<=sd_rtab(ipt-1)) then
        write (out,"('Radial grid point ',i0,' at r= ',g24.13,' is before the preceeding point at r= ',g24.13,'?!')") &
               ipt, sd_rtab(ipt), sd_rtab(ipt-1)
        stop 'spherical_data_initialize%initialize_radial_grid - radial grid is must be in increasing order'
      end if
      sd_drtab(ipt) = sd_rtab(ipt)-sd_rtab(ipt-1)
    end do distance_table
    !
    write (out,"(/'        Number of radial grid points = ',i0)") sd_nradial
    write (out,"( '                The world starts at r= ',g24.13)") sd_rgrid_rmin
    write (out,"( '    First explicit grid point is at r= ',g24.13)") sd_rtab(1)
    write (out,"( '      Last explicit grid point s at r= ',g24.13)") sd_rtab(sd_nradial)
    write (out,"( '                  The world ends at r= ',g24.13)") sd_rtab(sd_nradial+1)
    write (out,"( 'Effective nuclear charge for L=0 grid= ',g24.13)") sd_rgrid_zeta
    write (out,"()")
  end subroutine initialize_radial_grid
  !
  !  For the derivation of the implicit derivative expressions, see notes in:
  !
  !    hgm-operators-radial-derivatives-2014_Oct_31.pdf
  !    hgm-operators-involving-derivatives-2014_Nov_04.pdf
  !
  subroutine initialize_radial_gradient(repeat)
    logical, intent(in) :: repeat
    integer(ik)         :: ipt
    !
    !  Interior points; general expression.
    !  For sanity, all expressions are computer-generated.
    !
    grad_lx: do ipt=1,sd_nradial
      ! Diagonal
      sd_d1n_lx(ipt,1) = 1/sd_drtab(ipt) - 1/sd_drtab(ipt+1) + (-sd_drtab(ipt) + sd_drtab(ipt+1))/ &
                         (sd_drtab(ipt)**2 + sd_drtab(ipt)*sd_drtab(ipt+1) + sd_drtab(ipt+1)**2)
      sd_m1n_lx(ipt,1) = (sd_drtab(ipt) + sd_drtab(ipt+1))**2/ &
                         (2*(sd_drtab(ipt)**2 + sd_drtab(ipt)*sd_drtab(ipt+1) + sd_drtab(ipt+1)**2))
      if (ipt>=sd_nradial) cycle grad_lx
      ! Sub-diagonal
      sd_d1n_lx(ipt,2) = -((sd_drtab(ipt+2)**2*(2*sd_drtab(ipt+1) + sd_drtab(ipt+2)))/ &
                         (sd_drtab(ipt+1)*(sd_drtab(ipt+1) + sd_drtab(ipt+2))* &
                         (sd_drtab(ipt+1)**2 + sd_drtab(ipt+1)*sd_drtab(ipt+2) + sd_drtab(ipt+2)**2)))
      sd_m1n_lx(ipt,2) = sd_drtab(ipt+2)**2/(2*(sd_drtab(ipt+1)**2 + sd_drtab(ipt+1)*sd_drtab(ipt+2) + sd_drtab(ipt+2)**2))
      ! Super-diagonal
      sd_d1n_lx(ipt,3) = (sd_drtab(ipt)**2*(sd_drtab(ipt) + 2*sd_drtab(ipt+1)))/(sd_drtab(ipt+1)* &
                         (sd_drtab(ipt) + sd_drtab(ipt+1))*(sd_drtab(ipt)**2 + &
                         sd_drtab(ipt)*sd_drtab(ipt+1) + sd_drtab(ipt+1)**2))
      sd_m1n_lx(ipt,3) = sd_drtab(ipt)**2/(2*(sd_drtab(ipt)**2 + sd_drtab(ipt)*sd_drtab(ipt+1) + sd_drtab(ipt+1)**2))
    end do grad_lx
    sd_d1n_l0 = sd_d1n_lx
    sd_m1n_l0 = sd_m1n_lx
    !
    if (sd_rgrid_rmin<=0) then
      select case (sd_operators)
        case default
          write (out,"('spherical_data_initialize%initialize_radial_gradient - sd_operators=',a,' is not recognized')") &
                 trim(sd_operators)
          stop 'spherical_data_initialize%initialize_radial_gradient - bad sd_operators'
        case ('grad common','all common')
          !
          !  Use general-case operators also for L=0
          !
        case ('channel','lap common')
          !
          !  Special case: L=0.
          !
          sd_d1n_l0(1,1) = -(((sd_drtab(1) + sd_drtab(2))**2*(sd_drtab(1) + 2*sd_drtab(2))* &
                           (-6*sd_drtab(1)**2 + sd_drtab(2)**2 + 2*sd_drtab(1)**3*sd_rgrid_zeta + &
                           sd_drtab(1)*sd_drtab(2)*(3 - 2*sd_drtab(2)*sd_rgrid_zeta)))/ &
                           (sd_drtab(1)*sd_drtab(2)*(sd_drtab(1)**2 + sd_drtab(1)*sd_drtab(2) + sd_drtab(2)**2)* &
                           (-6*sd_drtab(1)**2 - 15*sd_drtab(1)*sd_drtab(2) - 8*sd_drtab(2)**2 + &
                           2*sd_drtab(1)*(sd_drtab(1) + sd_drtab(2))*(sd_drtab(1) + 2*sd_drtab(2))*sd_rgrid_zeta)))
          sd_m1n_l0(1,1) = ((sd_drtab(1) + sd_drtab(2))**2*(sd_drtab(1) + 2*sd_drtab(2))* (-3*sd_drtab(1) - sd_drtab(2) + &
                           sd_drtab(1)*(sd_drtab(1) + sd_drtab(2))*sd_rgrid_zeta))/((sd_drtab(1)**2 + sd_drtab(1)*sd_drtab(2) &
                           + sd_drtab(2)**2)* (-6*sd_drtab(1)**2 - 15*sd_drtab(1)*sd_drtab(2) - 8*sd_drtab(2)**2 + &
                           2*sd_drtab(1)*(sd_drtab(1) + sd_drtab(2))*(sd_drtab(1) + 2*sd_drtab(2))*sd_rgrid_zeta))
          sd_m1n_l0(1,3) = (sd_drtab(1)**2*(sd_drtab(1) + 2*sd_drtab(2))*(-3*sd_drtab(1) - 2*sd_drtab(2) + sd_drtab(1)* &
                           (sd_drtab(1) + sd_drtab(2))*sd_rgrid_zeta))/((sd_drtab(1)**2 + sd_drtab(1)*sd_drtab(2) + &
                           sd_drtab(2)**2)* (-6*sd_drtab(1)**2 - 15*sd_drtab(1)*sd_drtab(2) - 8*sd_drtab(2)**2 + 2*sd_drtab(1)* &
                           (sd_drtab(1) + sd_drtab(2))*(sd_drtab(1) + 2*sd_drtab(2))*sd_rgrid_zeta))
      end select
    else
      if (.not.repeat) write (out,"('Applying infinite-barrier boundary conditions for the gradient operator.')")
    end if
    !
    if (.not.repeat) then
      write (out,"('L>=1 channel d1(1,1)= ',g24.13,' m1(1,1)= ',g24.13)") sd_d1n_lx(1,1), sd_m1n_lx(1,1)
      write (out,"(' L=0 channel d1(1,1)= ',g24.13,' m1(1,1)= ',g24.13)") sd_d1n_l0(1,1), sd_m1n_l0(1,1)
    end if
    !
    !  Factor m1 matrices; strictly speaking this is not required for propagation,
    !  since we never take the first radial derivative by itself. However, it may
    !  be used in the gradient accuracy check, as well as for the dipole-velocity 
    !  operator.
    !
    call m3d_decompose(sd_m1n_l0,sd_m1nf_l0,sd_m1np_l0)
    call m3d_decompose(sd_m1n_lx,sd_m1nf_lx,sd_m1np_lx)
    !
    !  We also need transposes and matching inverses
    !
    call m3d_transpose(sd_d1n_l0,sd_d1t_l0)
    call m3d_transpose(sd_d1n_lx,sd_d1t_lx)
    call m3d_transpose(sd_m1n_l0,sd_m1t_l0)
    call m3d_transpose(sd_m1n_lx,sd_m1t_lx)
    !
    call m3d_decompose(sd_m1t_l0,sd_m1tf_l0,sd_m1tp_l0)
    call m3d_decompose(sd_m1t_lx,sd_m1tf_lx,sd_m1tp_lx)
  end subroutine initialize_radial_gradient
  !
  subroutine initialize_radial_laplacian(repeat)
    logical, intent(in) :: repeat
    integer(ik)         :: ipt
    !
    !  Begin with the interior points
    !
    lap_lx: do ipt=1,sd_nradial
      ! Diagonal
      sd_d2n_lx(ipt,1) = -2._rk/(sd_drtab(ipt)*sd_drtab(ipt+1))
      sd_m2n_lx(ipt,1) = (1._rk/6._rk)*(3._rk + sd_drtab(ipt)/sd_drtab(ipt+1) + sd_drtab(ipt+1)/sd_drtab(ipt))
      if (ipt>=sd_nradial) cycle lap_lx
      ! Sub-diagonal
      sd_d2n_lx(ipt,2) = 2._rk/(sd_drtab(ipt+1)*(sd_drtab(ipt+1)+sd_drtab(ipt+2)))
      sd_m2n_lx(ipt,2) = (1._rk/6._rk)*(2._rk - sd_drtab(ipt+2)/sd_drtab(ipt+1) &
                                             - sd_drtab(ipt+1)/(sd_drtab(ipt+1)+sd_drtab(ipt+2)))
      ! Super-diagonal
      sd_d2n_lx(ipt,3) = 2._rk/(sd_drtab(ipt+1)*(sd_drtab(ipt)+sd_drtab(ipt+1)))
      sd_m2n_lx(ipt,3) = (1._rk/6._rk)*(1._rk - sd_drtab(ipt)**2/((sd_drtab(ipt)+sd_drtab(ipt+1))*sd_drtab(ipt+1)))
    end do lap_lx
    !
    if (sd_rgrid_rmin<=0) then
      !
      !  For the Coulomb boundary conditions, the first row is special; first L>0
      !  Expressions are computer-generated to preserve my sanity.
      !
      sd_d2n_lx(1,1) = (-2*(3*sd_drtab(1) - sd_drtab(2))*(sd_drtab(1) + sd_drtab(2))**2)/ &
                       (sd_drtab(1)**3*sd_drtab(2)*(3*sd_drtab(1) + 4*sd_drtab(2)))
      sd_m2n_lx(1,1) = ((sd_drtab(1) + sd_drtab(2))**2*(sd_drtab(1) + 3*sd_drtab(2)))/ &
                       (3*sd_drtab(1)*sd_drtab(2)*(3*sd_drtab(1) + 4*sd_drtab(2)))
      sd_m2n_lx(1,3) = -((sd_drtab(1) - 2*sd_drtab(2))*(sd_drtab(1) + sd_drtab(2)))/ &
                       (3*sd_drtab(2)*(3*sd_drtab(1) + 4*sd_drtab(2)))
      !
      sd_d2n_l0 = sd_d2n_lx
      sd_m2n_l0 = sd_m2n_lx
      !
      select case (sd_operators)
        case default
          write (out,"('spherical_data_initialize%initialize_radial_laplacian - sd_operators=',a,' is not recognized')") &
                 trim(sd_operators)
          stop 'spherical_data_initialize%initialize_radial_laplacian - bad sd_operators'
        case ('all common','lap common')
          !
          !  Use general-case operators also for L=0
          !
        case ('channel','grad common')
          !
          !  Now L=0
          !
          sd_d2n_l0(1,1) = (2*(sd_drtab(1) + sd_drtab(2))*(6*sd_drtab(1) - (3*sd_drtab(1) - sd_drtab(2))* &
                           (sd_drtab(1) + sd_drtab(2))*sd_rgrid_zeta))/ (sd_drtab(1)**2*sd_drtab(2)* &
                           (-6*(sd_drtab(1) + sd_drtab(2)) + sd_drtab(1)*(3*sd_drtab(1) + 4*sd_drtab(2))*sd_rgrid_zeta))
          sd_m2n_l0(1,1) = ((sd_drtab(1) + sd_drtab(2))*(-3*(sd_drtab(1)**2 + 3*sd_drtab(1)*sd_drtab(2) + &
                           sd_drtab(2)**2) + sd_drtab(1)*(sd_drtab(1) + sd_drtab(2))*(sd_drtab(1) + 3*sd_drtab(2))* &
                           sd_rgrid_zeta))/  (3*sd_drtab(1)*sd_drtab(2)*(-6*(sd_drtab(1) + sd_drtab(2)) +  &
                           sd_drtab(1)*(3*sd_drtab(1) + 4*sd_drtab(2))*sd_rgrid_zeta))
          sd_m2n_l0(1,3) = (-3*sd_drtab(2)**2 - sd_drtab(1)**3*sd_rgrid_zeta + sd_drtab(1)**2* &
                           (3 + sd_drtab(2)*sd_rgrid_zeta) + sd_drtab(1)*sd_drtab(2)*(-3 + 2*sd_drtab(2)*sd_rgrid_zeta))/ &
                           (3*sd_drtab(2)*(-6*(sd_drtab(1) + sd_drtab(2)) + &
                           sd_drtab(1)*(3*sd_drtab(1) + 4*sd_drtab(2))*sd_rgrid_zeta))
      end select
    else
      sd_d2n_l0 = sd_d2n_lx
      sd_m2n_l0 = sd_m2n_lx
      if (.not.repeat) write (out,"('Applying infinite-barrier boundary conditions for the laplacian operator.')")
    end if
    !
    if (.not.repeat) then
      write (out,"('L>=1 channel d2(1,1)= ',g24.13,' m2(1,1)= ',g24.13)") sd_d2n_lx(1,1), sd_m2n_lx(1,1)
      write (out,"(' L=0 channel d2(1,1)= ',g24.13,' m2(1,1)= ',g24.13)") sd_d2n_l0(1,1), sd_m2n_l0(1,1)
    end if
    !
    !  Factor m2 matrices
    !
    call m3d_decompose(sd_m2n_l0,sd_m2nf_l0,sd_m2np_l0)
    call m3d_decompose(sd_m2n_lx,sd_m2nf_lx,sd_m2np_lx)
    !
    !  We also need transposes and matching inverses
    !
    call m3d_transpose(sd_d2n_l0,sd_d2t_l0)
    call m3d_transpose(sd_d2n_lx,sd_d2t_lx)
    call m3d_transpose(sd_m2n_l0,sd_m2t_l0)
    call m3d_transpose(sd_m2n_lx,sd_m2t_lx)
    !
    call m3d_decompose(sd_m2t_l0,sd_m2tf_l0,sd_m2tp_l0)
    call m3d_decompose(sd_m2t_lx,sd_m2tf_lx,sd_m2tp_lx)
  end subroutine initialize_radial_laplacian
  !
  subroutine initialize_channel_potentials(repeat)
    logical, intent(in) :: repeat
    integer(ik)         :: lval, ipt
    real(rk)            :: rdelta, rup, rdown, fup, fdown, fuplp, fdownlp
    !
    !  Potential itself, including the centrifugal term
    !
    scan_channels: do lval=0,sd_lmax
      scan_points: do ipt=1,sd_nradial
        sd_pottab(ipt,1,lval) = pt_evaluate_potential(lval,0,sd_rtab(ipt),centrifugal=.true.)
      end do scan_points
    end do scan_channels
    !
    !  Radial gradient of the potential (omitting the centrifugal term)
    !  We'll test for non-local potentials here: our dipole acceleration 
    !  expression is only valid for multiplicative potentials.
    !
    sd_pot_nonlocal = .false.
    gradient_points: do ipt=1,sd_nradial
      rdelta = sd_rgrid_grad_delta * sd_drtab(ipt)
      rup    = sd_rtab(ipt) + rdelta
      rdown  = sd_rtab(ipt) - rdelta
      fup    = pt_evaluate_potential(0,0,rup,  centrifugal=.false.)
      fdown  = pt_evaluate_potential(0,0,rdown,centrifugal=.false.)
      gradient_ltest: do lval=1,sd_lmax
        fuplp           = pt_evaluate_potential(lval,0,rup,  centrifugal=.false.)
        fdownlp         = pt_evaluate_potential(lval,0,rdown,centrifugal=.false.)
        sd_pot_nonlocal = sd_pot_nonlocal .or. (abs(fuplp  -fup  )>10._rk*spacing(abs(fup  ))) &
                                          .or. (abs(fdownlp-fdown)>10._rk*spacing(abs(fdown))) 
      end do gradient_ltest
      sd_dvdrtab(ipt) = (fup-fdown)/(rup-rdown)
    end do gradient_points
    if (sd_pot_nonlocal .and. .not.repeat) then
      write (out,"(/'WARNING: Atomic potential appears to be L-dependent.'/)") 
    end if
    !
    initialize_cap: do ipt=1,sd_nradial
      sd_captab(ipt) = cap_evaluate_potential(sd_rtab(ipt),sd_rtab_full(sd_nradial_max+1))
    end do initialize_cap
    sd_capped    = .false.
    sd_cap_start = sd_nradial + 1
    scan_cap: do ipt=1,sd_nradial
      if (sd_captab(ipt)==(0._rk,0._rk)) cycle scan_cap
      sd_capped    = .true.
      sd_cap_start = ipt
      exit scan_cap
    end do scan_cap
    if (sd_capped .and. .not.repeat) then
      write (out,"(/'Complex absorbing potential starts at grid point ',i0,' r= ',g24.13)") sd_cap_start, sd_rtab(sd_cap_start)
    end if
    !
    if (sd_rgrid_report/=' ' .and. .not.repeat) then
      open(iu_temp,form='formatted',action='write',position='rewind',status='replace',file=trim(sd_rgrid_report))
      write (iu_temp,"(('#',1x,a7,5(1x,a24)))") &
        ' IPT ',  ' R, Bohr ', ' V(L=0), Hartree ', ' dV/dR, Hartree/Bohr ', ' Re[Vcap] ', ' Im[Vcap] ', &
        '-----',  '---------', '-----------------', '---------------------', '----------', '----------'
      report_points: do ipt=1,sd_nradial
        write (iu_temp,"(1x,i8,5(1x,g24.16e3))") ipt, sd_rtab(ipt), sd_pottab(ipt,1,0), sd_dvdrtab(ipt), sd_captab(ipt)
      end do report_points
      close(iu_temp)
    end if
  end subroutine initialize_channel_potentials
end module spherical_data_initialize
