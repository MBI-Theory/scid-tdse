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
!   Time-propagation loop and related routines
!
module spherical_tdse_propagate
  use spherical_tdse_data
  use spherical_tdse_io
  implicit none
  private
  public rcsid_spherical_tdse_propagate
  public propagation
  !
  character(len=clen), save :: rcsid_spherical_tdse_propagate = "$Id: spherical_tdse_propagate.f90,v 1.1 2024/02/13 14:22:14 ps Exp $"
  !
  contains
  !
  !  Calculate the maximum "safe" extent of the wavefunction which current grid can handle
  !  We need to consider:
  !   a) kinetic-energy artifacts near the edge
  !   b) the absorbing boundary
  !   c) position of the iSURFV sensor sphere
  !   d) how much wavefunction could grow in a single time step
  !
  integer(ik) function safe_nr_max()
                    safe_nr_max = sd_nradial - sd_radial_edge
    if (sd_capped)  safe_nr_max = min(sd_cap_start-1,safe_nr_max)
    if (sts_volkov) safe_nr_max = min(sts_ipt_match-1,safe_nr_max)
    !
    safe_nr_max = safe_nr_max - 2*wt_adaptive_r_buffer
  end function safe_nr_max
  !
  !  Largest radial extent of the wavefunction across the ensemble
  !
  integer(ik) function ens_nradial()
    ens_nradial = max(maxval(wfns_l(:)%nradial),maxval(wfns_r(:)%nradial))
  end function ens_nradial
  !
  !  Largest angular extent of the wavefunction across the ensemble
  !
  integer(ik) function ens_lmax()
    ens_lmax = max(maxval(wfns_l(:)%lmax),maxval(wfns_r(:)%lmax))
  end function ens_lmax
  !
  subroutine grid_resize(rmax)
    integer(ik), intent(in) :: rmax ! Current maximum extent of the radial wavefunction in points
    !
    integer(ik) :: nr, ie
    !
    if (.not.sd_adaptive_r) return
    !
    !  Make sure the new radial grid extent is sensible [see comments above safe_nr_max()]
    !
    nr = int(sd_nradial_scl * rmax,kind=ik) + 1_ik     ! Grow by a factor
    nr = nr + sd_radial_edge + 2*wt_adaptive_r_buffer  ! Add safety buffers
    if (sts_volkov .and. nr>=sts_ipt_match) then       ! We'll reach the iSURF sensor surface, switch to the full grid
      nr = sd_nradial_max
    end if
    !
    !  Make sure we don't try to increment grid size by 1 - this adds too much overhead
    !
    if (sd_nradial<sd_nradial_max) then
      nr = max(nr,int(sd_nradial_scl * sd_nradial,kind=ik))
    end if
    !
    nr = min(nr,sd_nradial_max) ! Make sure we do not exceed the maximum
    if (nr==sd_nradial) return  ! Grid is already of the right size, nothing to do here
    !
    if (verbose>=0) then
      write (out,"('Resizing radial wavefunction to ',i8,' grid points')") nr
      call flush_wrapper(out)
    end if
    !
    resize_loop: do ie=1,ensemble_size
      call wt_resize(wfns_l(ie),nr)
      call wt_resize(wfns_r(ie),nr)
    end do resize_loop
    sd_nradial = nr
    call sd_initialize(repeat=.true.)
    call pt_reset_caches
  end subroutine grid_resize
  !
  subroutine merge_all_wavefunctions
    integer(ik) :: ie
    !
    merge_ensemble: do ie=1,ensemble_size
      call nt_merge_all(wfns_r(ie))
      if (.not.skip_left_propagation) then
        call nt_merge_all(wfns_l(ie))
      end if
    end do merge_ensemble
  end subroutine merge_all_wavefunctions
  !
  subroutine checkpoint_all_wavefunctions(its,go)
    integer(ik), intent(in) :: its
    logical, intent(out)    :: go
    integer(ik)             :: ie
    logical                 :: checkpoint_go(ensemble_size)
    !
    checkpoint_loop: do ie=1,ensemble_size
      call ckpt_do_checkpoint(ckpts(ie),checkpoint_go(ie),its,vpot_table(0:3,2*its),(/out,iu_detail(ie),iu_ens/), &
                wfns_l(ie),wfns_r(ie),tsurfs(ie))
    end do checkpoint_loop
    if (any(checkpoint_go) .and. .not.all(checkpoint_go)) then
      write (out,"('ERROR: Checkpoint files are out of sync.')") 
      write (out,"(50(1x,l1))") checkpoint_go
      stop 'spherical_tdse_propagate%checkpoint_all_wavefunctions - checkpoint mismatch'
    end if
    go = all(checkpoint_go)
  end subroutine checkpoint_all_wavefunctions
  !
  subroutine subdivide_timestep(its,lmax,nsub,vsub)
    integer(ik), intent(in)  :: its       ! "Macro" time step
    integer(ik), intent(in)  :: lmax      ! Largest angular momentum on this time step
    integer(ik), intent(out) :: nsub      ! Number of "small" time steps after subdivision 
    real(xk), allocatable    :: vsub(:,:) ! Subdivided time step parameters
    !
    real(xk)                 :: cdt       ! The "big" time (half-)step
    real(xk)                 :: sdt       ! The "small" time (half-)step
    real(xk)                 :: vp        ! Vector-potential at the current time step
    real(xk)                 :: ladt      ! The value of (lmax+1)*vector_potential*dt for the "big" step
    integer(ik)              :: alloc
    integer(ik)              :: isub, ic
    !
    cdt = 0.5_rk*(vpot_table(0,2*its+2)-vpot_table(0,2*its))
    vp  = abs(vpot_table(1,2*its+1))
    !
    select case (dt_subdivision)
      case default
        write (out,"('subdivide_timestep: Unrecognized value of dt_subdivision: ',a)") trim(dt_subdivision)
        stop 'spherical_tdse%subdivide_timestep - bad dt_subdivision'
      case ('off')
        ladt = dt_max_la
      case ('lmax','on')
        ladt = abs((lmax+1)*dt)
      case ('lmax-a')
        ladt = abs((lmax+1)*vp**dt)
      case ('lmax-a2')
        ladt = abs((lmax+1)*vp**2*dt)
    end select
    if (ladt<=dt_max_la) then
      nsub = 1
    else
      nsub = ceiling(ladt/dt_max_la)
    end if
    if (nsub>max_subdivision) then 
      max_subdivision = nsub
      if (verbose>=0) write (out,"('Maximum of the time step subdivision factor increased to ',i0)") max_subdivision
    end if
    if (nsub>1) timesteps_subdivided = timesteps_subdivided + 1
    timesteps_micro = timesteps_micro + nsub
    !
    !  Prepare vector-potential array as needed.
    !
    if (allocated(vsub)) then
      if (ubound(vsub,dim=2)<2*nsub) deallocate (vsub) ! Spurious gfortran warning here
    end if
    if (.not.allocated(vsub)) then
      allocate (vsub(0:3,0:2*nsub),stat=alloc)
      if (alloc/=0) then
        write (out,"('spherical_tdse%subdivide_timestep: allocation for nsub = ',i0,' failed. Error = ',i0)") nsub, alloc
        call flush_wrapper(out)
        stop 'spherical_tdse%subdivide_timestep - alloc'
      end if
    end if
    !
    !  Copy vector-potential points at the beginning and end of the interval.
    !  For the remaining points, we'll need to interpolate.
    !
    sdt = cdt / nsub
    vsub(:,0) = vpot_table(:,2*its)
    fill_vsub: do isub=1,2*nsub-1
      vsub(0,isub) = vsub(0,0) + sdt*isub
      fill_components: do ic=1,3
        vsub(ic,isub) = interpolate_vpot(ic,vsub(0,isub),its)
      end do fill_components
    end do fill_vsub
    vsub(:,2*nsub) = vpot_table(:,2*its+2)
    !
    if (verbose>=3) then
      write (out,"(/(2x,a8,4(1x,a20)))") &
                       ' TS ', ' T ', ' A ', ' Theta ', ' Phi ', &
                       '----', '---', '---', '-------', '-----'
      debug_print_subdivide: do isub=0,2*nsub
        if (isub==0*nsub) write (out,"('R ',i8,4(1x,g20.12e3))") 2*its+0, vpot_table(:,2*its+0)
        if (isub==1*nsub) write (out,"('R ',i8,4(1x,g20.12e3))") 2*its+1, vpot_table(:,2*its+1)
        if (isub==2*nsub) write (out,"('R ',i8,4(1x,g20.12e3))") 2*its+2, vpot_table(:,2*its+2)
                          write (out,"('S ',i8,4(1x,g20.12e3))") isub, vsub(:,isub)
      end do debug_print_subdivide
      write (out,"()")
    end if
    !
    contains
    !
    function interpolate_vpot(ic,tim,its) result(v)
      integer(ik), intent(in) :: ic  ! Component of vpot_table to interpolate: vpot_table(ic,:)
      real(xk), intent(in)    :: tim ! Time point to interpolate at
      integer(ik), intent(in) :: its ! Reference point in vpot_table: vpot_table(:,its)
      real(xk)                :: v   ! Interpolated value
      !
      integer(ik) :: lp, rp ! Left and right points in vpot_table
      !
      lp = max(lbound(vpot_table,dim=2),2*its-dt_interpolant_width/2)
      rp = min(ubound(vpot_table,dim=2),lp+dt_interpolant_width-1)
      v  = MathInterpolate(tim,vpot_table(0,lp:rp),vpot_table(ic,lp:rp))
    end function interpolate_vpot
  end subroutine subdivide_timestep
  !
  subroutine rotate_wavefunction(from,to,wfn_l,wfn_r)
    real(xk), intent(in)        :: from(:), to(:) ! Initial and final field orientation
    type(sd_wfn), intent(inout) :: wfn_l          ! Left wavefunction 
    type(sd_wfn), intent(inout) :: wfn_r          ! Right wavefunction 
    !
    select case (rotation_mode)
      case default
        write (out,"('spherical_tdse%propagation: Rotation mode ',a,' is not recognized')") trim(rotation_mode)
        stop 'spherical_tdse%propagation - bad rotation_mode'
      case ('none')
      case ('sparse')
        call rt_rotate(wfn_r,from=from,to=to,left=.false.)
        if (.not.skip_left_propagation) then
          call rt_rotate(wfn_l,from=from,to=to,left=.true.)
        end if
      case ('brute force')
        call rt_rotate_bruteforce(wfn_r,from=from,to=to,left=.false.)
        if (.not.skip_left_propagation) then
          call rt_rotate_bruteforce(wfn_l,from=from,to=to,left=.true.)
        end if
    end select
  end subroutine rotate_wavefunction
  !
  subroutine lab_dipole(wfn_l,wfn_r,rot,afield,efield,norm,dipole,velocity,acceleration,off_diagonal)
    type(sd_wfn), intent(in) :: wfn_l           ! Left wavefunction 
    type(sd_wfn), intent(in) :: wfn_r           ! Right wavefunction 
    real(rk), intent(in)     :: rot(3,3)        ! Rotation matrix (lab->local)
    real(rk), intent(in)     :: afield(3)       ! Vector-potential, needed for the dipole velocity
    real(rk), intent(in)     :: efield(3)       ! Electric field, needed for the dipole acceleration
    complex(rk), intent(in)  :: norm            ! Wavefunction norm, ditto
    complex(rk), intent(out) :: dipole(3)       ! <L|q r|R> expectation value in the lab frame
    complex(rk), intent(out) :: velocity(3)     ! (d/d t) <L|q r|R> in the lab frame
    complex(rk), intent(out) :: acceleration(3) ! (d^2/d t^2) <L|q r|R> in the lab frame
    logical, intent(in)      :: off_diagonal    ! .True. if <L| and |R> are not the same W.F.
                                                ! In this case, plasma correction must be disabled
    !
    complex(rk) :: loc_dipole(3,3)     ! Dipole and its derivatives in the local coordinate system
    !
    call wt_dipole(wfn_l,wfn_r,do_dipole,loc_dipole)
    !
    !  Rotate everything back into the lab system
    !
    !  gfortran will generate a temporary array below, promotiong rot to complex(rk).
    !  that matmul implementation really sucks ...
    ! 
    dipole       = matmul(transpose(rot),loc_dipole(:,1))
    velocity     = matmul(transpose(rot),loc_dipole(:,2))
    acceleration = matmul(transpose(rot),loc_dipole(:,3))
    !
    !  Add operator-derivative terms
    !
    if (do_dipole_plasma .and. .not.off_diagonal) then
      if (do_dipole(2)) velocity     = velocity     - (electron_charge**2/electron_mass) * afield
      if (do_dipole(3)) acceleration = acceleration + (electron_charge**2/electron_mass) * efield
    else
      if (do_dipole(2)) velocity     = velocity     - (electron_charge**2/electron_mass) * afield * norm
      if (do_dipole(3)) acceleration = acceleration + (electron_charge**2/electron_mass) * efield * norm
    end if
  end subroutine lab_dipole
  !
  subroutine write_detail_headers
    integer(ik) :: ie
    !
    headers: do ie=1,ensemble_size
      if (ens_detail_output(ie)==' ') cycle headers
      if (ensemble_size==1) then
        write (out,"(/'Saving detailed output to ',a/)") trim(ens_detail_output(ie))
      else
        write (out,"(/'Saving detailed output for w.f. ',i0,' to ',a/)") ie, trim(ens_detail_output(ie))
      end if
      open(newunit=iu_detail(ie),file=trim(ens_detail_output(ie)),form='formatted',status='replace', &
           position='rewind',recl=1050,pad='no')
      write (iu_detail(ie),"('# Field Columns Data')")
      write (iu_detail(ie),"('#  1    2,13    Timestep')")
      write (iu_detail(ie),"('#  2   15,46    Time, au[t]')")
      write (iu_detail(ie),"('#  3   48,79    Vector-potential magnitude')")
      write (iu_detail(ie),"('#  4   81,112   Vector-potential, lab theta')")
      write (iu_detail(ie),"('#  5  114,145   Vector-potential, lab phi')")
      write (iu_detail(ie),"('#  6  147,178   Re[<Psi_L|Psi_R>]')")
      write (iu_detail(ie),"('#  7  180,211   Im[<Psi_L|Psi_R>]')")
      write (iu_detail(ie),"('#  8  213,244   Re[<Psi_L|H_at+H_L+V_cap|Psi_R>], Hartree')")
      write (iu_detail(ie),"('#  9  246,277   Im[<Psi_L|H_at+H_L+V_cap|Psi_R>], Hartree')")
      write (iu_detail(ie),"('# 10  279,310   Re[<Psi_L|H_at+H_L|Psi_R>], Hartree')")
      write (iu_detail(ie),"('# 11  312,343   Im[<Psi_L|H_at+H_L|Psi_R>], Hartree')")
      write (iu_detail(ie),"('# 12  345,376   Re[<Psi_L|e . x|Psi_R>], e-Bohr')")
      write (iu_detail(ie),"('# 13  378,409   Im[<Psi_L|e . x|Psi_R>], e-Bohr')")
      write (iu_detail(ie),"('# 14  411,442   Re[<Psi_L|e . y|Psi_R>], e-Bohr')")
      write (iu_detail(ie),"('# 15  444,475   Im[<Psi_L|e . y|Psi_R>], e-Bohr')")
      write (iu_detail(ie),"('# 16  477,508   Re[<Psi_L|e . z|Psi_R>], e-Bohr')")
      write (iu_detail(ie),"('# 17  510,541   Im[<Psi_L|e . z|Psi_R>], e-Bohr')")
      write (iu_detail(ie),"('# 18  543,574   Re[(d^2/d t^2) <Psi_L|e . x|Psi_R>], e-Bohr/jiffy^2')")
      write (iu_detail(ie),"('# 19  576,607   Im[(d^2/d t^2) <Psi_L|e . x|Psi_R>], e-Bohr/jiffy^2')")
      write (iu_detail(ie),"('# 20  609,640   Re[(d^2/d t^2) <Psi_L|e . y|Psi_R>], e-Bohr/jiffy^2')")
      write (iu_detail(ie),"('# 21  642,673   Im[(d^2/d t^2) <Psi_L|e . y|Psi_R>], e-Bohr/jiffy^2')")
      write (iu_detail(ie),"('# 22  675,706   Re[(d^2/d t^2) <Psi_L|e . z|Psi_R>], e-Bohr/jiffy^2')")
      write (iu_detail(ie),"('# 23  708,739   Im[(d^2/d t^2) <Psi_L|e . z|Psi_R>], e-Bohr/jiffy^2')")
      write (iu_detail(ie),"('# 24  741,772   Re[(d/d t) <Psi_L|e . x|Psi_R>], e-Bohr/jiffy')")
      write (iu_detail(ie),"('# 25  774,805   Im[(d/d t) <Psi_L|e . x|Psi_R>], e-Bohr/jiffy')")
      write (iu_detail(ie),"('# 26  807,838   Re[(d/d t) <Psi_L|e . y|Psi_R>], e-Bohr/jiffy')")
      write (iu_detail(ie),"('# 27  840,871   Im[(d/d t) <Psi_L|e . y|Psi_R>], e-Bohr/jiffy')")
      write (iu_detail(ie),"('# 28  873,904   Re[(d/d t) <Psi_L|e . z|Psi_R>], e-Bohr/jiffy')")
      write (iu_detail(ie),"('# 29  906,937   Im[(d/d t) <Psi_L|e . z|Psi_R>], e-Bohr/jiffy')")
      write (iu_detail(ie),"('# 30  939,970   Electric field, lab X, atomic units')")
      write (iu_detail(ie),"('# 31  972,1003  Electric field, lab Y, atomic units')")
      write (iu_detail(ie),"('# 32 1005,1036  Electric field, lab Z, atomic units')")
      if (sd_pot_nonlocal) then
        write (iu_detail(ie),"('# WARNING: Non-local potential detected. Dipole acceleration results are incorrect')")
        write (out,      "(/ 'WARNING: Non-local potential detected. Dipole acceleration results are incorrect'/)")
      end if
      write (iu_detail(ie),"('#',1x,a12,31(1x,a32))") &
             ' i ', ' time ', ' vp ', ' theta ', ' phi ', ' re(norm) ', ' im(norm) ', &
             ' re(energy) ', ' im(energy) ', ' re(en-nocap) ', ' im(en-nocap) ', &
             ' re(dip_x) ', ' im(dip_x) ', ' re(dip_y) ', ' im(dip_y) ', ' re(dip_z) ', ' im(dip_z) ', &
             ' re(acc_x) ', ' im(acc_x) ', ' re(acc_y) ', ' im(acc_y) ', ' re(acc_z) ', ' im(acc_z) ', &
             ' re(vel_x) ', ' im(vel_x) ', ' re(vel_y) ', ' im(vel_y) ', ' re(vel_z) ', ' im(vel_z) ', &
             ' F_x ', ' F_y ', ' F_z '
      write (iu_detail(ie),"('#',1x,a12,31(1x,a32))") &
             ' 1 ', ' 2 ', ' 3 ', ' 4 ', ' 5 ', ' 6 ', ' 7 ', ' 8 ', ' 9 ', ' 10 ', ' 11 ', &
             ' 12 ', ' 13 ', ' 14 ', ' 15 ', ' 16 ', ' 17 ', &
             ' 18 ', ' 19 ', ' 20 ', ' 21 ', ' 22 ', ' 23 ', &
             ' 24 ', ' 25 ', ' 26 ', ' 27 ', ' 28 ', ' 29 ', &
             ' 30 ', ' 31 ', ' 32 '
    end do headers
  end subroutine write_detail_headers
  !
  subroutine write_ensemble_headers
    if (ensemble_output==' ' .or. ensemble_size<=1) return
    write (out,"(/'Saving ensemble output to ',a/)") trim(ensemble_output)
    open(newunit=iu_ens,file=trim(ensemble_output),form='formatted',status='replace', &
         position='rewind',recl=1080,pad='no')
    write (iu_ens,"('# Field Columns Data')")
    write (iu_ens,"('#  1    2,13    Timestep')")
    write (iu_ens,"('#  2   15,19    <Psi_L| index within the ensemble')")
    write (iu_ens,"('#  3   21,25    |Psi_R> index within the ensemble')")
    write (iu_ens,"('#  4   27,58    Time, au[t]')")
    write (iu_ens,"('#  5   60,91    Vector-potential magnitude')")
    write (iu_ens,"('#  6   93,124   Vector-potential, lab theta')")
    write (iu_ens,"('#  7  126,157   Vector-potential, lab phi')")
    write (iu_ens,"('#  8  159,190   Re[<Psi_L|Psi_R>]')")
    write (iu_ens,"('#  9  192,223   Im[<Psi_L|Psi_R>]')")
    write (iu_ens,"('# 10  225,256   Re[<Psi_L|H_at+H_L+V_cap|Psi_R>], Hartree')")
    write (iu_ens,"('# 11  258,289   Im[<Psi_L|H_at+H_L+V_cap|Psi_R>], Hartree')")
    write (iu_ens,"('# 12  291,322   Re[<Psi_L|H_at+H_L|Psi_R>], Hartree')")
    write (iu_ens,"('# 13  324,355   Im[<Psi_L|H_at+H_L|Psi_R>], Hartree')")
    write (iu_ens,"('# 14  357,388   Re[<Psi_L|e . x|Psi_R>], e-Bohr')")
    write (iu_ens,"('# 15  390,421   Im[<Psi_L|e . x|Psi_R>], e-Bohr')")
    write (iu_ens,"('# 16  423,454   Re[<Psi_L|e . y|Psi_R>], e-Bohr')")
    write (iu_ens,"('# 17  456,487   Im[<Psi_L|e . y|Psi_R>], e-Bohr')")
    write (iu_ens,"('# 18  489,520   Re[<Psi_L|e . z|Psi_R>], e-Bohr')")
    write (iu_ens,"('# 19  522,553   Im[<Psi_L|e . z|Psi_R>], e-Bohr')")
    write (iu_ens,"('# 20  555,586   Re[(d^2/d t^2) <Psi_L|e . x|Psi_R>], e-Bohr/jiffy^2')")
    write (iu_ens,"('# 21  588,619   Im[(d^2/d t^2) <Psi_L|e . x|Psi_R>], e-Bohr/jiffy^2')")
    write (iu_ens,"('# 22  621,652   Re[(d^2/d t^2) <Psi_L|e . y|Psi_R>], e-Bohr/jiffy^2')")
    write (iu_ens,"('# 23  654,685   Im[(d^2/d t^2) <Psi_L|e . y|Psi_R>], e-Bohr/jiffy^2')")
    write (iu_ens,"('# 24  687,718   Re[(d^2/d t^2) <Psi_L|e . z|Psi_R>], e-Bohr/jiffy^2')")
    write (iu_ens,"('# 25  720,751   Im[(d^2/d t^2) <Psi_L|e . z|Psi_R>], e-Bohr/jiffy^2')")
    write (iu_ens,"('# 26  753,784   Re[(d/d t) <Psi_L|e . x|Psi_R>], e-Bohr/jiffy')")
    write (iu_ens,"('# 27  786,817   Im[(d/d t) <Psi_L|e . x|Psi_R>], e-Bohr/jiffy')")
    write (iu_ens,"('# 28  819,850   Re[(d/d t) <Psi_L|e . y|Psi_R>], e-Bohr/jiffy')")
    write (iu_ens,"('# 29  852,883   Im[(d/d t) <Psi_L|e . y|Psi_R>], e-Bohr/jiffy')")
    write (iu_ens,"('# 30  885,916   Re[(d/d t) <Psi_L|e . z|Psi_R>], e-Bohr/jiffy')")
    write (iu_ens,"('# 31  918,949   Im[(d/d t) <Psi_L|e . z|Psi_R>], e-Bohr/jiffy')")
    write (iu_ens,"('# 32  951,982   Electric field, lab X, atomic units')")
    write (iu_ens,"('# 33  984,1015  Electric field, lab Y, atomic units')")
    write (iu_ens,"('# 34 1017,1048  Electric field, lab Z, atomic units')")
    if (sd_pot_nonlocal) then
      write (iu_ens,"('# WARNING: Non-local potential detected. Dipole acceleration results are incorrect')")
      write (out,      "(/ 'WARNING: Non-local potential detected. Dipole acceleration results are incorrect'/)")
    end if
    write (iu_ens,"('#',1x,a12,2(1x,a5),31(1x,a32))") &
           ' i ', ' <l| ', ' |r> ', ' time ', ' vp ', ' theta ', ' phi ', ' re(norm) ', ' im(norm) ', &
           ' re(energy) ', ' im(energy) ', ' re(en-nocap) ', ' im(en-nocap) ', &
           ' re(dip_x) ', ' im(dip_x) ', ' re(dip_y) ', ' im(dip_y) ', ' re(dip_z) ', ' im(dip_z) ', &
           ' re(acc_x) ', ' im(acc_x) ', ' re(acc_y) ', ' im(acc_y) ', ' re(acc_z) ', ' im(acc_z) ', &
           ' re(vel_x) ', ' im(vel_x) ', ' re(vel_y) ', ' im(vel_y) ', ' re(vel_z) ', ' im(vel_z) ', &
           ' F_x ', ' F_y ', ' F_z '
    write (iu_ens,"('#',1x,a12,2(1x,a5),31(1x,a32))") &
           ' 1 ', ' 2 ', ' 3 ', ' 4 ', ' 5 ', ' 6 ', ' 7 ', ' 8 ', ' 9 ', ' 10 ', ' 11 ', &
           ' 12 ', ' 13 ', ' 14 ', ' 15 ', ' 16 ', ' 17 ', &
           ' 18 ', ' 19 ', ' 20 ', ' 21 ', ' 22 ', ' 23 ', &
           ' 24 ', ' 25 ', ' 26 ', ' 27 ', ' 28 ', ' 29 ', &
           ' 30 ', ' 31 ', ' 32 ', ' 33 ', ' 34 '
  end subroutine write_ensemble_headers
  !
  subroutine propagation
    integer(ik)           :: its, ie           ! Time step and ensemble member
    integer(ik)           :: its_from          ! Time step (possibly half-step) from which we want to step with t-SURF
    integer(ik)           :: its_to            ! Time step (possibly half-step) to which we want to step with t-SURF
    integer(ik)           :: nsub, isub        ! Number of micro time steps after subdivision, and the counter
    real(xk), allocatable :: vsub(:,:)         ! Same as vpot_table(), but for the subdivided time step.
                                               ! The first index in 0=time, 1=signed magnitude, 2=theta, 3=phi
                                               ! The last index runs from 0 to 2*nsub.
                                               ! The array is maintaned by subdivide_timestep
    real(xk)              :: time              ! Time at the beginning of the macro time step
    real(rk)              :: afield(3)         ! The vector potential
    real(rk)              :: efield(3)         ! The electric field
    real(xk)              :: rdt, rdt2
    real(xk)              :: vp1, vp2          ! Magnitude of the vector-potential for the first and second half-steps
    real(xk)              :: th1, th2          ! vector-potential direction, ditto
    real(xk)              :: ph1, ph2          ! vector-potential direction, ditto
    complex(xk)           :: cdt1, cdt2        ! Time steps for the first and second half-steps
    logical               :: checkpoint_go
    !
    !  Propagation routine can handle varying time steps, at least in principle
    !
    if (all(ens_detail_output==' ') .and. (ensemble_output==' '.or.ensemble_size<=1)) then
      if (any(do_dipole(2:3))) then
        write (out,"(/'Detailed output is not enabled; will not calculate dipole velocity or acceleration'/)")
      end if
      do_dipole(2:3) = .false.
    end if
    if (skip_left_propagation) then
      if (fake_left_propagation) then
        write (out,"(/'WARNING: Left wavefunction approximated by conjugate of the right wavefunction.')")
        write (out,"( 'WARNING: Expect reduced accuracy for the observables during propagation.'/)")
      else
        write (out,"(/'WARNING: Left wavefunction is not evaluated. Observables are incorrect.'/)")
      end if
    end if
    if (nts%this_node==1) then
      call write_detail_headers
      call write_ensemble_headers
    end if
    !
    !  Upon entry, we have no idea how the grids are sized, and whether things are balanced.
    !
    call grid_resize(ens_nradial())
    call merge_all_wavefunctions
    call nt_rebalance(ens_lmax())
    !
    time_steps: do its=0,timesteps-1
      time   = vpot_table(0,2*its)
      !
      !  vector-potential at the beginning of the macro step
      !
      vp1    = vpot_table(1,2*its)
      th1    = vpot_table(2,2*its)
      ph1    = vpot_table(3,2*its)
      !
      call checkpoint_all_wavefunctions(its,checkpoint_go)
      if (.not.checkpoint_go) cycle time_steps ! Restart has been requested, but we are not at the point
                                               ! indicated in the checkpoint file yet.
      !
      !  Cartesian vector-potential and electric field are at the beginning of the time step
      !
      afield = real(vp1*(/sin(th1)*cos(ph1),sin(th1)*sin(ph1),cos(th1)/),kind=kind(afield))
      efield = real(efield_table(:,2*its),kind=kind(efield))
      call write_output(force=.false.)  ! Internal subroutine
      call visualize   (force=.false.)  ! Internal subroutine
      !
      !  See whether this time step needs to be subdivided to keep things stable
      !
      call subdivide_timestep(its,ens_lmax(),nsub,vsub)
      !
      !  The subdivided timestep.
      !
      subdivided_timestep: do isub=0,nsub-1
        cdt1 = 0.5_xk*(vsub(0,2*isub+1) - vsub(0,2*isub+0))
        cdt2 = 0.5_xk*(vsub(0,2*isub+2) - vsub(0,2*isub+1))
        !  vector-potential at the beginning of the micro step
        vp1    = vsub(1,2*isub)
        th1    = vsub(2,2*isub)
        ph1    = vsub(3,2*isub)
        !  vector-potential at the end of the micro step
        vp2    = vsub(1,2*isub+2)
        th2    = vsub(2,2*isub+2)
        ph2    = vsub(3,2*isub+2)
        !
        !  Each function within the ensemble is independent; we can propagate them
        !  one at a time, improving cache locality. We only need to sync the time
        !  at the end of the full time step.
        !
        ensemble_loop: do ie=1,ensemble_size
          call nt_merge_borders(wfns_r(ie))
          call pt_fwd_laser_n(wfns_r(ie),vp1,cdt1)
          if (.not.skip_left_propagation) then
            call nt_merge_borders(wfns_l(ie))
            call pt_fwd_laser_t(wfns_l(ie),vp1,cdt1)
          end if
          !
          !  We are at the mid-step; apply atomic propagator for the first half-step
          !
          call pt_fwd_atomic_n(wfns_r(ie),cdt1)
          call pt_rev_atomic_n(wfns_r(ie),cdt1)
          if (.not.skip_left_propagation) then
            call pt_fwd_atomic_t(wfns_l(ie),cdt1)
            call pt_rev_atomic_t(wfns_l(ie),cdt1)
          end if
          !
          !  Rotation; nominally a full step. It will be broken into still smaller sub-steps
          !  if this turns out to be necessary.
          !
          call rotate_wavefunction(from=(/th1,ph1/),to=(/th2,ph2/),wfn_l=wfns_l(ie),wfn_r=wfns_r(ie))
          !
          !  Atomic propagator for the second half-step
          !
          call pt_fwd_atomic_n(wfns_r(ie),cdt2)
          call pt_rev_atomic_n(wfns_r(ie),cdt2)
          if (.not.skip_left_propagation) then
            call pt_fwd_atomic_t(wfns_l(ie),cdt2)
            call pt_rev_atomic_t(wfns_l(ie),cdt2)
          end if
          !
          !  Second half of the time step; vector potential is at the <i>next</i> time-step
          !
          call nt_merge_borders(wfns_r(ie))
          call pt_rev_laser_n(wfns_r(ie),vp2,cdt2)
          if (.not.skip_left_propagation) then
            call nt_merge_borders(wfns_l(ie))
            call pt_rev_laser_t(wfns_l(ie),vp2,cdt2)
          end if
        end do ensemble_loop
        !
        call grid_adapt ! Internal subroutine here
        !
      end do subdivided_timestep
      !
      !  Node rebalancing is potentially an expensive operation!
      !
      if (nt_rebalance_needed(ens_lmax())) then
        call merge_all_wavefunctions
        call nt_rebalance(ens_lmax())
      end if
      !
      !  Photoelectron spectrum accumulation with t-SURF. 
      !
      !  We do not necessarily wish to run t-SURF at each time step (it is expensive,
      !  compared to the propagator). However, we also do not want to overshoot the
      !  end of the simulation. As the result, we need to do a bit of fiddling with
      !  the effective t-SURF time step here.
      !
      !  The test below is deliberately for the time step index "its", not (its+1)
      !  where the actual time step will occur. We wish to perform (an asymmetric) 
      !  t-SURF time step at the very beginning of the simulation as well.
      !
      if (modulo(its,sts_step_frequency)==0) then
         !
         !  We are good to do a t-SURF time step centered at the (full) time step (its+1)
         !  The corresponding time and vector potential are at the offset (2*its+2) in
         !  vpot_table. 
         !
         its_from = max(          0,(2*its+2)-sts_step_frequency)  ! Don't step beyond beginning of the simulation
         its_to   = min(2*timesteps,(2*its+2)+sts_step_frequency)  ! Don't step beyond end of the simulation
         !
         !  If this is the last time step to trigger t-SURF, there is a possibility of missing 
         !  a partial time-step at the very end. Extend the time step to the end of the simulation.
         !
         if (its+sts_step_frequency>timesteps-1) its_to = 2*timesteps
         !
         !  Now we can calculate the duration of the t-SURF time step, and finally do it. 
         !  Because we try to use large time steps in t-SURF, it is not a good idea to assume
         !  that vector-potential remains constant; thankfully, we know the time derivative
         !  of the vector-potential (the electric field; except for the sign change).
         !
         rdt  = vpot_table(0,2*its+2)-vpot_table(0,its_from)  ! Time step to the expansion point
         rdt2 = vpot_table(0,its_to) -vpot_table(0,2*its+2)   ! Time step to the end of the interval
         surf_loop: do ie=1,ensemble_size
           call sts_timestep(wfns_l(ie),wfns_r(ie),tsurfs(ie),rdt,rdt2, &
                             vpot_table(1:3,2*its+2),efield_table(1:3,2*its+2))
         end do surf_loop
      end if
    end do time_steps
    if (allocated(vsub)) deallocate (vsub)
    !
    time   = vpot_table(0,2*its)
    vp1    = vpot_table(1,2*its)
    th1    = vpot_table(2,2*its)
    ph1    = vpot_table(3,2*its)
    afield = real(vp1*(/sin(th1)*cos(ph1),sin(th1)*sin(ph1),cos(th1)/),kind=kind(afield))
    efield = real(efield_table(:,2*its),kind=kind(efield))
    !
    call write_output(force=.true.)  ! Internal subroutine
    call visualize   (force=.true.)  ! Internal subroutine
    !
    close_detail_loop: do ie=1,ensemble_size
      if (ens_detail_output(ie)/=' ' .and. nts%this_node==1) then
        close(iu_detail(ie))
      end if
    end do close_detail_loop
    if (ensemble_output/=' ' .and. ensemble_size>1 .and. nts%this_node==1) then
      close(iu_ens)
    end if
    !
    !  Rotate wavefunction back into the lab orientation for analysis
    !
    rotate_at_end: do ie=1,ensemble_size
      call rotate_wavefunction(from=vpot_table(2:3,2*its),to=(/0._xk,0._xk/),wfn_l=wfns_l(ie),wfn_r=wfns_r(ie))
    end do rotate_at_end
    if (rotation_mode/='none') then
      write (out,"(/'Wavefunction restored back to laboratory frame'/)")
    end if
    !
    call grid_resize(sd_nradial_max)
    !
    contains 
    !
    subroutine write_output(force)
      logical, intent(in) :: force
      !
      integer(ik)         :: ie, je
      real(rk)            :: abs_dipole
      complex(rk)         :: energy      (2,ensemble_size,ensemble_size)
      complex(rk)         :: norm        (  ensemble_size,ensemble_size)
      complex(rk)         :: dipole      (3,ensemble_size,ensemble_size)
      complex(rk)         :: velocity    (3,ensemble_size,ensemble_size)
      complex(rk)         :: acceleration(3,ensemble_size,ensemble_size)
      real(rk)            :: rm(3,3)
      logical             :: normal_out     ! Produce output on the standard output
      logical             :: detail_out     ! Produce output on the detailed output
      logical             :: ens_out        ! Produce output on the ensemble output
      !
      normal_out = force .or. mod(its,output_each)==0
      detail_out = (detail_output  /=' ') .and. (force .or. mod(its,detail_frequency)==0)
      ens_out    = (ensemble_output/=' ') .and. (force .or. mod(its,detail_frequency)==0) .and. ensemble_size>1
      if (.not.normal_out .and. .not.detail_out .and. .not.ens_out) return
      !
      call MathRotationMatrix(real((/ph1,th1,0._xk/),kind=kind(rm)),rm)
      compute_diagonal: do ie=1,ensemble_size
        ! The diagonal part
        if (skip_left_propagation .and. fake_left_propagation) then
          call wt_fake_left(wfns_l(ie),wfns_r(ie))
        end if
        call nt_merge_borders(wfns_l(ie))
        call nt_merge_borders(wfns_r(ie))
        call wt_energy(wfns_l(ie),wfns_r(ie),vp1,energy(:,ie,ie),norm(ie,ie))
        call lab_dipole(wfns_l(ie),wfns_r(ie),rm,afield,efield,norm(ie,ie), &
                        dipole(:,ie,ie),velocity(:,ie,ie),acceleration(:,ie,ie),off_diagonal=.false.) 
      end do compute_diagonal
      if (ens_out) then
        compute_offdiagonal_ie: do ie=1,ensemble_size
          compute_offdiagonal_je: do je=1,ensemble_size
            if (je==ie) cycle compute_offdiagonal_je
            call wt_energy(wfns_l(ie),wfns_r(je),vp1,energy(:,ie,je),norm(ie,je))
            call lab_dipole(wfns_l(ie),wfns_r(je),rm,afield,efield,norm(ie,je), &
                            dipole(:,ie,je),velocity(:,ie,je),acceleration(:,ie,je),off_diagonal=.true.) 
          end do compute_offdiagonal_je
        end do compute_offdiagonal_ie
      end if
      !
      report_diagonal: do ie=1,ensemble_size
        if (detail_out .and. nts%this_node==1) then
          write (iu_detail(ie),"(1x,i12,31(1x,g32.22e4))") &
                 its, time, vp1, th1, ph1, norm(ie,ie), energy(:,ie,ie), dipole(:,ie,ie), &
                 acceleration(:,ie,ie), velocity(:,ie,ie), efield
        end if
        !
        abs_dipole = sqrt(sum(abs(dipole(:,ie,ie))**2))
        if (normal_out) then
          if (ensemble_size<=1) then
            write (out,"('@ ',i9,' t= ',f12.4,' a= ',f14.6,' th= ',f10.3,' ph= ',f10.3,' <l|r>= ',g23.12,1x,g17.6," // &
                       "' <l|h|r>= ',g23.12,1x,g17.6,' |d|= ',g20.9)") &
                   its, time, vp1, th1, ph1, norm(ie,ie), energy(1,ie,ie), abs_dipole
          else
            write (out,"('@w=',i0,' ',i9,' t= ',f12.4,' a= ',f14.6,' th= ',f10.3,' ph= ',f10.3,' <l|r>= ',g23.12,1x,g17.6," // &
                       "' <l|h|r>= ',g23.12,1x,g17.6,' |d|= ',g20.9)") &
                   ie, its, time, vp1, th1, ph1, norm(ie,ie), energy(1,ie,ie), abs_dipole
          end if
          call flush_wrapper(out)
        end if
      end do report_diagonal
      if (ens_out) then
        report_offdiagonal_ie: do ie=1,ensemble_size
          report_offdiagonal_je: do je=1,ensemble_size
            write (iu_ens,"(1x,i12,2(1x,i5),31(1x,g32.22e4))") &
                   its, ie, je, time, vp1, th1, ph1, norm(ie,je), energy(:,ie,je), dipole(:,ie,je), &
                   acceleration(:,ie,je), velocity(:,ie,je), efield
            end do report_offdiagonal_je
        end do report_offdiagonal_ie
      end if
    end subroutine write_output
    !
    subroutine visualize(force)
      logical, intent(in) :: force ! Force output
      integer(ik)         :: ie
      !
      if (mod(its,visualize_each)/=0 .and. .not.force) return
      !
      visualize_loop: do ie=1,ensemble_size
        call visualize_wavefunctions(wfns_l(ie),wfns_r(ie),tsurfs(ie), &
                                     ens_visualize_prefix(ie),ens_visualize_1stseq(ie),its,th1,ph1)
      end do visualize_loop
    end subroutine visualize
    !
    subroutine grid_adapt
      integer(ik) :: ie
      integer(ik) :: rmax  ! Current extent of the left/right wavefunctions
      !
      if (.not.sd_adaptive) return
      !
      !  Update radial extent and LMAX
      !
      radial_update: do ie=1,ensemble_size
        call wt_update_lrmax(verbose,wfns_r(ie))
        if (.not.skip_left_propagation) then
          call wt_update_lrmax(verbose,wfns_l(ie))
        end if
      end do radial_update
      !
      !  If radial grid is already at maximum, there is nothing to do.
      !
      if (sd_nradial==sd_nradial_max) return
      !
      !  If we are still within the safe range, there is nothing to do
      !
      rmax = ens_nradial()
      if (rmax<safe_nr_max()) return
      !
      !  We need to resize now!
      !
      call grid_resize(rmax)
    end subroutine grid_adapt
  end subroutine propagation
end module spherical_tdse_propagate
