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
!  Checkpointing of long-running calculations.
!
!  The checkpoint will attempt do to the following:
!  1) It will flush all currently open files, bringing them to date
!  2) It will write out the essential state of the calculation, which can not be easily recovered from
!     the output files.
!  3) If the number of checkpoints is growing too large, it will delete some of the older checkpoints
!
!  The checkpoint file is binary; it is not intended to be transportable between architectures or 
!  compiler versions, nor for long-term preservation of the results.
!
!  The checkpoint code will NOT try to preserve or restore any external files; this is up to the user.
!
module checkpoint_tools
  use accuracy
  use constants
  use spherical_data
  use spherical_data_initialize
  use spherical_tsurf_data
  use wavefunction_tools
  use node_tools
  use timer
  implicit none
  private
  ! Entry points
  public ckpt_do_checkpoint
  ! Exposed global data
  public rcsid_checkpoint_tools
  public ckpt_data
  public ckpt_save_basename, ckpt_load_filename, ckpt_max_checkpoints, ckpt_interval
  !
  integer(ik), parameter    :: ckpt_version         = 104_ik    ! Version of the checkpoint file format. Please increment each
                                                                ! time modifications are made to the checkpoint file format!
                                                                ! Versions:
                                                                !  101 = initisl version
                                                                !  102 = added support for adaptive angular momentum
                                                                !  103 = added support for adaptive radial grid
                                                                !  104 = full support for adaptive radial grid
  integer, parameter        :: iu_temp              = 32        ! Some unique unit number to use here.
  character(len=clen), save :: ckpt_save_basename   = ' '       ! Base name for checkpoints. Timestep index and extension '.ckpt'
                                                                ! will be appended. Blank disables checkpoints.
  character(len=clen), save :: ckpt_load_filename   = ' '       ! Name of a checkpoint to restart from. Blank means no restart.
  integer(ik), save         :: ckpt_max_checkpoints = 3_ik      ! Maximum number of checkpoints to keep.
  integer(ik), save         :: ckpt_interval        = 10000_ik  ! Number of timesteps between checkpoints.
  !
  character(len=clen), save :: rcsid_checkpoint_tools = "$Id: checkpoint_tools.f90,v 1.17 2022/10/08 17:24:26 ps Exp $"
  !
  type ckpt_data
    private                                                         ! This data is not for exterrrnal consumption
    logical                          :: initialized      = .false.  ! .true. if we've been here before
    integer(ik)                      :: checkpoint_count            ! Number of checkpoints in the checkpoint_list()
    character(len=clen), allocatable :: checkpoint_list(:)          ! List of files containing checkpoint data. We'll
                                                                    ! allocate ckpt_max_checkpoints+1 entries, which 
                                                                    ! simplifies the cleanup logic a bit.
    logical                          :: do_restart                  ! .true. if we are restarting
    logical                          :: do_checkpoint               ! .true. if we are checkpointing
    !
    !  Entries below are only relevant if do_restart==.true.
    !
    integer(ik)                      :: its                         ! Time step at which we restart
    real(xk)                         :: t_vpot(0:3)                 ! Time and vector-potential at the restart time
                                                                    ! (consistency check)
  end type ckpt_data
  !
  contains
  !
  !  This routine is called once per time step. Both checkpoints and restarts are its responsibility.
  !
  subroutine ckpt_do_checkpoint(ck,checkpoint_go,its,t_vpot,units,wfn_l,wfn_r,tsurf)
    type(ckpt_data), intent(inout) :: ck            ! Checkpointer state. All global state goes here
    logical, intent(out)           :: checkpoint_go ! .false. - We are waiting to reach time step specified in the checkpoint file
                                                    !           Skip all processing steps for this timestep
                                                    ! .true.  - Perform normal propagation and property accumulation for this time step
    integer(ik), intent(in)        :: its           ! Current time step
    real(xk), intent(in)           :: t_vpot(0:3)   ! Current values of time and vector-potential (spherical coordinates)
    integer(ik), intent(in)        :: units(:)      ! Curently open output units which need to be flushed.
                                                    ! Negtative unit numbers will be ignored.
    type(sd_wfn), intent(inout)    :: wfn_l         ! Left wavefunction
    type(sd_wfn), intent(inout)    :: wfn_r         ! Right wavefunction
    type(sts_data), intent(inout)  :: tsurf         ! tSURFF/iSURF state
    !
    !  Early return if checkpoints/restarts are not desired.
    !
    if (ckpt_save_basename==' ' .and. ckpt_load_filename==' ') then
      checkpoint_go = .true.
      return
    end if
    ! 
    call TimerStart('Checkpoint')
    if (.not.ck%initialized) then
      !
      !  This is the first time we are called. Initialize state and optionally load checkpointed data
      !
      call ckpt_initialize(ck)
      if (ckpt_load_filename/=' ') then
        ck%do_restart = .true.
        if (nts%this_node==1) then
          call ckpt_load_checkpoint(ck,wfn_l,wfn_r,tsurf)
        end if
        !
        !  Checkpoint reload may have changed sd_nradial on the master node. Do some extra
        !  initialization to compensate.
        !
        if (sd_adaptive_r) then
          call nt_broadcast(sd_nradial)
          call wt_resize(wfn_l,sd_nradial)  ! Will do nothing if the size didn't change
          call wt_resize(wfn_r,sd_nradial)
          call sd_initialize(repeat=.true.)
        end if
        !
        call nt_broadcast(ck%its)
        call nt_broadcast(ck%t_vpot)
        call nt_broadcast(wfn_l)
        call nt_broadcast(wfn_r)
        call nt_broadcast(tsurf)
        call nt_rebalance(max(wfn_l%lmax,wfn_r%lmax))
      end if
      if (ckpt_save_basename/=' ') then
        ck%do_checkpoint = .true.
      end if
    end if
    !
    if (ck%do_restart) then
      !
      !  We are trying to restart. Have we reached the expected time step yet?      
      !
      if (its==ck%its) then
        if (any(abs(t_vpot-ck%t_vpot)>=10._rk*spacing(max(abs(t_vpot),abs(ck%t_vpot))))) then
          !
          !  Oops. The time step is at the correct point, but the time and/or vector-potential are not
          !
          write (out,"(/'Data in the checkpoit file does not match the simulation at time step ',i0)") its
          write (out,"('     t = ',g42.32e3,' vs ',g42.32e3)") t_vpot(0), ck%t_vpot(0)
          write (out,"('     a = ',g42.32e3,' vs ',g42.32e3)") t_vpot(1), ck%t_vpot(1)
          write (out,"(' theta = ',g42.32e3,' vs ',g42.32e3)") t_vpot(2), ck%t_vpot(2)
          write (out,"('   phi = ',g42.32e3,' vs ',g42.32e3)") t_vpot(3), ck%t_vpot(3)
          call flush_wrapper(out)
          stop 'checkpoint_tools%ckpt_do_checkpoint - checkpoint mismatch'
        end if
        !
        !  Everything matches; resume normal calculations
        !
        checkpoint_go = .true.
        ck%do_restart = .false.
      else ! its/=ck%its
        checkpoint_go = .false.
      end if
    else ! .not. ck%do_restart
      !
      !  We are in normal operation. Let's see whether a checkpoint is due.
      !
      checkpoint_go = .true.
      if (ck%do_checkpoint) then
        if (modulo(its+1,ckpt_interval)==0) then
          !
          !  Checkpoint is needed. Do not forget to synchronize across nodes!
          !
          call nt_merge_all(wfn_l)
          call nt_merge_all(wfn_r)
          call nt_merge_all(tsurf)
          call ckpt_flush_all(units)
          if (nts%this_node==1) then
            call ckpt_save_checkpoint(ck,its,t_vpot,wfn_l,wfn_r,tsurf)
            call ckpt_cleanup_checkpoints(ck)
          end if
        end if
      end if
    end if
    !
    call TimerStop('Checkpoint')
  end subroutine ckpt_do_checkpoint
  !
  !
  !  Internal routines below this point
  ! 
  subroutine ckpt_initialize(ck)
    type(ckpt_data), intent(inout) :: ck   ! Checkpointer state
    !
    integer(ik) :: alloc
    !
    allocate (ck%checkpoint_list(ckpt_max_checkpoints+1),stat=alloc)
    if (alloc/=0) then
      write (out,"('checkpoint_tools%ckpt_initialize - allocation failed with code ',i0)") alloc
      stop 'checkpoint_tools%ckpt_initialize - allocation failed'
    end if
    ck%initialized      = .true.
    ck%checkpoint_count = 0
    ck%do_checkpoint    = .false.
    ck%do_restart       = .false.
  end subroutine ckpt_initialize
  !
  subroutine ckpt_load_checkpoint(ck,wfn_l,wfn_r,tsurf)
    type(ckpt_data), intent(inout) :: ck            ! Checkpointer state. All global state goes here
    type(sd_wfn), intent(inout)    :: wfn_l         ! Left wavefunction
    type(sd_wfn), intent(inout)    :: wfn_r         ! Right wavefunction
    type(sts_data), intent(inout)  :: tsurf         ! tSURFF/iSURF state
    !
    character(len=clen) :: action
    integer(ik)         :: ios
    logical             :: options (2) ! Optionally present data
    character(len=9)    :: tag
    integer(ik)         :: version
    logical             :: options2(2) ! Copy of options from the checkpoint file
    integer(ik)         :: nr          ! Currently-active sd_nradial from checkpoint. May differ from
                                       ! the currently-active grid in memory
    !
    !  Which optional items do we have?
    !
    options(1) = allocated(tsurf%vphase)  ! Volkov block
    options(2) = allocated(tsurf%fcgr)    ! Coulomb block
    !
    !  All errors in checkpoint loading are fatal, but must be reported
    !
    error_block: do
      action = 'opening'
      open (iu_temp,file=trim(ckpt_load_filename),form='unformatted',action='read',status='old',iostat=ios)
      if (ios/=0) exit error_block
      action = 'reading header'
      read (iu_temp,iostat=ios) tag, version, options2, ck%its, ck%t_vpot, nr
      if (ios/=0) exit error_block
      !
      !  Header verification can get a little hairy: most Fortran compilers represent .false.
      !  as zero; however, the value use for .true. differs compiler-by-compiler. This can
      !  lead to puzzling output below.
      !
      action = 'verifying header'
      if (tag/='SCID CKPT' .or. version/=ckpt_version .or. any(options.neqv.options2)) then
        write (out,"('Checkpoint header mismatch')")
        ! gfortran introduces an array temporary for the transfer() below. There
        ! does not seem to be anything we could do about it: adding an explicit
        ! temporary still creates a temporary array ...
        write (out,"('Expected: ',a,' V.',i0,' options = ',2l1,'(',2z8,')')") &
               'SCID CKPT', ckpt_version, options, transfer(options,1_ik,2)
        write (out,"('Received: ',a,' V.',i0,' options = ',2l1,'(',2z8,')')") &
               tag, version, options2, transfer(options2,1_ik,2)
        exit error_block
      end if
      !
      !  It looks like we are restoring an adaptive checkpoint; resize wavefunction grid.
      !  We may need to re-initialize data arrays later as well
      !
      if (nr/=sd_nradial) then
        write (out,"('WARNING: Checkpoint is saved for the ',i0,'-point radial grid.')") nr
        call wt_resize(wfn_l,nr)
        call wt_resize(wfn_r,nr)
      end if
      !
      action = 'reading left wavefunction'
      read (iu_temp,iostat=ios) wfn_l%lmax, wfn_l%nradial, wfn_l%wfn
      if (ios/=0) exit error_block
      action = 'reading right wavefunction'
      read (iu_temp,iostat=ios) wfn_r%lmax, wfn_r%nradial, wfn_r%wfn
      if (ios/=0) exit error_block
      if (options(1)) then
        action = 'reading Volkov tSURF/iSURF data'
        read (iu_temp,iostat=ios) tsurf%vphase, tsurf%pref, tsurf%amplitude
        if (ios/=0) exit error_block
      end if
      if (options(2)) then
        action = 'reading Coulomb iSURF data'
        read (iu_temp,iostat=ios) tsurf%fcgr, tsurf%fcamp, tsurf%coulamp
        if (ios/=0) exit error_block
      end if
      action = 'closing'
      close (iu_temp,status='keep',iostat=ios)
      if (ios/=0) exit error_block
      !
      if (.not.sd_adaptive_l .and.  (wfn_l%lmax<sd_lmax .or. wfn_r%lmax<sd_lmax)) then
        wfn_l%lmax = sd_lmax
        wfn_r%lmax = sd_lmax
        write (out,"('WARNING: Restating adaptive-L checkpoint. Adaptive Lmax reset to sd_lmax')")
      end if
      wfn_l%lmax_top = wfn_l%lmax
      wfn_r%lmax_top = wfn_r%lmax
      if (sd_adaptive_r) then
        sd_nradial = nr
        call sd_initialize(repeat=.true.)
      end if
      if (.not.sd_adaptive_r .and. (wfn_l%nradial<sd_nradial_max .or. wfn_r%nradial<sd_nradial_max)) then
        wfn_l%nradial = sd_nradial_max
        wfn_r%nradial = sd_nradial_max
        call wt_resize(wfn_l,sd_nradial_max)
        call wt_resize(wfn_r,sd_nradial_max)
        write (out,"('WARNING: Restating adaptive-R checkpoint. Adaptive Nradial reset to sd_nradial')")
      end if
      !
      write (out,"('State at the beginning of time step ',i0,' restored from ',a)") ck%its, trim(ckpt_load_filename)
      return
    end do error_block
    write (out,"('FATAL: Encountered error ',i0,' while ',a,' checkpoint file ',a)") ios, trim(action), trim(ckpt_load_filename)
    call flush_wrapper(out)
    stop 'checkpoint_tools%ckpt_load_checkpoint'
  end subroutine ckpt_load_checkpoint
  !
  subroutine ckpt_flush_all(units)
    integer(ik), intent(in)        :: units(:)      ! Curently open output units which need to be flushed.
                                                    ! Negative unit numbers will be ignored.
    !
    integer(ik) :: i
    logical     :: isopen
    !
    scan_units: do i=1,size(units)
      if (units(i)<0) cycle scan_units
      inquire (unit=units(i),opened=isopen)
      if (isopen) call flush_wrapper(units(i))
    end do scan_units
  end subroutine ckpt_flush_all
  !
  subroutine ckpt_save_checkpoint(ck,its,t_vpot,wfn_l,wfn_r,tsurf)
    type(ckpt_data), intent(inout) :: ck            ! Checkpointer state. All global state goes here
    integer(ik), intent(in)        :: its           ! Current time step
    real(xk), intent(in)           :: t_vpot(0:3)   ! Current values of time and vector-potential (spherical coordinates)
    type(sd_wfn), intent(inout)    :: wfn_l         ! Left wavefunction
    type(sd_wfn), intent(inout)    :: wfn_r         ! Right wavefunction
    type(sts_data), intent(inout)  :: tsurf         ! tSURFF/iSURF state
    !
    character(len=clen) :: fname
    character(len=clen) :: action
    integer(ik)         :: ios
    logical             :: options(2) ! Optionally present data
    !
    write (fname,"(a,'-',i0,'.ckpt')") trim(ckpt_save_basename), its
    !
    !  Which optional items do we have?
    !
    options(1) = allocated(tsurf%vphase)  ! Volkov block
    options(2) = allocated(tsurf%fcgr)    ! Coulomb block
    !
    !  All I/O errors in checkpoint creation are non-fatal (they may be transients), but must be reported
    !
    error_block: do
      action = 'creating'
      open (iu_temp,file=trim(fname),form='unformatted',action='write',status='new',iostat=ios)
      if (ios/=0) exit error_block
      action = 'writing header'
      write (iu_temp,iostat=ios) 'SCID CKPT', ckpt_version, options, its, t_vpot, sd_nradial
      if (ios/=0) exit error_block
      action = 'writing left wavefunction'
      write (iu_temp,iostat=ios) wfn_l%lmax, wfn_l%nradial, wfn_l%wfn
      if (ios/=0) exit error_block
      action = 'writing right wavefunction'
      write (iu_temp,iostat=ios) wfn_r%lmax, wfn_r%nradial, wfn_r%wfn
      if (ios/=0) exit error_block
      if (options(1)) then
        action = 'writing Volkov tSURF/iSURF data'
        write (iu_temp,iostat=ios) tsurf%vphase, tsurf%pref, tsurf%amplitude
        if (ios/=0) exit error_block
      end if
      if (options(2)) then
        action = 'writing Coulomb iSURF data'
        write (iu_temp,iostat=ios) tsurf%fcgr, tsurf%fcamp, tsurf%coulamp
        if (ios/=0) exit error_block
      end if
      action = 'closing'
      close (iu_temp,status='keep',iostat=ios)
      if (ios/=0) exit error_block
      !
      write (out,"('State at the beginning of time step ',i0,' saved to ',a)") its, trim(fname)
      ck%checkpoint_count = ck%checkpoint_count + 1
      if (ck%checkpoint_count > size(ck%checkpoint_list)) stop 'checkpoint_tools%ckpt_save_checkpoint - checkpoint_list overflow'
      ck%checkpoint_list(ck%checkpoint_count) = fname
      return
    end do error_block
    if (ios/=0) then
      write (out,"('WARNING: Encountered error ',i0,' while ',a,' checkpoint file ',a)") ios, trim(action), trim(fname)
      close (iu_temp,status='delete',iostat=ios) ! Ignore errors at this point
    end if
  end subroutine ckpt_save_checkpoint
  !
  subroutine ckpt_cleanup_checkpoints(ck)
    type(ckpt_data), intent(inout) :: ck            ! Checkpointer state. All global state goes here
    !
    integer(ik) :: ios
    !
    if (ck%checkpoint_count<=ckpt_max_checkpoints) return
    !
    open (iu_temp,file=trim(ck%checkpoint_list(1)),status='old',iostat=ios)
    if (ios==0) close (iu_temp,status='delete',iostat=ios)
    if (ios==0) then
      write (out,"('Deleted old checkpoint file ',a)") trim(ck%checkpoint_list(1))
    else
      write (out,"('WARNING: Error ',i0,' deleting old checkpoint file ',a)") ios, trim(ck%checkpoint_list(1))
    end if
    ck%checkpoint_list(1:ckpt_max_checkpoints) = ck%checkpoint_list(2:ckpt_max_checkpoints+1)
    ck%checkpoint_count = ck%checkpoint_count - 1
  end subroutine ckpt_cleanup_checkpoints
end module checkpoint_tools
