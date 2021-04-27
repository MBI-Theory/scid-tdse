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
!  Simple post-processing tool for building simulated photoelectron spectra
!  This is basically the Fortran version of the ./examples/plot-pes.sh
!
!  WARNING: This routine is obsolete. Please use STS_VOLKOV, STS_VOLKOV_ATEND,
!  WARNING: and STS_COULOMB_ATEND options in the main code instead.
!
module build_pes
  use accuracy
  use timer
  implicit none
  private
  public start
  public rcsid_build_pes
  !
  character(len=clen), save :: rcsid_build_pes = "$Id: build_pes.f90,v 1.5 2021/04/26 15:44:44 ps Exp ps $"
  !
  integer, parameter       :: iu_temp               = 25           ! Some arbitrary I/O unit
  !
  integer(ik)              :: verbose               = 2_ik         ! How verbose do we need to be?
  character(len=clen)      :: comment               = ' '          ! Descriptive string, to be copied to the output
  character(len=clen)      :: scid_output           = 'scid.out'   ! Output file from a scid-tdse run
  integer(ik)              :: mmin                  = -256         ! Smallest M value to look for
  integer(ik)              :: mmax                  =  256         ! Largest M value to look for
  integer(ik)              :: lmax                  =  256         ! Largest L value to look for
  integer(ik)              :: max_records           = 10000000     ! Maximum number of states to look for
  real(rk)                 :: e_min                 =  0._rk       ! Minimal energy to include in the spectrum (Hartree)
  real(rk)                 :: e_max                 = 10._rk       ! Maximal energy to include in the spectrum (Hartree)
  integer(ik)              :: e_count               = 10000_ik     ! Number of energy values to include in the spectrum
  character(len=clen)      :: spec_file             = 'spec.table' ! Calculated spectrum
  character(len=clen)      :: peak_file             = 'peak.table' ! Peaks in the spectrum
  !
  integer(ik)              :: n_tot                 = 0            ! Actual number of records processed
  integer(ik)              :: n_records             = 0            ! Actual number of records accepted
  real(rk), allocatable    :: en_tab(:)                            ! Table of eigenstate energies (Hartree)
  real(rk), allocatable    :: gam_tab(:)                           ! Table of eigenstate widths (Hartree)
  real(rk), allocatable    :: wgt_tab(:)                           ! Table of eigenstate amplitudes 
  real(rk), allocatable    :: spec_en(:)                           ! Table of energies in the spectrum (Hartree)
  real(rk), allocatable    :: spec_val(:)                          ! Table of probability density in the spectrum (1/Hartree)
  !
  namelist /pes_input/ verbose, comment, &
                       scid_output, mmin, mmax, lmax, max_records, &
                       e_min, e_max, e_count, &
                       spec_file, peak_file
  !
  contains
  !
  !  This would have been -so- much easier in awk or perl ...
  !
  subroutine read_state_expansion
    integer(ik)         :: ios
    integer(ik)         :: line
    integer(ik)         :: trigger_line
    integer(ik)         :: lval, mval, ival
    real(rk)            :: ree, ime, rew
    character(len=20)   :: action
    character(len=20)   :: state
    character(len=clen) :: buf
    !
    call TimerStart('Read data')
    if (verbose>=0) then
      write (out,"('Reading ',a)") trim(scid_output)
    end if
    line  = 0
    buf   = ' '
    state = 'waiting'
    error_block: do
      action = 'opening'
      open (iu_temp,form='formatted',recl=256,action='read',position='rewind',status='old',file=trim(scid_output),iostat=ios)
      if (ios/=0) exit error_block
      !
      read_loop: do
        action = 'reading'
        line   = line + 1
        read (iu_temp,"(a)",iostat=ios) buf
        if (ios/=0) exit error_block
        select case (state)
          case default
            stop 'build_pes%read_state_expansion - state machine corrupt'
          case ('waiting')
            if (index(buf,'Large amplitudes of individual field-free states')>0) then
              state = 'header'
              trigger_line = line + 3 ! Skip three lines following the header
            end if
          case ('header')
            if (line>=trigger_line) state = 'scanning'
          case ('scanning')
            if (buf==" ") exit read_loop
            action = 'parsing'
            read (buf,"(1x,i3,1x,i4,1x,i5,2(2x,g24.13e3,1x,g24.13e3),2x,(t120,2(g24.13e3,1x,g24.13e3)))",iostat=ios) &
                  lval, mval, ival, ree, ime, rew
            if (ios/=0) exit error_block
            n_tot = n_tot + 1
            if (verbose>=3) then
              write (out,"('Peak ',i0,': l=',i0,' m=',i0,' seq=',i0,' re[e]=',g24.13e3,' im[e]=',g24.13e3,' wgt=',g24.13e3)") &
                     n_tot, lval, mval, ival, ree, ime, rew
            end if
            if (mval>=mmin .and. mval<=mmax .and. lval<=lmax) then
              n_records = n_records + 1
              if (n_records>max_records) stop 'build_pes%read_state_expansion - increase max_records'
              en_tab (n_records) = ree
              gam_tab(n_records) = max(-2*ime,spacing(1._rk))
              wgt_tab(n_records) = rew
            end if
        end select
      end do read_loop
      !
      close (iu_temp,iostat=ios)
      if (ios/=0) exit error_block
      if (verbose>=0) then
        write (out,"('  Total number of state records = ',i0)") n_tot
        write (out,"('States included in the spectrum = ',i0)") n_records
      end if
      call TimerStop('Read data')
      return
    end do error_block
    write (out,"(/'Fatal error ',a,' file ',a)") trim(action), trim(scid_output)
    write (out,"( '   Error code = ',i0)") ios
    write (out,"( ' Current line = ',i0)") line
    write (out,"( 'Current state = ',a)") trim(state)
    write (out,"( '  Line buffer = ',a)") trim(buf)
    stop 'build_pes%read_state_expansion - fatal error'
  end subroutine read_state_expansion
  !
  subroutine build_spectrum
    integer(ik)         :: ipt, ipeak
    real(rk)            :: de
    real(rk)            :: max_de   ! Maximum energy difference where peak still matters
    integer(ik)         :: ip1, ipn ! First and last buckets affected by this peak
    real(rk)            :: alp, bet ! Gaussian parameters
    real(rk)            :: e0       ! Central energy of the Gaussian
    real(rk)            :: e1, en   ! Beginning and end of the affected range
    ! Constants needed in evaluating the pre-factor and the exponent parts of each peak
    real(rk), parameter :: scl_alp = sqrt(4._rk * log(2._rk) / (4._rk*atan2(1._rk,1._rk)))
    real(rk), parameter :: scl_gam =      4._rk * log(2._rk) 
    !
    call TimerStart('Build spectrum')
    !
    de = (e_max-e_min)/(e_count-1)
    fill_energies: do ipt=1,e_count
      spec_en(ipt) = e_min + (ipt-1)*de
    end do fill_energies
    spec_val(:) = 0
    !
    process_peaks: do ipeak=1,n_records
      ! Our peak is given by alp*exp(-bet*(e-e0)**2)
      e0     = en_tab(ipeak)
      alp    = wgt_tab(ipeak) * scl_alp / gam_tab(ipeak)
      bet    = scl_gam / gam_tab(ipeak)**2
      ! Figure out how far this peak will reach on the grid
      max_de = sqrt((log(abs(alp))-log(tiny(alp)))/bet)
      e1  = e0-max_de
      en  = e0+max_de
      if (verbose>=3) then
        write (out,"('Peak: ',i0,' e0=',g23.10e3,' alp=',g23.10e3,' bet=',g23.10e3)") ipeak, e0, alp, bet
        write (out,"('Range: ',g23.10e3,1x,g23.10e3)") e1, en
      end if
      if (e1>e_max .or. en<e_min) cycle process_peaks ! Peak is off-grid
      !
      !  Find all points affected by this peak
      !
      ip1 = 1
      if (e1>e_min) ip1 = max(1+floor  ((e1-e_min)/de),1)
      ipn = e_count
      if (en<e_max) ipn = min(1+ceiling((en-e_min)/de),e_count)
      !
      evaluate_peak: do ipt=ip1,ipn
        spec_val(ipt) = spec_val(ipt) + alp * exp(-bet*(spec_en(ipt)-e0)**2)
      end do evaluate_peak
    end do process_peaks
    call TimerStop('Build spectrum')
  end subroutine build_spectrum
  !
  subroutine print_spectrum
    integer(ik)       :: ios
    integer(ik)       :: ipt
    integer(ik)       :: line
    character(len=20) :: action
    !
    call TimerStart('Print spectrum')
    !
    if (verbose>=0) then
      write (out,"('Saving spectrum to ',a)") trim(spec_file)
    end if
    !
    line = 0
    error_block: do
      action = 'opening'
      open (iu_temp,form='formatted',recl=256,action='write',position='rewind',status='replace',file=trim(spec_file),iostat=ios)
      if (ios/=0) exit error_block
      !
      action = 'writing'
      write_loop: do ipt=1,e_count
        line   = line + 1
        write (iu_temp,"(1x,g28.15e3,1x,g28.15e3)",iostat=ios) spec_en(ipt), spec_val(ipt)
        if (ios/=0) exit error_block
      end do write_loop
      !
      action = 'closing'
      close (iu_temp,iostat=ios)
      if (ios/=0) exit error_block
      call TimerStop('Print spectrum')
      return
    end do error_block
    write (out,"(/'Fatal error ',a,' file ',a)") trim(action), trim(spec_file)
    write (out,"( '   Error code = ',i0)") ios
    write (out,"( ' Current line = ',i0)") line
    stop 'build_pes%print_spectrum - fatal error'
  end subroutine print_spectrum
  !
  subroutine find_peaks
    integer(ik)       :: ios
    character(len=20) :: action
    integer(ik)       :: ipt
    integer(ik)       :: npt
    real(rk)          :: prev
    !
    call TimerStart('Find peaks')
    !
    if (verbose>=0) then
      write (out,"('Saving peaks to ',a)") trim(peak_file)
    end if
    !
    npt  = 0
    prev = 0
    error_block: do
      action = 'opening'
      open (iu_temp,form='formatted',recl=256,action='write',position='rewind',status='replace',file=trim(peak_file),iostat=ios)
      if (ios/=0) exit error_block
      !
      action = 'writing'
      find_local_maxima: do ipt=2,e_count-1
        if (spec_val(ipt-1)>=spec_val(ipt)) cycle find_local_maxima
        if (spec_val(ipt+1)>=spec_val(ipt)) cycle find_local_maxima
        !
        npt = npt + 1
        write (iu_temp,"(1x,g28.15e3,1x,g28.15e3,1x,g28.15e3)") spec_en(ipt), spec_val(ipt), spec_en(ipt) - prev
        prev = spec_en(ipt)
      end do find_local_maxima
      !
      action = 'closing'
      close (iu_temp,iostat=ios)
      if (ios/=0) exit error_block
      write (out,"('Found ',i0,' peaks')") npt
      call TimerStop('Find peaks')
      return
    end do error_block
    write (out,"(/'Fatal error ',a,' file ',a)") trim(action), trim(peak_file)
    write (out,"( '   Error code = ',i0)") ios
    write (out,"( ' Current line = ',i0)") npt
    stop 'build_pes%fibd_peaks - fatal error'
  end subroutine find_peaks
  !
  subroutine start
    write (out,"('Version: ',a/)") __BUILD_ID__
    write (out,"('  Integer kind = ',i0,' (',i0,' decimals)')") kind(1_ik), int(log10(huge(1_ik)+1._rk))
    write (out,"('     Real kind = ',i0,' (',i0,' decimals)')") kind(1._rk), precision(1._rk)
    write (out,"('Aux. real kind = ',i0,' (',i0,' decimals)')") kind(1._xk), precision(1._xk)
    write (out,"()")
    !
    call TimerStart('start')
    !
    read (input,nml=pes_input)
    write (out,"(' ===== begin simulation parameters ===== ')")
    write (out,nml=pes_input)
    write (out,"(' ====== end simulation parameters ====== ')")
    !
    write (out,"(/a/)") trim(comment)
    !
    allocate (en_tab(max_records),gam_tab(max_records),wgt_tab(max_records))
    allocate (spec_en(e_count),spec_val(e_count))
    !
    call read_state_expansion
    !
    call build_spectrum
    !
    call print_spectrum
    !
    call find_peaks
    !
    call TimerStop('start')
    call TimerReport
  end subroutine start
end module build_pes
!
program main
  use build_pes

  call start
end program main
