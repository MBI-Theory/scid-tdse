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
!  Analysis aid: projects wavefunction on the field-free solutions,
!  and reports the results. Long-time corrections to the t-SURF 
!  PES and scattering-wave decomposition are also calculated in
!  ca_analyze()
!
module composition_analysis
  use accuracy
  use constants
  use timer
  use spherical_data
  use spherical_tsurf
  use spherical_tsurf_data
  use wavefunction_tools
  use node_tools
!$ use OMP_LIB
  implicit none
  private
  public ca_analyze
  public ca_maxram
  public rcsid_composition_analysis
  !
  character(len=clen) :: rcsid_composition_analysis = "$Id: composition_analysis.f90,v 1.20 2025/07/11 15:08:35 ps Exp $"
  !
  real(rk), save :: ca_maxram           = 0._rk       ! Maximum amount of memory which can be used during the analysis step
                                                      ! This limit does NOT include the memory needed to compute atomic 
                                                      ! solutions on the fly.
                                                      ! 0 means no limit.
  !
  contains
  !
  !  Calculation of the field-free eigenspectrum
  !
  subroutine ca_analyze(verbose,threshold,emax,wfn_l,wfn_r,tsurf)
    integer(ik), intent(in)       :: verbose   ! Debugging level
    real(rk), intent(in)          :: threshold ! Reporting threshold for the amplitudes of the field-free solutions
    real(rk), intent(in)          :: emax      ! Only include eigenvalues with real part not exceeding emax
    type(sd_wfn), intent(in)      :: wfn_l     ! Left/right wavefunction pair to analyze
    type(sd_wfn), intent(in)      :: wfn_r     ! 
    type(sts_data), intent(inout) :: tsurf     ! A side-effect of the analysis allows calculation of photoelectron
                                               ! spectra; this is done here.
    !
    integer(ik)              :: mval, lval, ispin, iev, mmin, mmax
    integer(ik)              :: alloc
    complex(rk), allocatable :: block_evec(:,:,:)         ! Eigenvectors (unused here)
    complex(rk)              :: block_eval(sd_nradial)    ! Eigenvalues
    complex(rk), allocatable :: amp(:,:,:,:,:), en(:,:)   ! Energies and amplitudes
    complex(rk)              :: wgt
    complex(rk)              :: pop_total, pop_bound, pop_all_bound, pop_all
    complex(rk)              :: lm_norms(sd_mmin:sd_mmax,0:sd_lmax)
    real(rk)                 :: ram_global, ram_thread    ! Memory requirements, in Mbytes
    real(rk)                 :: ecut                      ! Actual value of the energy cut-off, could be
                                                          ! emax or huge(1._rk)
    integer(ik)              :: iemax                     ! Index of the top eigenvalue included in the analysis.
    integer(ik)              :: max_threads
    !
    call TimerStart('Field-free analysis')
    !
    !  Trap for incompatible input parameters: spectral-form tsurf can only be done if 
    !  the full set of eigenstate amplitudes is calculated.
    !
    ecut = emax
    if (sts_active_atend .and. ecut<huge(1._rk)) then
      ecut = huge(1._rk)
      write (out,"(/'WARNING: Spectral-form iSURF requires full eigendecomposition.')") 
      write (out,"( 'WARNING: The value of composition_max_energy has been reset.'/)")
      call flush_wrapper(out)
    end if
    !
    !  Begin by calculating norms of each L,M channel contribution; this is useful check
    !  on the quality of the analysis later.
    !
    lm_norms = 0
    !
    !$omp parallel do default(none) &
    !$omp& private(lval,mmin,mmax,mval,ispin) &
    !$omp& shared(sd_lmax,sd_mmin,sd_mmax,sd_nspin,lm_norms,wfn_l,wfn_r,nts)
    lm_norm_l: do lval=0,sd_lmax
      if (nts%dist(lval)>0) cycle lm_norm_l ! Hook for distributed-memory parallelization
      mmin = max(sd_mmin,-lval)
      mmax = min(lval,sd_mmax)
      lm_norm_m: do mval=mmin,mmax
        lm_norm_s: do ispin=1,sd_nspin
          lm_norms(mval,lval) = lm_norms(mval,lval) + sum(wfn_l%wfn(:,ispin,lval,mval)*wfn_r%wfn(:,ispin,lval,mval))
        end do lm_norm_s
      end do lm_norm_m
    end do lm_norm_l
    !$omp end parallel do
    call nt_add(lm_norms)
    !
    write (out,"(/t5,'Final norm, by total angular momentum'/)")
    write (out,"(1x,a5,1x,a6,2x,a34,2x,a34)") ' L ', '  M ', '      Re(norm)    ', '     Im(norm)     ', &
                                              '---', ' ---', '------------------', '------------------'
    lm_print_l: do lval=0,sd_lmax
      mmin = max(sd_mmin,-lval)
      mmax = min(lval,sd_mmax)
      lm_print_m: do mval=mmin,mmax
        write (out,"(1x,i5,1x,i6,2x,g34.23e3,2x,g34.23e3)") lval, mval, lm_norms(mval,lval)
      end do lm_print_m
    end do lm_print_l
    write (out,"()")
    call flush_wrapper(out)
    !
    !  The rest is expensive; skip it if threshold is negative
    !  We always need the decomposition if we are going to do t-SURF infinite-time terms
    !
    if (threshold<0 .and. .not.sts_active_atend) then
      call TimerStop('Field-free analysis')
      return
    end if
    write (out,"(/t5,'Analyzing the final wavefunction in terms of field-free eigenstates'/)")
    !
    call sts_atend_prepare(tsurf)
    !
    ram_global = ((2._rk)*rk_bytes()/1024._rk**2) * sd_nradial*(sd_lmax+1)*(1+2*sd_nspin*(sd_mmax-sd_mmin+1))
    ram_thread = ((2._rk)*rk_bytes()/1024._rk**2) * 2*real(sd_nradial,kind=rk)**2
    !
    write (out,"(/'    Global memory requirements for the analysis step: ',f20.3,' Mbytes')") ram_global
    write (out,"( 'Per-thread memory requirements for the analysis step: ',f20.3,' Mbytes'/)") ram_thread
    write (out,"('This estimate assumes that the atomic solutions are precomputed.')")
    write (out,"('Computing atomic solutions on the fly will approximately triple per-thread memory requirements'/)")
    !
    max_threads = 0
!$  max_threads = omp_get_max_threads()
    if (ca_maxram>0) then
      max_threads = min(int(floor((ca_maxram-ram_global)/ram_thread),kind=ik),max_threads)
      !
      write (out,"('Maximum amount of RAM to be used in composition analysis: ',f20.3,' Mbytes')") ca_maxram
      write (out,"('                        Maximum number of threads to use: ',i5/)") max_threads
      if (max_threads<1) then
        write (out,"('Not enough memory. Abort.')") 
        call flush_wrapper(out)
        stop 'composition_analysis%ca_analyze - need more memory to proceed'
      end if
    end if
    !
    allocate (amp(sd_nradial,2,sd_nspin,sd_mmin:sd_mmax,0:sd_lmax),en(sd_nradial,0:sd_lmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('ca_analyze: No memory for analysis. Error code ',i0)") alloc
      stop 'composition_analysis%ca_analyze - no memory for analysis'
    end if
    !
    !  The block below, which is the expensive part of the calculation, is
    !  parallelized over both distributed and shared memory. The rest of the
    !  routine is shared-memory parallel, and is repeated on each node.
    !
    amp = 0
    en  = 0
    !$omp parallel default(none) num_threads(max_threads) &
    !$omp& private(lval,mval,mmin,mmax,ispin,block_evec,block_eval,alloc,iemax) &
    !$omp& shared(sd_nradial,sd_lmax,sd_mmin,sd_mmax,sd_nspin,verbose,emax) &
    !$omp& shared(wfn_l,wfn_r,amp,en,tsurf,nts)
    allocate (block_evec(sd_nradial,sd_nradial,2),stat=alloc)
    if (alloc/=0) then
      write (out,"('ca_analyze: No per-thread memory for analysis. Error code ',i0)") alloc
      stop 'composition_analysis%ca_analyze - no per-thread memory for analysis'
    end if
    !$omp do schedule(guided,1)
    scan_l_channels: do lval=0,sd_lmax
      if (nts%dist(lval)>0) cycle scan_l_channels ! Hook for distributed-memory parallelization
      mmin = max(sd_mmin,-lval)
      mmax = min(lval,sd_mmax)
      if (mmin>mmax) cycle scan_l_channels ! This L channel is not present; skip
      !
      call wt_atomic_solutions(verbose,lval,block_eval,block_evec)
      en(:,lval) = block_eval(:)
      iemax = find_iemax(lval)
      scan_spin: do ispin=1,sd_nspin
        ! gfortran is having trouble generating good code for this array assignment
        ! there does not seem to be anything we could do about it though ...
        amp(:iemax,1,ispin,mmin:mmax,lval) = matmul(transpose(block_evec(:,:iemax,1)),wfn_r%wfn(:,ispin,lval,mmin:mmax))
        amp(:iemax,2,ispin,mmin:mmax,lval) = matmul(transpose(block_evec(:,:iemax,2)),wfn_l%wfn(:,ispin,lval,mmin:mmax))
      end do scan_spin
      !
      call sts_atend_terms(tsurf,lval,en(:,lval),block_evec(:,:,2),amp(:,1,1,sd_mmin:sd_mmax,lval))
    end do scan_l_channels
    !$omp end do
    deallocate (block_evec)
    !$omp end parallel
    !
    call nt_add(en)
    !
    !  Our MPI reduction routine is extremely memory-inefficient. Since amp() array can get
    !  quite large, we'll trade off a bit of performance for the decreased memory usage.
    !
    merge_amp_l: do lval=0,sd_lmax
      mmin = max(sd_mmin,-lval)
      mmax = min(lval,sd_mmax)
      merge_amp_m: do mval=mmin,mmax
        merge_amp_spin: do ispin=1,sd_nspin
          call nt_add(amp(:,:,ispin,mval,lval))
        end do merge_amp_spin
      end do merge_amp_m
    end do merge_amp_l
    !
    !  Reporting part, in increasingly excruciating detail
    !
    if (ecut<huge(1._rk)) then
      write (out,"(/'WARNING: Populations listed below include only eigenstates with Re[E]<=',g0.8/)") ecut
    end if
    write (out,"(/t5,'Final populations, by total angular momentum'/)")
    write (out,"((1x,a5,3(2x,a24,1x,a24)))") &
           ' L ', ' Total population ', ' ', ' Bound population ', ' ', ' Continuum population ', ' ', &
           '---', '------------------', ' ', '------------------', ' ', '----------------------', ' '
    pop_all       = 0
    pop_all_bound = 0
    l_resolved_l_channels: do lval=0,sd_lmax
      pop_total = 0
      pop_bound = 0
      iemax     = find_iemax(lval)
      !$omp parallel do default(none) &
      !$omp& private(mval,iev,wgt) reduction(+:pop_total,pop_bound) &
      !$omp& shared(sd_mmin,sd_mmax,sd_nradial,lval,amp,en,iemax)
      l_resolved_m_channels: do mval=max(sd_mmin,-lval),min(lval,sd_mmax)
        l_resolved_evs: do iev=1,iemax
          ! wgt = abs(sum(amp(iev,1,:,mval,lval)*amp(iev,2,:,mval,lval)))
          wgt = sum(amp(iev,1,:,mval,lval)*amp(iev,2,:,mval,lval))
          pop_total = pop_total + wgt
          if (real(en(iev,lval),kind=rk)<0) pop_bound = pop_bound + wgt
        end do l_resolved_evs
      end do l_resolved_m_channels
      !$omp end parallel do
      write (out,"((1x,i5,3(2x,g24.13e3,1x,g24.13e3)))") lval, pop_total, pop_bound, pop_total-pop_bound
      pop_all       = pop_all + pop_total
      pop_all_bound = pop_all_bound + pop_bound
    end do l_resolved_l_channels
    write (out,"()")
    !
    write (out,"()")
    if (ecut<huge(1._rk)) &
    write (out,"('                      Energy cut-off: ',g24.13)") ecut
    write (out,"('Total population across all channels: ',g24.13,1x,g24.13)") pop_all
    write (out,"('                        Bound states: ',g24.13,1x,g24.13)") pop_all_bound
    write (out,"('                    Continuum states: ',g24.13,1x,g24.13)") pop_all-pop_all_bound
    write (out,"()")
    !
    write (out,"(/t5,'Final populations, by total angular momentum and angular momentum projection'/)")
    write (out,"((1x,a5,1x,a6,3(2x,a24,1x,a24)))") &
           ' L ', ' M ', ' Total population ', ' ', ' Bound population ', ' ', ' Continuum population ', ' ', &
           '---', '---', '------------------', ' ', '------------------', ' ', '----------------------', ' '
    lm_resolved_l_channels: do lval=0,sd_lmax
      lm_resolved_m_channels: do mval=max(sd_mmin,-lval),min(lval,sd_mmax)
        pop_total = 0
        pop_bound = 0
        iemax     = find_iemax(lval)
        !$omp parallel do default(none) &
        !$omp& private(iev,wgt) reduction(+:pop_total,pop_bound) &
        !$omp& shared(mval,sd_nradial,lval,amp,en,iemax)
        lm_resolved_evs: do iev=1,iemax
          ! wgt = abs(sum(amp(iev,1,:,mval,lval)*amp(iev,2,:,mval,lval)))
          wgt = sum(amp(iev,1,:,mval,lval)*amp(iev,2,:,mval,lval))
          pop_total = pop_total + wgt
          if (real(en(iev,lval),kind=rk)<0) pop_bound = pop_bound + wgt
        end do lm_resolved_evs
        !$omp end parallel do
        write (out,"((1x,i5,1x,i6,3(2x,g24.13e3,1x,g24.13e3)))") lval, mval, pop_total, pop_bound, pop_total-pop_bound
      end do lm_resolved_m_channels
    end do lm_resolved_l_channels
    write (out,"()")
    !
    write (out,"(/t5,'Large amplitudes of individual field-free states'/)")
    write (out,"(1x,a5,1x,a6,1x,a7,2x,a24,1x,a24,2x,a24,1x,a24,2x,a24,1x,a24,2x,a24,1x,a24)") &
      ' L ', ' M ', ' I ', ' Re[E(i)], H ', ' Im[E(i)] ', '  Re[Wgt]  ', '  Im[Wgt]  ', &
                            ' Re[<I|W>] ', ' Im[<I|W>] ', ' Re[<W|I>] ', ' Im[<W|I>] ', &
      '---', '---', '---', '-------------', '----------', '-----------', '-----------', &
                              '-----------', '-----------', '-----------', '---------'
    state_resolved_l_channels: do lval=0,sd_lmax
      iemax = find_iemax(lval)
      state_resolved_m_channels: do mval=max(sd_mmin,-lval),min(lval,sd_mmax)
        state_resolved_evs: do iev=1,iemax
          if (all(abs(amp(iev,:,:,mval,lval))<threshold)) cycle state_resolved_evs
          ! wgt = abs(sum(amp(iev,1,:,mval,lval)*amp(iev,2,:,mval,lval)))
          wgt = sum(amp(iev,1,:,mval,lval)*amp(iev,2,:,mval,lval))
          write (out,"(1x,i5,1x,i6,1x,i7,2(2x,g24.13e3,1x,g24.13e3),2x,(t126,2(g24.13e3,1x,g24.13e3)))") &
                 lval, mval, iev, en(iev,lval), wgt, amp(iev,:,:,mval,lval)
        end do state_resolved_evs
      end do state_resolved_m_channels
    end do state_resolved_l_channels
    write (out,"()")
    !
    deallocate (amp,en)
    call flush_wrapper(out)
    !
    call TimerStop('Field-free analysis')
    !
    contains

    function find_iemax(lval) result(iemax)
      integer(ik), intent(in) :: lval
      integer(ik)             :: iemax, iev
      !
      iemax = sd_nradial
      if (emax>real(en(sd_nradial,lval),kind=rk)) return
      do_find_iemax: do iev=1,sd_nradial
        if (real(en(iev,lval),kind=rk)>emax) then
          iemax = iev - 1
          exit do_find_iemax
        end if
      end do do_find_iemax
    end function find_iemax
  end subroutine ca_analyze
end module composition_analysis
