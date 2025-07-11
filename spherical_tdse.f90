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
!
!  Solving atomic TDSE in the presence of a strong laser field.
!  The structure of this code is strongly inspired by the paper:
!
!   H.G. Muller, An Efficient Propagation Scheme for the Time-Dependent 
!   Schroedinger Equation in the Velocity Gauge, Laser Physics, 9, 138-148 (1999)
!
!  However, the capabilities of the code, and some of the numerical 
!  details are quite different from HGM's. They are described in:
!
!   S. Patchkovskii and H.G. Muller, "Simple, accurate, and efficient 
!   implementation of 1-electron atomic time-dependent Schrödinger equation 
!   in spherical coordinates", Comp. Phys. Comm. 199, 153-169 (2016). 
!   doi:10.1016/j.cpc.2015.10.014
!
!  In parallel runs, this code may benefit significantly from large page support,
!  especially on systems with 2K default pages and small TLB (x86_64!). Please
!  refer to your system documentation on how to enable large page support. Most
!  recent Linux systems will come with anonymous huge page support enabled.
!
module spherical_tdse
  use spherical_tdse_data
  use spherical_tdse_field
  use spherical_tdse_initialwf
  use spherical_tdse_io
  use spherical_tdse_propagate
  implicit none
  private
  public start
  public rcsid_spherical_tdse
  !
  character(len=clen), save :: rcsid_spherical_tdse = "$Id: spherical_tdse.f90,v 1.142 2025/07/11 15:08:35 ps Exp $"
  !
  contains
  !
  subroutine version_header
    use ISO_FORTRAN_ENV
    external :: versions
    if (verbose<0) return
    !
    write (out,"(t5,a)") trim(rcsid_spherical_tdse)
    write (out,"(t5,a)") trim(rcsid_spherical_tdse_data)
    write (out,"(t5,a)") trim(rcsid_spherical_tdse_field)
    write (out,"(t5,a)") trim(rcsid_spherical_tdse_initialwf)
    write (out,"(t5,a)") trim(rcsid_spherical_tdse_io)
    call versions  ! Expect gfortran warning
    write (out,"()")
    write (out,"('Compiler: ',a)") compiler_version()
    write (out,"('Build flags: ',a)") compiler_options()
    write (out,"()")
    write (out,"('    Integer kind = ',i0,' (',i0,' decimals)')") kind(1_ik), int(log10(huge(1_ik)+1._rk))
    write (out,"('       Real kind = ',i0,' (',i0,' decimals)')") kind(1._rk), precision(1._rk)
    write (out,"('  Aux. real kind = ',i0,' (',i0,' decimals)')") kind(1._xk), precision(1._xk)
    write (out,"('Max. LAPACK kind = ',i0,' (',i0,' decimals)')") kind(1._lrk), precision(1._lrk)
    write (out,"()")
  end subroutine version_header
  !
  subroutine math_init
    real(rk) :: dummy
    !
    if (log10(huge(1._rk))>=4930._rk) then
      ! Quad precision; don't forget "factorial_slack"
      dummy = MathFactorial(1750_ik-5_ik)
      dummy = MathLogFactorial(20000_ik)
    else if (log10(huge(1._rk))>=308._rk) then
      ! Double precision
      dummy = MathFactorial(170_ik-5_ik)
      dummy = MathLogFactorial(10000_ik)
    else
      ! Single precision
      dummy = MathFactorial(34_ik-5_ik)
      dummy = MathLogFactorial(500_ik)
    end if
  end subroutine math_init
  !
  subroutine prepare_initial_wavefunctions
    integer(ik) :: alloc, ie, ios
    logical     :: scoreboard(ensemble_size)
    !
    allocate (wfns_l(ensemble_size),wfns_r(ensemble_size),ckpts(ensemble_size), &
              tsurfs(ensemble_size),iu_detail(ensemble_size), ens_visualize_1stseq(ensemble_size), &
              ens_detail_output(ensemble_size), ens_final_wf_dump_prefix(ensemble_size), &
              ens_initial_wf_dump_prefix(ensemble_size), ens_visualize_prefix(ensemble_size), &
              ens_ckpt_save_basename(ensemble_size), ens_ckpt_load_filename(ensemble_size), &
              ens_volkov_opendx(ensemble_size), ens_volkov_table(ensemble_size), &
              ens_coulomb_waves(ensemble_size), ens_coulomb_table(ensemble_size), &
              ens_coulomb_opendx(ensemble_size), stat=alloc)
    if (alloc/=0) then
      write (out,"('prepare_initial_wavefunctions: Error ',i0,' allocating memory')") alloc
      stop 'spherical_tdse%prepare_initial_wavefunctions - No memory for descriptors'
    end if
    !
    iu_detail(:) = -1_ik
    scoreboard(:) = .false.
    call mark_and_call_init
    prepare_ensemble: do ie=2,ensemble_size
      read (input,nml=sph_multi,iostat=ios)
      write (out,"()")
      write (out,"(' ==== Ensemble input ',i0,' (',i0,'-th &sph_multi namelist) ====')") ie, ie-1
      write (out,nml=sph_multi)
      write (out,"(' ==== End of ensemble input ',i0,' ====')") ie-1
      if (ios/=0) then
        write (out,"('ERROR: Reading &sph_multi namelist #',i0,' failed with error code ',i0)") ie-1, ios
        call flush_wrapper(out)
        stop 'spherical_tdse%prepare_initial_wavefunctions - Bad &sph_multi'
      end if
      call mark_and_call_init
    end do prepare_ensemble
    if (.not.all(scoreboard)) then
      write (out,"('ERROR: Some wavefunctions are not initialized')")
      write (out,"('scoreboard: ',(t14,(25(l1,1x))))") scoreboard
      stop 'spherical_tdse%prepare_initial_wavefunctions - Missing ensemble index'
    end if
    !
    contains 
    subroutine mark_and_call_init
      if (ensemble_index<1 .or. ensemble_index>ensemble_size) then
        write (out,"('ERROR: ensemble_index=',i0,' is not between 1 and ',i0)") &
               ensemble_index, ensemble_size
        stop 'spherical_tdse%prepare_initial_wavefunctions - Bad ensemble index'
      end if
      if (scoreboard(ensemble_index)) then
        write (out,"('ERROR: ensemble_index=',i0,' appears more than once in the input')") &
               ensemble_index
        stop 'spherical_tdse%prepare_initial_wavefunctions - Duplicate ensemble index'
      end if
      scoreboard(ensemble_index) = .true.
      ens_detail_output         (ensemble_index) = detail_output
      ens_final_wf_dump_prefix  (ensemble_index) = final_wf_dump_prefix
      ens_initial_wf_dump_prefix(ensemble_index) = initial_wf_dump_prefix
      ens_visualize_prefix      (ensemble_index) = visualize_prefix
      ens_ckpt_save_basename    (ensemble_index) = ckpt_save_basename
      ens_ckpt_load_filename    (ensemble_index) = ckpt_load_filename
      ens_volkov_opendx         (ensemble_index) = sts_volkov_opendx
      ens_volkov_table          (ensemble_index) = sts_volkov_table
      ens_coulomb_waves         (ensemble_index) = sts_coulomb_waves
      ens_coulomb_table         (ensemble_index) = sts_coulomb_table
      ens_coulomb_opendx        (ensemble_index) = sts_coulomb_opendx
      ens_visualize_1stseq      (ensemble_index) = visualize_1stseq
      ! function from spherical_tdse_initialwf
      call prepare_one_initial_wavefunction(wfns_l(ensemble_index),wfns_r(ensemble_index))
      ! initialize checkpoint filenames from global variables.
      call ckpt_initialize(ckpts(ensemble_index))
    end subroutine mark_and_call_init
  end subroutine prepare_initial_wavefunctions
  !
  subroutine surf_init_all_instances
    integer(ik) :: ie
    !
    surf_instances: do ie=1,ensemble_size
      call sts_initialize_instance(tsurfs(ie),wfns_l(ie),wfns_r(ie), &
                                   volkov_opendx =ens_volkov_opendx (ie), &
                                   volkov_table  =ens_volkov_table  (ie), &
                                   coulomb_waves =ens_coulomb_waves (ie), &
                                   coulomb_table =ens_coulomb_table (ie), &
                                   coulomb_opendx=ens_coulomb_opendx(ie))
    end do surf_instances
  end subroutine surf_init_all_instances
  !
  subroutine choose_rotation_code
    select case (rotation_mode)
      case default
        write (out,"('spherical_tdse%choose_rotation_mode: rotation mode ',a,' is not recognized')") trim(rotation_mode)
        stop 'spherical_tdse%choose_rotation_mode - bad rotation_mode'
      case ('none','brute force','sparse')
      case ('auto')
        ! If we use a range of M values, presumably we want rotation
        if (sd_mmin/=sd_mmax) rotation_mode = 'sparse'
        ! If we know that VP is along Z, we do not need rotation
        ! This includes 'z Gaussian', 'z Sin2', etc
        if (vp_shape(1:2)=='z '  ) rotation_mode = 'none'
        if (vp_shape(1:3)=='zz ' ) rotation_mode = 'none'
        if (vp_shape(1:4)=='zzz ') rotation_mode = 'none'
        if (vp_shape(1:4)=='zero') rotation_mode = 'none'
        ! If we still do not know, we choose the general case (rotation)
        if (rotation_mode=='auto') rotation_mode = 'sparse'
        write (out,"(/'Chosen rotation_mode = ''',a,''''/)") trim(rotation_mode)
    end select 
    !
    !  Issue warnings based on the final choice of rotation_mode
    !
    select case (rotation_mode)
      case default
        stop 'spherical_tdse%choose_rotation_mode - logic error'
      case ('none')
        if (sd_mmin/=sd_mmax) then
          write (out,"(/'WARNING: Multiple, uncoupled M values are propagated.'/)") 
        end if
        if (vp_shape(1:2)/='z ' .and. vp_shape(1:3)/='zz ' .and. vp_shape(1:4)/='zzz ' .and. vp_shape(1:4)/='zero') then
          write (out,"(/'WARNING: Vector-potential may contain components along X/Y, but the propagator')")
          write (out,"( 'WARNING: assumes it is along Z. The results are likely incorrect.'/)")
        end if
      case ('brute force','sparse')
        if (sd_mmin>-sd_lmax .or. sd_mmax<sd_lmax) then
          write (out,"(/'WARNING: Not all possible projections of angular momentum are included in propagation.')")
          write (out,"( 'WARNING: The results are likely incorrect.'/)")
        end if
    end select 
    call flush_wrapper(out)
  end subroutine choose_rotation_code
  !
  subroutine report_adaptive
    integer(ik) :: ie
    !
    if (.not.sd_adaptive) return
    ensemble_report: do ie=1,ensemble_size
      if (ensemble_size<=1) then
        write (out,"()") 
      else
        write (out,"(/'For ensemble wavefunction ',i0,':'/)") ie
      end if
      if (sd_adaptive_l) write (out,"(' Largest left/right Lmax ',2i4)") wfns_l(ie)%lmax_top, wfns_r(ie)%lmax_top
      if (sd_adaptive_l) write (out,"('   Final left/right Lmax ',2i4)") wfns_l(ie)%lmax,     wfns_r(ie)%lmax
      if (sd_adaptive_r) write (out,"('Final left/right nradial ',2i8)") wfns_l(ie)%nradial,  wfns_r(ie)%nradial
      write (out,"()")
    end do ensemble_report
  end subroutine report_adaptive
  !
  subroutine reconstruct_all_left
    integer(ik) :: ie
    if (.not.reconstruct_left) return
    !
    call TimerStart('Left reconstruction')
    write (out,"(/'Reconstructing left wavefunction(s)'/)")
    reconstruct: do ie=1,ensemble_size
      call wt_reconstruct_left(verbose,wfns_l(ie),wfns_r(ie))
    end do reconstruct
    call TimerStop('Left reconstruction')
  end subroutine reconstruct_all_left
  !
  subroutine analyze_and_surf_all
    integer(ik) :: ie

    call TimerStart('Final analysis')
    write (out,"(/'Analyzing wavefunction composition'/)")
    call flush_wrapper(out)
    analyze: do ie=1,ensemble_size
      if (ensemble_size>1) then
        write (out,"(/'Ensemble component ',i0/)") ie
      end if
      call ca_analyze(verbose,composition_threshold,composition_max_energy,wfns_l(ie),wfns_r(ie),tsurfs(ie))
      !
      call sts_atend_direct(tsurfs(ie),wfns_r(ie))
      !
      call sts_report(tsurfs(ie),'At end')
    end do analyze
    call TimerStop('Final analysis')
  end subroutine analyze_and_surf_all
  !
! subroutine memory_trace
!   use ISO_C_BINDING
!   interface
!     subroutine mtrace() bind(C,name="mtrace")
!     end subroutine mtrace
!   end interface
!   !
!   !  Environment variable MALLOC_TRACE must also be set, pointing to a writeable filename
!   !
!   call mtrace
! end subroutine memory_trace
  !
  subroutine start
    !$ use OMP_LIB
    logical             :: have_openmp
    integer(ik)         :: ie, ios
    character(len=clen) :: iom
    !
!   call memory_trace  ! this module
    !
    write (out,"('Version: ',a/)") __BUILD_ID__
    !
    call TimerStart('start')
    !
    call TimerStart('Initialization')
    !
    !  Some MPI libraries (OpenMPI - I am looking at you!) do not support redirection of
    !  stdin to all ranks. We will therefore try to get the name of the stdin file from
    !  the command-line arguments.
    !
    call nt_force_stdin
    !
    !  The sequence below is a little awkward. We can't complete initialization of a multi-node
    !  run until we've red in the input parameters.
    !
    read (input,nml=sph_tdse,iostat=ios,iomsg=iom)
    !
    call nt_initialize
    !
    call version_header  ! This module
    !
    write (out,"(' ===== begin simulation parameters ===== ')")
    write (out,nml=sph_tdse)
    write (out,"(' ====== end simulation parameters ====== ')")
    if (ios/=0) then
      write (out,"('ERROR: Reading of the input namelist &sph_tdse failed with error code ',i0)") ios
      write (out,"('ERROR: ',a)") trim(iom)
      call flush_wrapper(out)
      stop 'spherical_tdse%start - Bad input namelist'
    end if
    !
    write (out,"(/a/)") trim(comment)
    !
    have_openmp = .false.
    !$ have_openmp = .true.
    !$ if (omp_num_threads/=0) then
    !$   write (out,"('Forcing number of OpenMP threads to ',i0)") omp_num_threads
    !$   call omp_set_num_threads(omp_num_threads)
    !$ end if
    !$ write (out,"('Maximum number of OpenMP threads for this run is ',i0)") omp_get_max_threads()
    if (omp_num_threads>0 .and. .not.have_openmp) then
      write (out,"(/'WARNING: omp_num_threads input parameter has been specified, but the code was built')")
      write (out,"( 'WARNING: without OpenMP support. Execution will continue on a single CPU core.'/)")
    end if
    if (nts%n_nodes>1) then
      write (out,"('Executing on ',i0,' nodes. This is node ',i0)") nts%n_nodes, nts%this_node
    end if
    !
    if (skip_left_propagation) then
      write (out,"(/'WARNING: Left wavefunction will not be propagated.')")
      write (out,"( 'WARNING: MOST OBSERVABLES WILL BE INCORRECT.')")
      write (out,"( 'WARNING: Don''t do this unless you know exactly what you are doing!'/)")
    end if
    !
    call math_init  ! This module
    !
    !  Set up potential and CAP evaluation
    !
    call pt_initialize
    call cap_initialize
    !
    !  Set up grids and operators. This is the first call, produce more output
    !
    call sd_initialize(repeat=.false.)
    !
    !  Check accuracy of numerical derivatives
    !
    if (.not.skip_tests) call derivatives_test(verbose)
    !
    !  Construct field-free states using direct diagonalization
    !
    if (.not.skip_tests) call fieldfree_test(verbose)
    !
    ! Caches must be initialized before prepare_initial_wavefunction(), 
    ! which may refer to them.
    !
    call wt_init_memory_caches(verbose)
    !
    ! Calculate transition matrix elements between bound states if asked for
    ! Only makes sense if memory caches for atomic solutions is enabled and cap disabled
    !
    call wt_transition_matrix_elements(cap_name)
    !
    !  We need a wavefunction to start from. If running an ensemble simulation,
    !  prepare_intial_wavefunctions will also read the additional input for all
    !  ensemble components.
    !
    if (ensemble_size<=0) then
      write (out,"('ERROR: ensemble_size must be at least 1. It was ',i0)") ensemble_size
      call flush_wrapper(out)
      stop 'spherical_tdse%start - bad ensemble_size'
    end if
    call prepare_initial_wavefunctions ! This module
    !
    !  Prepare for a t-SURF calculation
    !
    call sts_initialize_global(dt,sqrt(vp_scale**2+vp_scale_x**2))
    !  Warning: visualize_wavefunctions() depends on sts_initialize_instance() calculating
    !           the correct value of tsurfs(1)%wfn_scale, even though it is not a tSURF routine.
    call surf_init_all_instances ! This module
    !
    call dump_all_wavefunctions('initial',ens_initial_wf_dump_prefix) ! spherical_tdse_io.f90
    !
    call TimerStop('Initialization')
    write (out,"(/'Done with initialization'/)")
    call TimerReport
    !
    select case (task)
      case default
        write (out,"('Task ',a,' is not recognized')") trim(task)
        stop 'spherical_tdse - bad task'
      case ('real time')
        write (out,"(/'Real-time propagation'/)")
        call choose_rotation_code  ! This module
        call initialize_field      ! spherical_tdse_field.f90
        call propagation           ! spherical_tdse_propagate.f90
      case ('imaginary time')
        write (out,"(/'Imaginary-time propagation, A = ',g24.13/)") vp_scale
        call flush_wrapper(out)
        if (ensemble_size>1) stop 'spherical_tdse: Imaginary-time not supported for ensembles'
        call imaginary_propagation_test(verbose,wfns_l(1),wfns_r(1),apot=vp_scale,dt=dt,nstep=timesteps)
    end select
    if (dt_subdivision/='off') then
      write (out,"(/'Largest number of microsteps per time step was ',i0)") max_subdivision
      write (out,"( '    Number of time steps requiring subdivision ',i0)") timesteps_subdivided
      write (out,"( '                    Total number of microsteps ',i0)") timesteps_micro
      write (out,"( '    Average number of microsteps per time step ',f0.3/)") real(timesteps_micro,kind=rk)/timesteps
    end if
    !
    call report_adaptive ! This module
    !
    call flush_wrapper(out)
    !
    tsurf_report_loop: do ie=1,ensemble_size
      call sts_report(tsurfs(ie),'In TDSE')
    end do tsurf_report_loop
    !
    call dump_all_wavefunctions('final',ens_final_wf_dump_prefix) ! spherical_tdse_io.f90
    !
    call reconstruct_all_left  ! This module
    !
    call analyze_and_surf_all  ! This module
    !
    call TimerStop('start')
    call TimerReport
    !
    call nt_finalize
  end subroutine start
end module spherical_tdse
!
program main
  use spherical_tdse
  call start
end program main
