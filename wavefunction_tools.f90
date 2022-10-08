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
!  Sundry tools for manipulating wavefunction
!
!  The routines in this modules are shared-memory parallelized (OpenMP). They
!  are aware of a possible distributed-memory parallelization, but rely on the 
!  upstream routines to set up the distributed-memory calculations.
!
module wavefunction_tools
  use accuracy
  use constants
  use hacks
  use timer
  use spherical_data
  use spherical_data_initialize
  use tridiagonal_tools
  use node_tools
  use sort_tools
  use lapack
  use math
  implicit none
  private
  public wt_atomic_cache_prefix, wt_iterative_improvement, wt_disable_orthogonalization
  public wt_max_solution_iterations
  public wt_atomic_solutions, wt_one_atomic_solution
  public wt_normalize, wt_energy, wt_dipole
  public wt_adaptive_r_buffer
  public wt_reset_lrmax, wt_update_lrmax
  public wt_resize
  public wt_r2r_scale
  public rcsid_wavefunction_tools
  !
  character(len=clen), save :: rcsid_wavefunction_tools = "$Id: wavefunction_tools.f90,v 1.62 2022/10/06 17:14:31 ps Exp ps $"
  !
  integer, parameter  :: iu_temp                = 24  ! An arbitrary unit number, which can be used here
  integer, parameter  :: wt_adaptive_r_buffer   = 8   ! Maximum number of points to examine for adaptive nradial determination
  character(len=clen) :: wt_atomic_cache_prefix = ' ' ! Atomic solution for each L will be cached
                                                      ! here; these quantities are expensive to
                                                      ! recompute, especially if we run a 2D simulation
                                                      ! Blank means no caching
  logical             :: wt_iterative_improvement = .true.
                                                      ! Set to .true. if the atomic eigenstates are to be 
                                                      ! iteratively refined
  logical             :: wt_disable_orthogonalization = .false.
                                                      ! Set to .true. to skip explicit Gram-Schmidt orthogonalization
                                                      ! for the atomic solutions
  integer(ik)         :: wt_max_solution_iterations = 20_ik
                                                      ! Max number of iterations in trying to find the eigenvectors
                                                      ! in wt_one_atomic_solution()
  !
  contains
  !
  !  Calculate a block of atomic solutions. We return the non-zero part of the
  !  eigenvectors only, instead of the much larger total wavefunction.
  !
  subroutine wt_atomic_solutions(verbose,lval,eval,evec)
    integer(ik), intent(in)  :: verbose     ! Reporting level
    integer(ik), intent(in)  :: lval        ! Desired angular momentum
    complex(rk), intent(out) :: eval(:)     ! Eigenvalues of the Hamiltonian matrix block for L=lval
    complex(rk), intent(out) :: evec(:,:,:) ! Left (:,:,1) and right (:,:,2) eigenvectors  
    !
    character(len=clen) :: cache
    integer(ik)         :: ios
    !
    call TimerStart('Atomic solutions')
    if (size(eval)/=sd_nradial .or. ubound(evec,1)/=sd_nradial .or. ubound(evec,2)/=sd_nradial .or. ubound(evec,3)/=2 &
        .or. lval<0 .or. lval>sd_lmax) then
      stop 'wavefunction_tools%wt_atomic_solutions - bad arguments'
    end if
    !
    !  Try to load from cache first
    !
    if (wt_atomic_cache_prefix/=' ') then
      !
      !  This function may be called from a parallel region; make sure it
      !  does not try to access the same file from multiple threads!
      !
      !$omp critical
      write (cache,"(a,'-L=',i5.5)") trim(wt_atomic_cache_prefix), lval
      open (iu_temp,action='read',status='old',form='unformatted',iostat=ios,file=trim(cache))
      if (ios==0) then
        read (iu_temp) eval, evec
        close (iu_temp)
      end if
      !$omp end critical
      if (ios==0) then
        call TimerStop('Atomic solutions')
        return
      end if
      !
      ! Open was unsuccessful; proceed through calculation, and try to save the 
      ! results once we have them
      !
    end if
    !
    !  No cached result available; compute it from scratch. 
    !
    !  Begin by computing the guess eigendecomposition using LAPACK.
    !  Unfortunately, LAPACK sometimes produces very inaccurate results for general
    !  complex matrices, so that these eigenvalues/eigenvectors cannot be used as is
    !
    call wt_atomic_solution_guess(verbose,lval,eval,evec)
    !
    call wt_atomic_solution_biorthogonalize(evec)
    !
    if (wt_iterative_improvement) then
      call wt_atomic_solution_improve(verbose,lval,eval,evec)
    end if
    !
    call wt_atomic_solution_verify(lval,evec)
    !
    if (wt_atomic_cache_prefix/=' ') then
      !$omp critical
      open (iu_temp,action='write',status='new',form='unformatted',iostat=ios,file=trim(cache))
      if (ios/=0) then
        write (out,"('Error ',i0,' creating new cache file ',a)") ios, trim(cache)
        write (out,"('Skipping wavefunction save for L=',i0,' and continuing')") lval
      else
        write (iu_temp) eval, evec
        close (iu_temp)
      end if
      !$omp end critical
    end if
    call TimerStop('Atomic solutions')
  end subroutine wt_atomic_solutions
  !
  !  Compute a single atomic solution using inverse iterations.
  !  Requires approximate energy of the solution as the starting guess.
  !
  subroutine wt_one_atomic_solution(verbose,lval,eval,evec)
    integer(ik), intent(in)    :: verbose   ! Reporting level
    integer(ik), intent(in)    :: lval      ! Desired angular momentum
    complex(rk), intent(inout) :: eval      ! In: guess eigenvalue
                                            ! Out: improved eigenvalue
    complex(rk), intent(out)   :: evec(:,:) ! Left (:,1) and right (:,2) eigenvectors  
    !
    real(rk)    :: guess_buf(sd_nradial)
    complex(rk) :: eval_new, delta
    integer(ik) :: iter
    logical     :: converged
    !
    call TimerStart('Atomic solution (single)')
    if (ubound(evec,1)/=sd_nradial .or. ubound(evec,2)/=2 .or. lval<0 .or. lval>sd_lmax) then
      stop 'wavefunction_tools%wt_one_atomic_solution - bad arguments'
    end if
    !
    !  Start with a random real guess. Don't forget to normalize!
    !
    call random_number(guess_buf)
    evec(:,1) = guess_buf
    call random_number(guess_buf)
    evec(:,2) = guess_buf
    call wt_normalize_one_atomic_solution(evec)
    !
    write (out,"('Solving for atomic eigenvector L=',i0,' E(guess)=',2(g26.12,1x),'using inverse iterations.')") lval, eval
    converged = .false.
    inverse_iteration_passes: do iter=1,wt_max_solution_iterations
      eval_new = eval
      call wt_improve_one_atomic_eigenstate(verbose,lval,eval_new,evec)
      delta = eval_new - eval
      eval  = eval_new
      if (verbose>=0) then
        write (out,"(1x,'Iteration ',i0,' E= ',2(g35.24e3,1x),' change = ',2(g15.6e3,1x))") iter, eval, delta
      end if
      if (abs(delta)<=100._rk*spacing(abs(eval_new))) then
        !
        !  We are converged, but do couple more iterations for a good measure
        !
        eval_new = eval
        call wt_improve_one_atomic_eigenstate(verbose,lval,eval_new,evec)
        call wt_improve_one_atomic_eigenstate(verbose,lval,eval_new,evec)
        delta = eval_new - eval
        eval  = eval_new
        if (verbose>=0) then
          write (out,"(1x,'Final iteration E= ',2(g35.24e3,1x),' change = ',2(g15.6e3,1x))") eval, delta
        end if
        converged = .true.
        exit inverse_iteration_passes
      end if
    end do inverse_iteration_passes
    !
    write (out,"('Final energy ',2(g35.24e3,1x))") eval
    if (.not.converged) then
      write (out,"('Inverse iterations failed to converge after ',i0,' passes.')") wt_max_solution_iterations
      write (out,"('Continuing with possibly unconverged eigenvectors.')")
    end if
    !
    call TimerStop('Atomic solution (single)')
  end subroutine wt_one_atomic_solution
  !
  subroutine wt_atomic_solution_guess(verbose,lval,eval,evec)
    integer(ik), intent(in)  :: verbose     ! Reporting level
    integer(ik), intent(in)  :: lval        ! Desired angular momentum
    complex(rk), intent(out) :: eval(:)     ! Eigenvalues of the Hamiltonian matrix block for L=lval
    complex(rk), intent(out) :: evec(:,:,:) ! Left (:,:,1) and right (:,:,2) eigenvectors  
    !
    real(rk), allocatable     :: tmat(:,:)
    complex(lrk), allocatable :: evec_guess(:,:,:), eval_guess(:)
    real(rk), allocatable     :: eval_guess_real(:)
    integer(ik)               :: order(sd_nradial)
    integer(ik)               :: alloc, ipt
    !
    call TimerStart('Atomic solutions: Guess')
    allocate (tmat(sd_nradial,sd_nradial),evec_guess(sd_nradial,sd_nradial,2), &
              eval_guess(sd_nradial),eval_guess_real(sd_nradial),stat=alloc)
    if (alloc/=0) stop 'wavefunction_tools%wt_atomic_solution_guess - no memory for tmat'
    !
    if (lval==0) then
      call sd_expand_implicit_operator(sd_d2n_l0,sd_m2n_l0,sd_m2nf_l0,sd_m2np_l0,tmat)
    else
      call sd_expand_implicit_operator(sd_d2n_lx,sd_m2n_lx,sd_m2nf_lx,sd_m2np_lx,tmat)
    end if
    evec_guess(:,:,1) = cmplx((-1._rk/(2._rk*electron_mass))*tmat,kind=lrk)
    !
    deallocate (tmat)
    !
    add_potential: do ipt=1,sd_nradial
      evec_guess(ipt,ipt,1) = evec_guess(ipt,ipt,1) + real(sd_pottab(ipt,1,lval),kind=lrk)
    end do add_potential
    if (sd_capped) then
      add_cap: do ipt=sd_cap_start,sd_nradial
        evec_guess(ipt,ipt,1) = evec_guess(ipt,ipt,1) + cmplx(sd_captab(ipt),kind=lrk)
      end do add_cap
    end if
    !
    if (verbose>=2) then
      write (out,"('       Largest matrix element of H for L=',i0,' is ',g24.13)") &
             lval, maxval(abs(evec_guess(:,:,1)))
      write (out,"('Largest deviation from hermiticity for L=',i0,' is ',g24.13)") &
             lval, maxval(abs(evec_guess(:,:,1)-conjg(transpose(evec_guess(:,:,1)))))
    end if
    !
    call lapack_geev(sd_nradial,evec_guess,eval_guess)
    eval_guess_real = real(eval_guess,kind=rk)
    call order_keys(eval_guess_real,order)
    eval = eval_guess(order)
    !
    !  LAPACK _GEEV conjugates the left eigenvectors. Let's begin by undoing the damage.
    !
    evec(:,:,1) = conjg(evec_guess(:,order,1))
    evec(:,:,2) =       evec_guess(:,order,2)
    !
    deallocate (evec_guess,eval_guess,eval_guess_real)
    call TimerStop('Atomic solutions: Guess')
  end subroutine wt_atomic_solution_guess
  !
  subroutine wt_atomic_solution_biorthogonalize(evec)
    complex(rk), intent(inout) :: evec(:,:,:) ! Left (:,:,1) and right (:,:,2) eigenvectors  
    !
    integer(ik)  :: imo, jmo
    complex(rk)  :: scl
    ! complex(rk)  :: fact_l(sd_nradial), fact_r(sd_nradial)
    !
    if (wt_disable_orthogonalization) return
    call TimerStart('Atomic solutions: Orthogonalize')
    !
    !  LAPACK _GEEV routines do not give biorthogonal eigenvectors
    !  in the present of multiple roots. Although this should not
    !  occur in our case, it is better to be safe than sorry. Since
    !  we also do not like LAPACK's normalization convention for the
    !  eigenvectors (Cartesian norm = 1 separately for the left and
    !  right vectors), we have to do a bit of work here
    !
    normalize_eigenvectors: do imo=1,sd_nradial
      !
      !  Enforce bi-orthogonality to all lower eigenvectors,
      !  using a variation of the Gram-Schmidt process.
      !  We assume that all lower eigenvectors are already normalized to
      !  our conventions (L.R = 1)
      !
      orthogonalize_lower: do jmo=1,imo-1
        ! Right eigenvector
        scl = sum(evec(:,jmo,1)*evec(:,imo,2))
        evec(:,imo,2) = evec(:,imo,2) - scl*evec(:,jmo,2)
        ! Left eigenvector
        scl = sum(evec(:,jmo,2)*evec(:,imo,1))
        evec(:,imo,1) = evec(:,imo,1) - scl*evec(:,jmo,1)
      end do orthogonalize_lower
      ! The matrix-vector code below should be faster; in fact, it is
      ! actually (much) slower. Oops.
      ! fact_l(:imo-1) = matmul(evec(:,imo,2),evec(:,:imo-1,1))
      ! fact_r(:imo-1) = matmul(evec(:,imo,1),evec(:,:imo-1,2))
      ! evec(:,imo,1) = evec(:,imo,1) - matmul(evec(:,:imo-1,1),fact_r(:imo-1))
      ! evec(:,imo,2) = evec(:,imo,2) - matmul(evec(:,:imo-1,2),fact_l(:imo-1))
      !
      call wt_normalize_one_atomic_solution(evec(:,imo,:))
    end do normalize_eigenvectors
    call TimerStop('Atomic solutions: Orthogonalize')
  end subroutine wt_atomic_solution_biorthogonalize
  !
  !  Normalize right eigenvector to Euclidian norm of 1,
  !  simultaneously making the largest coefficient positive real.
  !
  subroutine wt_normalize_one_atomic_solution(vecx)
    complex(rk), intent(inout) :: vecx(:,:)
    !
    real(rk)    :: nrm
    complex(rk) :: scl
    integer(ik) :: imax
    !
    if (ubound(vecx,1)/=sd_nradial .or. ubound(vecx,2)/=2) then
      stop 'wavefunction_tools%wt_normalize_one_atomic_solution - bad arguments'
    end if
    !
    !  Normalize right eigenvector. Iterative refinement tends to produce very large eigenvectors,
    !  with norms which can easily exceed the dynamics range of the underlying real type. We therefore
    !  have to be very careful here to avoid overflowing our data types.
    !
    !  Start by making the largest element of the right eigenvector equal to 1.
    !
    imax      = maxloc(abs(vecx(:,2)),dim=1)
    scl       = 1._rk / vecx(imax,2)
    vecx(:,2) = scl * vecx(:,2)
    !
    !  Do the same for the left eigenvector
    !
    imax      = maxloc(abs(vecx(:,1)),dim=1)
    scl       = 1._rk / vecx(imax,1)
    vecx(:,1) = scl * vecx(:,1)
    !
    !  Now set the 2-norm of the right eigenvector to 1.
    !
    nrm       = sum(abs(vecx(:,2))**2)
    vecx(:,2) = (1._rk/sqrt(nrm)) * vecx(:,2)
    !
    !  Finally, normalize the left eigenvector
    !
    scl       = sum(vecx(:,1)*vecx(:,2))
    vecx(:,1) = (1._rk/scl) * vecx(:,1)
    !
    !  Debugging code: check for NaNs. Should not happen, but being paranoid costs us very little
    !
    ! if (any(vecx/=vecx)) then
    !   write (out,"(/'FATAL: atomic eigenvectors contain NaNs in wt_normalize_one_atomic_solution')")
    !   write (out,"('Left eigenvector:')")
    !   write (out,"(6(1x,g24.14))") vecx(:,1)
    !   write (out,"('Right eigenvector:')")
    !   write (out,"(6(1x,g24.14))") vecx(:,2)
    !   stop 'wavefunction_tools%wt_normalize_one_atomic_solution - got NaN'
    ! end if
  end subroutine wt_normalize_one_atomic_solution
  !
  subroutine wt_atomic_solution_improve(verbose,lval,eval,evec)
    integer(ik), intent(in)    :: verbose     ! Reporting level
    integer(ik), intent(in)    :: lval        ! Desired angular momentum
    complex(rk), intent(inout) :: eval(:)     ! Eigenvalues of the Hamiltonian matrix block for L=lval
    complex(rk), intent(inout) :: evec(:,:,:) ! Left (:,:,1) and right (:,:,2) eigenvectors  
    !
    integer(ik) :: ipass
    integer(ik) :: iev
    complex(rk) :: eval_initial(size(eval))
    !
    call TimerStart('Atomic solutions: Improve')
    eval_initial = eval
    improvement_passes: do ipass=1,1
      scan_eigenvalues: do iev=1,sd_nradial
        call wt_improve_one_atomic_eigenstate(verbose,lval,eval(iev),evec(:,iev,:))
        if (verbose>=2) then
          write (out,"('Iterative improvement pass ',i0,' L= ',i0,' I= ',i0,' E= ',2g26.18)") ipass, lval, iev, eval(iev)
        end if
      end do scan_eigenvalues
      call wt_atomic_solution_biorthogonalize(evec)
    end do improvement_passes
    !
    if (verbose>=1) then
      eval_initial = eval_initial - eval
      iev = maxloc(abs(eval_initial),dim=1)
      write (out,"(/'Iterative update of L=',i0,' solutions complete')") lval
      write (out,"('        ground-state energy is ',g34.22,1x,g34.22,' change = ',g18.7,1x,g18.7)") eval(1), eval_initial(1)
      write (out,"('   most affected eigenvalue is ',g34.22,1x,g34.22,' change = ',g18.7,1x,g18.7)") eval(iev), eval_initial(iev)
    end if
    !
    call TimerStop('Atomic solutions: Improve')
  end subroutine wt_atomic_solution_improve
  !
  !  Improve a single eigenstate using inverse iterations
  !
  subroutine wt_improve_one_atomic_eigenstate(verbose,lval,eval,evec)
    integer(ik), intent(in)    :: verbose   ! Reporting level
    integer(ik), intent(in)    :: lval      ! Desired angular momentum
    complex(rk), intent(inout) :: eval      ! Eigenvalue to improve
    complex(rk), intent(inout) :: evec(:,:) ! Left (:,1) and right (:,2) eigenvectors to improve
    !
    integer(ik) :: iter_evec, iter_eval, imax
    logical     :: fail(2)
    logical     :: upd_nan(3)
    real(rk)    :: mn   (sd_nradial,3), mt   (sd_nradial,3)
    complex(rk) :: leqn (sd_nradial,3), leqt (sd_nradial,3)
    complex(rk) :: leqnf(sd_nradial,m3d_dc_size)
    complex(rk) :: leqtf(sd_nradial,m3d_dc_size)
    logical     :: leqnp(sd_nradial), leqtp(sd_nradial)
    complex(rk) :: tmp (sd_nradial)
    complex(rk) :: vecx(sd_nradial,2)
    complex(rk) :: sclx      ! Largest element of vecx(:,2)
    complex(rk) :: scr (sd_nradial,m3d_sc_size)
    complex(rk) :: delta_m1  ! Correction to the eigenvalue^{-1}
    complex(rk) :: corr      ! Correction to the eigenvalue
    complex(rk) :: eval0     ! Eigenvalue from the previous iteration
    !
    improve_eval: do iter_eval=1,5
      eval0 = eval
      call build_right_linear_system(eval,leqn,leqnf,leqnp,mn,fail(2))
      call build_left_linear_system (eval,leqt,leqtf,leqtp,mt,fail(1))
      if (any(fail)) then 
        if (iter_eval<=1 .or. verbose>2) then
          write (out,"('WARNING: wt_improve_one_atomic_eigenstate: Update iteration ',i0,' failed, leave solutions unchanged')") &
                 iter_eval
        end if
        return
      end if
      !
      !  There is non-negligible cost to building linear system solutions;
      !  perform a few eigenvector updates before updating the eigenvalue
      !
      !  There is a possibility of blowing the dynamic range of the floating-
      !  point representation during the inverse iterations. We need to to be
      !  a little careful with the results.
      !
      improve_evec: do iter_evec=1,3
        ! Inverse iteration: right vector
        call m3d_multiply(mn,evec(:,2),tmp)
        call m3d_solve(leqn,leqnf,leqnp,tmp,vecx(:,2),scr)
        ! Inverse iteration: left vector
        call m3d_solve(leqt,leqtf,leqtp,evec(:,1),tmp,scr)
        call m3d_multiply(mt,tmp,vecx(:,1))
        ! Adjust the range of vecx(:,1) to prevent overflows
        ! If we already had an infinity somewhere, it should become a NaN, and trigger an error exit below.
        imax = maxloc(abs(vecx(:,1)),dim=1)
        sclx = 1._rk / vecx(imax,1)
        vecx(:,1) = sclx * vecx(:,1)
        ! Adjust the range of vecx(:,2) to prevent overflows
        imax = maxloc(abs(vecx(:,2)),dim=1)
        sclx = 1._rk / vecx(imax,2)
        vecx(:,2) = sclx * vecx(:,2)
        ! Compute correction to the eigenvalue. Don't forget to apply the scale factor later!
        delta_m1  = sum(evec(:,1)*vecx(:,2))
        corr      = sclx/delta_m1
        !
        call wt_normalize_one_atomic_solution(vecx)
        !
        ! Make sure none of the updated quantities are NaN - this could happen if we have
        ! insufficient dynamic range.
        !
        upd_nan(1) = isnan_wrapper(corr)
        upd_nan(2) = any(isnan_wrapper(vecx(:,1)))
        upd_nan(3) = any(isnan_wrapper(vecx(:,2)))
        !
        if (verbose>2) then
          write (out,"('Update ',i0,',',i0,' upd_nan = ',3l1,' eval = ',2(g32.20,1x),' corr. = ',2(g32.20,1x))") &
                 iter_eval, iter_evec, upd_nan, eval, corr
        end if
        !
        if (any(upd_nan)) then
          if (iter_eval<=1 .or. verbose>2) then
            write (out,"('WARNING: wt_improve_one_atomic_eigenstate: Update ',i0,',',i0,' failed, taking last update.')") &
                   iter_eval, iter_evec
          end if
          return
        end if
        ! Store the updated eigenvectors
        evec(:,:) = vecx(:,:)
        ! Update eigenvalue, in case we'll need to do an early exit later
        eval = eval0 + corr
      end do improve_evec
    end do improve_eval
    !
    contains
    subroutine build_right_linear_system(eps,leqn,leqnf,leqnp,mn,fail)
      complex(rk), intent(in)  :: eps        ! Current eigenvalue guess
      complex(rk), intent(out) :: leqn (:,:) ! Linear system
      complex(rk), intent(out) :: leqnf(:,:) ! Factorization of the linear system
      logical,     intent(out) :: leqnp(  :) ! Pivot list for the linear system
      real(rk), intent(out)    :: mn   (:,:) ! Scaling factor for the RHS
      logical, intent(out)     :: fail
      !
      complex(rk) :: vtmp(  sd_nradial)
      !
      !  Prepare linear system matrices for solving right linear system: 
      !
      !     [-(1/(2m))M^{-1}Delta+V-eps] Y = Z.
      !
      !  This is done as:
      !
      !     [-(1/(2m))Delta+M(V-eps)] Y = M Z
      !
      vtmp(:) = sd_pottab(:,1,lval) - eps
      if (sd_capped) then
        vtmp(sd_cap_start:) = vtmp(sd_cap_start:) + sd_captab(sd_cap_start:)
      end if
      if (lval==0) then
        call m3d_right_scale(sd_m2n_l0,vtmp,leqn)
        leqn(:,:) = leqn + (-1._rk/(2._rk*electron_mass))*sd_d2n_l0
        mn  (:,:) = sd_m2n_l0
      else
        call m3d_right_scale(sd_m2n_lx,vtmp,leqn)
        leqn(:,:) = leqn + (-1._rk/(2._rk*electron_mass))*sd_d2n_lx
        mn  (:,:) = sd_m2n_lx
      end if
      call m3d_decompose(leqn,leqnf,leqnp,fail=fail)
    end subroutine build_right_linear_system
    !
    subroutine build_left_linear_system(eps,leqt,leqtf,leqtp,mt,fail)
      complex(rk), intent(in)  :: eps        ! Current eigenvalue guess
      complex(rk), intent(out) :: leqt (:,:) ! Linear system
      complex(rk), intent(out) :: leqtf(:,:) ! Factorization of the linear system
      logical,     intent(out) :: leqtp(  :) ! Pivot list for the linear system
      real(rk), intent(out)    :: mt   (:,:) ! Scaling factor for the solution
      logical, intent(out)     :: fail
      !
      complex(rk) :: vtmp(  sd_nradial)
      !
      !  Prepare linear system matrices for solving left linear system: 
      !
      !     [-(1/(2m))Delta^T M^{-T}+V-eps] Y = Z.
      !
      !  This is done as:
      !
      !     [-(1/(2m))Delta^T+(V-eps)M^T] Q = Z
      !     Y = M^T Q
      !
      vtmp(:) = sd_pottab(:,1,lval) - eps
      if (sd_capped) then
        vtmp(sd_cap_start:) = vtmp(sd_cap_start:) + sd_captab(sd_cap_start:)
      end if
      if (lval==0) then
        call m3d_left_scale(vtmp,sd_m2t_l0,leqt)
        leqt(:,:) = leqt(:,:) + (-1._rk/(2._rk*electron_mass))*sd_d2t_l0
        mt  (:,:) = sd_m2t_l0
      else
        call m3d_left_scale(vtmp,sd_m2t_lx,leqt)
        leqt(:,:) = leqt(:,:) + (-1._rk/(2._rk*electron_mass))*sd_d2t_lx
        mt  (:,:) = sd_m2t_lx
      end if
      call m3d_decompose(leqt,leqtf,leqtp,fail=fail)
    end subroutine build_left_linear_system
  end subroutine wt_improve_one_atomic_eigenstate
  !
  subroutine wt_atomic_solution_verify(lval,evec)
    integer(ik), intent(in) :: lval
    complex(rk), intent(in) :: evec(:,:,:) ! Left (:,:,1) and right (:,:,2) eigenvectors  
    !
    integer(ik)              :: imo, jmo, alloc
    real(rk)                 :: nrm
    complex(rk), allocatable :: norm(:,:)
    integer(ik)              :: warning_count
    complex(rk)              :: err, err_max
    integer(ik), parameter   :: max_warnings = 5
    !
    call TimerStart('Atomic solutions: Verify')
    allocate (norm(sd_nradial,sd_nradial),stat=alloc)
    if (alloc/=0) stop 'wavefunction_tools%wt_atomic_solution_verify - allocation failed'
    !
    norm          = matmul(transpose(evec(:,:,1)),evec(:,:,2))
    warning_count = 0
    err_max       = 0
    test_norm_j: do jmo=1,sd_nradial
      test_norm_i: do imo=1,sd_nradial
        nrm = 0
        if (imo==jmo) nrm = 1
        err = norm(imo,jmo)-nrm
        if (abs(err)>1e5_rk*spacing(1._rk)) then
          warning_count = warning_count + 1
          if (abs(err)>abs(err_max)) err_max = err
          if (warning_count<=max_warnings) then
            write (out,"('WARNING: _GEEV eigenvectors ',i0,',',i0,' L=',i0,' are not biorthogonal. Error =',2(1x,g24.13))") &
                   imo, jmo, lval, err
          end if
        end if
      end do test_norm_i
    end do test_norm_j
    !
    if (warning_count>max_warnings) then
      write (out,"('WARNING: ',i0,' _GEEV biorthogonality warnings for L=',i0,' suppressed. Max error was ',2(1x,g24.13))") &
             warning_count-max_warnings, lval, err_max
    end if
    !
    deallocate (norm)
    call TimerStop('Atomic solutions: Verify')
  end subroutine wt_atomic_solution_verify
  !
  !  Choose wavefunction lmax and nradial. This is only done once; do not bother
  !  parallizing this routine.
  !
  subroutine wt_reset_lrmax(wfn)
    type(sd_wfn), intent(inout) :: wfn   ! Wavefunction to examine
    !
    integer(ik) :: lval, mval, ir
    real(rk)    :: max_val_l(0:sd_lmax)
    real(rk)    :: max_val_r(sd_nradial)
    real(rk)    :: eps_l, eps_r          ! For the initial reset, use tigher tolerances.
                                         ! This will stop the wavefunction from growing spuriously later on
                                         ! due to the edge-reflection effects
    !
    wfn%lmax    = sd_lmax
    wfn%nradial = sd_nradial
    if (.not.sd_adaptive) return
    !
    call TimerStart('WT reset_lrmax')
    !
    max_val_l = 0
    max_val_r = 0
    !
    !  Collect maximum wavefunction values across L channels
    !
    sense_loop_m: do mval=sd_mmin,sd_mmax
      l_sense_loop_l: do lval=abs(mval),sd_lmax
        max_val_l(lval) = max(max_val_l(lval),maxval(abs(wfn%wfn(:,:,lval,mval))))
      end do l_sense_loop_l
      sense_loop_r: do ir=1,sd_nradial
        max_val_r(ir) = max(max_val_r(ir),maxval(abs(wfn%wfn(ir,:,abs(mval):sd_lmax,mval))))
      end do sense_loop_r
    end do sense_loop_m
    !
    if (sd_adaptive_l) then
      eps_l = 1e-2_rk * wfn%sd_tol_l
      choose_lmax: do lval=sd_lmax,0,-1
        if (max_val_l(lval)>eps_l) then
          wfn%lmax = min(lval+1,sd_lmax)
          exit choose_lmax
        end if
      end do choose_lmax
    end if
    if (sd_adaptive_r) then
      eps_r = 1e-2_rk * wfn%sd_tol_r
      choose_rmax: do ir=sd_nradial,1,-1
        if (max_val_r(ir)>eps_r) then
          wfn%nradial = min(ir+1,sd_nradial)
          exit choose_rmax
        end if
      end do choose_rmax
    end if
    call TimerStop('WT reset_lrmax')
  end subroutine wt_reset_lrmax
  !
  !  Update wavefunction lmax if necessary
  !
  subroutine wt_update_lrmax(verbose,wfn)
    integer(ik), intent(in)     :: verbose ! Verbosity level
    type(sd_wfn), intent(inout) :: wfn     ! Wavefunction to examine
    !
    integer(ik) :: mval, mmin, mmax, ir, lval, old_lmax
    integer(ik) :: nrtop, old_nradial
    real(rk)    :: max_val_r, max_val_l(2), max_t
    real(rk)    :: max_r_buf(wt_adaptive_r_buffer+1)
    !
    if (.not.sd_adaptive) return
    call TimerStart('WT update_lrmax')
    !
    !  Determine the radial extent to examine
    !
    nrtop = min(sd_nradial,wfn%nradial+wt_adaptive_r_buffer)
    !
    !  Adjust angular momentum. We may adjust both up and down!
    !
    if (sd_adaptive_l) then
      !
      !  max_val_l(1)=max_l1 is the largest element for L=lmax
      !  max_val_l(2)=max_l2 is the largest element for L=lmax-1
      !
      !  Older versions of Gfortran have a bug in max: reduction for
      !  arrays; we'll work around it by using explicit variables.
      !
      max_val_l = 0
      !
      if (nts%dist(wfn%lmax)==0) then
        mmin  = max(sd_mmin,-wfn%lmax)
        mmax  = min(sd_mmax, wfn%lmax)
        max_t = 0
        !$omp parallel do if(mmax>mmin) default(none) &
        !$omp& shared(mmin,mmax,sd_nspin,wfn,nrtop,nts) &
        !$omp& private(mval) reduction(max:max_t)
        l_sense_loop_m1: do mval=mmin,mmax
          if (abs(mval)>wfn%lmax) cycle l_sense_loop_m1
          max_t = max(max_t,maxval(abs(wfn%wfn(:nrtop,1:sd_nspin,wfn%lmax,mval))))
        end do l_sense_loop_m1
        !$omp end parallel do
        max_val_l(1) = max_t
      end if
      if (wfn%lmax>0) then
        if (nts%dist(wfn%lmax-1)==0) then
          mmin  = max(sd_mmin,-(wfn%lmax-1))
          mmax  = min(sd_mmax, (wfn%lmax-1))
          max_t = 0
          !$omp parallel do if(mmax>mmin) default(none) &
          !$omp& shared(mmin,mmax,sd_nspin,wfn,nrtop,nts) &
          !$omp& private(mval) reduction(max:max_t)
          l_sense_loop_m2: do mval=mmin,mmax
            if (abs(mval)>wfn%lmax-1) cycle l_sense_loop_m2
            max_t = max(max_t,maxval(abs(wfn%wfn(:nrtop,1:sd_nspin,wfn%lmax-1,mval))))
          end do l_sense_loop_m2
          !$omp end parallel do
          max_val_l(2) = max_t
        end if
      end if
      call nt_max(max_val_l)
      !
      if (max_val_l(1)>wfn%sd_tol_l) then
        !
        !  Grow maximum L value by 1
        !
        if (verbose>0 .and. wfn%lmax<sd_lmax) then
          write (out,"('wt_update_lmax: L increased from ',i0,' max_val = ',2g24.12)") wfn%lmax, max_val_l
        end if
        wfn%lmax     = min(sd_lmax,wfn%lmax+1_ik)
        wfn%lmax_top = max(wfn%lmax_top,wfn%lmax)
      end if
      !
      if (all(max_val_l<=0.1_rk*wfn%sd_tol_l)) then
        !
        !  Shrink maximum L value by 1
        !
        if (verbose>0 .and. wfn%lmax>0) then
          write (out,"('wt_update_lmax: L decreased from ',i0,' max_val = ',2g24.12)") wfn%lmax, max_val_l
        end if
        old_lmax = wfn%lmax
        wfn%lmax = max(0_ik,wfn%lmax-1_ik)
        !
        !  If lmax decreased, clear the wavefunction array: we assume that parts of the
        !  wavefunction with L>wfn%lmax is rigorously zero - or we'll start picking up
        !  errors as we grow the array again.
        !
        if (old_lmax/=wfn%lmax .and. nts%dist(old_lmax)==0) then
          mmin = max(sd_mmin,-old_lmax)
          mmax = min(sd_mmax, old_lmax)
          !$omp parallel do if(mmax>mmin) default(none) &
          !$omp& private(mval) shared(mmin,mmax,wfn,sd_nspin,old_lmax)
          l_clear_loop_m: do mval=mmin,mmax
            if (abs(mval)>old_lmax) cycle l_clear_loop_m
            wfn%wfn(:,1:sd_nspin,old_lmax,mval) = 0
          end do l_clear_loop_m
          !$omp end parallel do
        end if
      end if
    end if
    !
    !  Adjust radial extent
    !
    if (sd_adaptive_r .and. wfn%nradial<sd_nradial) then
      mmin = max(sd_mmin,-wfn%lmax)
      mmax = min(sd_mmax, wfn%lmax)
      old_nradial = wfn%nradial
      r_sense_loop_r: do ir=old_nradial,nrtop
        max_val_r = 0
        !$omp parallel do collapse(2) default(none) &
        !$omp& shared(mmin,mmax,sd_nspin,wfn,nrtop,ir,nts) &
        !$omp& private(lval,mval) &
        !$omp& reduction(max:max_val_r)
        r_sense_loop_m: do mval=mmin,mmax
          r_sense_loop_l: do lval=0,wfn%lmax
            if (lval<abs(mval)) cycle r_sense_loop_l
            if (nts%dist(lval)>0) cycle r_sense_loop_l
            max_val_r = max(max_val_r,maxval(abs(wfn%wfn(ir,1:sd_nspin,lval,mval))))
          end do r_sense_loop_l
        end do r_sense_loop_m
        !$omp end parallel do
        max_r_buf(ir-old_nradial+1) = max_val_r
      end do r_sense_loop_r
      call nt_max(max_r_buf)
      !
      r_adjust_loop: do ir=old_nradial,nrtop
        max_val_r = max_r_buf(ir-old_nradial+1)
        if (max_val_r>wfn%sd_tol_r) then
          wfn%nradial = min(sd_nradial,ir+1)
        end if
        ! write (out,"('wfn%nradial updated to ',i8)") wfn%nradial
      end do r_adjust_loop
    end if
    call TimerStop('WT update_lrmax')
  end subroutine wt_update_lrmax
  !
  !  Resize wavefunction arrays. If the wavefunction is truncated, data will be discarded
  !  for larger values of the indices. If it is grown, new elements will be initialized to
  !  zero. There is no data conversion or resampling of any kind goind on.
  !
  !  wt_resize() does not need to be aware of distributed-memory parallization at this point.
  !
  subroutine wt_resize(wfn,nradial)
    type(sd_wfn), intent(inout) :: wfn     ! Wavefunction to reallocate
    integer(ik), intent(in)     :: nradial ! New nradial
    !
    complex(rk), allocatable :: tmp_wfn(:,:,:,:)
    integer(ik)              :: alloc
    integer(ik)              :: old_nr, min_nr
    integer(ik)              :: lv, mv
    !
    if (nradial==size(wfn%wfn,dim=1)) return
    !
    call TimerStart('WT resize')
    !
    old_nr = size(wfn%wfn,dim=1)
    min_nr = min(old_nr,nradial)
    allocate (tmp_wfn(old_nr,sd_nspin,0:sd_lmax,sd_mmin:sd_mmax),stat=alloc)
    if (alloc/=0) stop 'wavefunction_tools%wt_resize - allocate failed (1)'
    !
    !$omp parallel default(none) &
    !$omp& shared(wfn,tmp_wfn,alloc,old_nr,min_nr,sd_nspin,sd_lmax,sd_mmin,sd_mmax,nradial) &
    !$omp& private(lv,mv)
    !
    !  Make a copy
    !
    !$omp do collapse(2)
    copy_out_m: do mv=sd_mmin,sd_mmax
      copy_out_l: do lv=0,sd_lmax
        tmp_wfn(:,:,lv,mv) = wfn%wfn(:,:,lv,mv)
      end do copy_out_l
    end do copy_out_m
    !$omp end do
    !
    !$omp single
    deallocate (wfn%wfn,stat=alloc)
    if (alloc/=0) stop 'wavefunction_tools%wt_resize - deallocate failed (1)'
    allocate (wfn%wfn(nradial,sd_nspin,0:sd_lmax,sd_mmin:sd_mmax),stat=alloc)
    if (alloc/=0) stop 'wavefunction_tools%wt_resize - allocate failed (2)'
    !$omp end single
    !
    !$omp do collapse(2)
    copy_in_m: do mv=sd_mmin,sd_mmax
      copy_in_l: do lv=0,sd_lmax
        wfn%wfn(:min_nr,:,lv,mv) = tmp_wfn(:min_nr,:,lv,mv)
        wfn%wfn(min_nr+1:nradial,:,lv,mv) = 0
      end do copy_in_l
    end do copy_in_m
    !$omp end do
    !
    !$omp end parallel
    !
    deallocate(tmp_wfn,stat=alloc)
    if (alloc/=0) stop 'wavefunction_tools%wt_resize - deallocate failed (1)'
    call TimerStop('WT resize')
  end subroutine wt_resize
  !
  !  Reset wavefunction norm.
  !  We will set Cartesian norm of the right wavefunction to 1, while maintaining its phase.
  !  The left wavefunction will be adjusted so that <l|r> product is also 1.
  !  Keep in mind that the left wavefunction does not need to be conjugated!
  !
  subroutine wt_normalize(wfn_l,wfn_r,norm)
    type(sd_wfn), intent(inout) :: wfn_l   ! Left wavefunction
    type(sd_wfn), intent(inout) :: wfn_r   ! Right wavefunction
                                           ! Since we are dealing with (potentially) non-Hermitian
                                           ! operators, our wavefunctions always come as a pair of
                                           ! the left and right vectors.
    complex(rk), intent(out)    :: norm(2) ! norm(1) - Input wavefunction norm <l|r> 
                                           ! norm(2) - Cartesian norm of the input right wavefunction <r|r>
                                           ! norm(2) is guaranteed to be real.
    !
    integer(ik) :: lval, mval, sval, my_lmax
    real(rk)    :: sum_rr, scl_r
    complex(rk) :: sum_lr, scl_l
    !
    call TimerStart('WF normalize')
    sum_rr  = 0
    sum_lr  = 0
    my_lmax = max(wfn_l%lmax,wfn_r%lmax)
    !$omp parallel default(none) &
    !$omp& shared(sd_mmin,sd_mmax,my_lmax,sd_nspin,wfn_l,wfn_r,nts) &
    !$omp& private(mval,lval,sval) &
    !$omp& shared(sum_rr,sum_lr,scl_r,scl_l)
    sense_loop_m: do mval=sd_mmin,sd_mmax
      !$omp do reduction(+:sum_rr,sum_lr)
      sense_loop_l: do lval=abs(mval),my_lmax
        if (nts%dist(lval)>0) cycle sense_loop_l ! Hook for distributed-memory parallelization
        sense_loop_s: do sval=1,sd_nspin
          sum_rr = sum_rr + sum(abs(wfn_r%wfn(:,sval,lval,mval))**2)
          sum_lr = sum_lr + sum(wfn_l%wfn(:,sval,lval,mval)*wfn_r%wfn(:,sval,lval,mval))
        end do sense_loop_s
      end do sense_loop_l
      !$omp end do nowait
    end do sense_loop_m
    !$omp barrier
    !$omp single
    call nt_add(sum_rr)
    call nt_add(sum_lr)
    scl_r = 1._rk / sqrt(sum_rr)
    scl_l = sqrt(sum_rr) / sum_lr
    !$omp end single
    !$omp barrier
    scale_loop_m: do mval=sd_mmin,sd_mmax
      !$omp do
      scale_loop_l: do lval=abs(mval),my_lmax
        if (nts%dist(lval)>0) cycle scale_loop_l ! Hook for distributed-memory parallelization
        scale_loop_s: do sval=1,sd_nspin
          wfn_r%wfn(:,sval,lval,mval) = scl_r*wfn_r%wfn(:,sval,lval,mval)
          wfn_l%wfn(:,sval,lval,mval) = scl_l*wfn_l%wfn(:,sval,lval,mval)
        end do scale_loop_s
      end do scale_loop_l
      !$omp end do nowait
    end do scale_loop_m
    !$omp end parallel
    norm(1) = sum_lr
    norm(2) = sum_rr
    call TimerStop('WF normalize')
  end subroutine wt_normalize
  !
  !  Calculate expectation value of the Hamiltonian
  !
  !  For distributed-memory execution, this routine requires a border region, of
  !  at least one extra L value on each side. There are no checks of whether this
  !  is true!
  !
  subroutine wt_energy(wfn_l,wfn_r,apot,energy,norm)
    type(sd_wfn), intent(in) :: wfn_l     ! Left wavefunction 
    type(sd_wfn), intent(in) :: wfn_r     ! Right wavefunction 
    real(xk), intent(in)     :: apot      ! Current value of the vector-potential; assumed to be along Z
    complex(rk), intent(out) :: energy(2) ! Expectation value of the total energy;
                                          ! [1] = including all terms in the Hamiltonian
                                          ! [2] = excluding the CAP term; only collected if sd_capped is .true.
    complex(rk), intent(out) :: norm      ! Norm of the input wavefunction; <psi|psi>
    !
    integer(ik) :: lval, mval, my_lmax, nr
    complex(rk) :: hpsi(sd_nradial), tmp(sd_nradial)
    complex(rk) :: scr (sd_nradial,m3d_sc_size) 
    complex(rk) :: ap_energy, en_nocap
    real(rk)    :: a2_potential, a_factor
    !
    call TimerStart('WF energy')
    !
    if (sd_nspin/=1) stop 'wavefunction_tools%wt_energy - spinorbit not implemented'
    !
    !  Coupling coefficients and prefactors due to the vector potential
    !
    a_factor     = real((electron_charge/electron_mass) * apot,kind=rk)
    a2_potential = real(0.5_rk * (electron_charge**2 / electron_mass) * apot**2,kind=rk)
    !
    energy(:) = 0
    norm      = 0
    my_lmax   = max(wfn_l%lmax,wfn_r%lmax)
    nr        = max(wfn_l%nradial,wfn_r%nradial)
    !$omp parallel default(none) &
    !$omp& shared(sd_mmin,sd_mmax,my_lmax,nr) &
    !$omp& shared(sd_d2n_l0,sd_m2n_l0,sd_m2nf_l0,sd_m2np_l0) &
    !$omp& shared(sd_d2n_lx,sd_m2n_lx,sd_m2nf_lx,sd_m2np_lx) &
    !$omp& shared(sd_d1n_l0,sd_m1n_l0,sd_m1nf_l0,sd_m1np_l0) &
    !$omp& shared(sd_d1n_lx,sd_m1n_lx,sd_m1nf_lx,sd_m1np_lx) &
    !$omp& shared(sd_pottab,sd_capped,sd_captab,sd_cap_start) &
    !$omp& shared(a2_potential,a_factor,sd_rminus) &
    !$omp& shared(wfn_l,wfn_r,nts) &
    !$omp& private(mval,lval,tmp,hpsi,ap_energy,scr,en_nocap) &
    !$omp& reduction(+:energy,norm)
    energy(:) = 0
    norm      = 0
    sense_loop_m: do mval=sd_mmin,sd_mmax
      !$omp do
      sense_loop_l: do lval=abs(mval),my_lmax
        if (nts%dist(lval)>0) cycle sense_loop_l  ! Hook for distributed-memory parallelization
        !
        !  L-diagonal part of the Hamiltonian: field-free and A^2 terms
        !
        if (lval==0) then
          call m3d_multiply(sd_d2n_l0,wfn_r%wfn(:,1,lval,mval),tmp)
          call m3d_solve(sd_m2n_l0,sd_m2nf_l0,sd_m2np_l0,tmp,hpsi,scr)
        else
          call m3d_multiply(sd_d2n_lx,wfn_r%wfn(:,1,lval,mval),tmp)
          call m3d_solve(sd_m2n_lx,sd_m2nf_lx,sd_m2np_lx,tmp,hpsi,scr)
        end if
        hpsi(:nr) = (-0.5_rk/electron_mass)*hpsi(:nr)
        hpsi(:nr) = hpsi(:nr) + (a2_potential+sd_pottab(:nr,1,lval))*wfn_r%wfn(:nr,1,lval,mval)
        !
        en_nocap  = sum(wfn_l%wfn(:nr,1,lval,mval)*hpsi(:nr))
        energy(:) = energy(:) + en_nocap
        norm      = norm      + sum(wfn_l%wfn(:nr,1,lval,mval)*wfn_r%wfn(:nr,1,lval,mval))
        if (sd_capped .and. sd_cap_start<=nr) then
          energy(1) = energy(1) &
                    + sum(wfn_l%wfn(sd_cap_start:nr,1,lval,mval)*sd_captab(sd_cap_start:nr)*wfn_r%wfn(sd_cap_start:nr,1,lval,mval))
        end if
        !
        !  A.p terms; these only need to be calculated if vector potential is not zero
        !
        if (a_factor==0) cycle sense_loop_l
        !
        !  Act on the RHS with the radial-gradient part of Hmix
        !
        if (lval==0) then
          call m3d_multiply(sd_d1n_l0,wfn_r%wfn(:,1,lval,mval),tmp)
          call m3d_solve(sd_m1n_l0,sd_m1nf_l0,sd_m1np_l0,tmp,hpsi,scr)
        else
          call m3d_multiply(sd_d1n_lx,wfn_r%wfn(:,1,lval,mval),tmp)
          call m3d_solve(sd_m1n_lx,sd_m1nf_lx,sd_m1np_lx,tmp,hpsi,scr)
        end if
        !
        !  Act of the RHS with the (1/r) term in Hang
        !
        tmp(:nr) = sd_rminus(:nr)*wfn_r%wfn(:nr,1,lval,mval)
        !
        !  Assemble the field-coupling term
        !
        ap_energy = 0
        if (lval>abs(mval)) then ! We can couple down
          ap_energy = ap_energy &
                    + real(sd_laser_clm(lval,  mval),kind=rk)*sum(wfn_l%wfn(:nr,1,lval-1,mval)*(hpsi(:nr)+(lval  )*tmp(:nr)))
        endif
        if (lval<my_lmax) then ! We can couple up
          ap_energy = ap_energy &
                    + real(sd_laser_clm(lval+1,mval),kind=rk)*sum(wfn_l%wfn(:nr,1,lval+1,mval)*(hpsi(:nr)-(lval+1)*tmp(:nr)))
        end if
        energy(:) = energy(:) + (0._rk,1._rk)*a_factor*ap_energy
      end do sense_loop_l
      !$omp end do nowait
    end do sense_loop_m
    !$omp barrier
    !$omp end parallel
    call nt_add(energy)
    call nt_add(norm)
    call TimerStop('WF energy')
  end subroutine wt_energy
  !
  !  Calculate dipole expectation and dipole acceleration 
  !
  !  WARNING: Dipole acceleration is valid ONLY for multiplicative potentials
  !
  !  The calculated dipole velocity does not include the vector-potential contribution;
  !  this term is best evaluated later, in the laboratory frame. Specifically, we need
  !  to add:
  !
  !     -(e^2/m) A <L|R>
  !
  !  The calculated acceleration term does not include the time derivative of the
  !  vector potential; this term is best evaluated later, in the laboratory frame.
  !  Specifically, one needs to add:
  !
  !      -(e^2/m) (d A/d t) <L|R> = (e^2/m) E <L|R>
  !
  !  where E is the instantaneous value of the laser electric field.
  !
  subroutine wt_dipole(wfn_l,wfn_r,do_term,dipole)
    type(sd_wfn), intent(in) :: wfn_l           ! Left wavefunction 
    type(sd_wfn), intent(in) :: wfn_r           ! Right wavefunction 
    logical, intent(in)      :: do_term(3)      ! do_term(1) - .true. if dipole is needed
                                                ! do_term(2) - .true. if dipole velocity is needed
                                                ! do_term(3) - .true. if dipole acceleration is needed
    complex(rk), intent(out) :: dipole(3,3)     ! dipole(:,1) - <L|q r|R> expectation value
                                                ! dipole(:,2) - [(d/d t) <L|q r|R>] + (e^2/m) A <L|R>
                                                ! dipole(:,3) - [(d^2/d t^2) <L|q r|R>] + (e^2/m) (d A/d t) <L|R>
                                                ! All dipole terms are in the local (rotating) frame
    !
    integer(ik) :: l_left, m_left, l_right, m_right, m_op, ispin, my_lmax, nr
    complex(rk) :: dip_sph(-1:1), vel_sph(-1:1), acc_sph(-1:1), wgt_dip, wgt_acc
    complex(rk) :: wgt_vel
    real(rk)    :: wgt_ang
    ! Arrays needed for dipole velocity calculation
    complex(rk) :: grad_m(sd_nradial,sd_nspin)   !  (d/dr) Psi - ((l+1)/r)*psi
    complex(rk) :: grad_p(sd_nradial,sd_nspin)   !  (d/dr) Psi +     (l/r)*psi
    complex(rk) :: scr (sd_nradial,m3d_sc_size), tmp(sd_nradial)
    !
    !  Quick return if possible
    !
    if (.not.any(do_term)) then
      dipole = 0
      return
    end if
    !
    call TimerStart('WF dipole')
    my_lmax = max(wfn_l%lmax,wfn_r%lmax)
    nr      = max(wfn_l%nradial,wfn_r%nradial)
    !
    !  Compute spherical-tensor form for the dipole
    !
    dip_sph(:) = 0
    vel_sph(:) = 0
    acc_sph(:) = 0
    !$omp parallel default(none) &
    !$omp& shared(sd_mmin,sd_mmax,my_lmax,nr,sd_nspin,wfn_l,wfn_r,sd_rtab,sd_dvdrtab) &
    !$omp& shared(do_term,sd_rminus,nts) &
    !$omp& shared(sd_d1n_l0,sd_m1n_l0,sd_m1nf_l0,sd_m1np_l0) &
    !$omp& shared(sd_d1n_lx,sd_m1n_lx,sd_m1nf_lx,sd_m1np_lx) &
    !$omp& private(m_left,l_left,m_right,l_right,m_op,ispin,wgt_ang,wgt_dip,wgt_acc) &
    !$omp& private(grad_m,grad_p,scr,tmp) &
    !$omp& private(wgt_vel) &
    !$omp& reduction(+:dip_sph,vel_sph,acc_sph)
    dip_sph(:) = 0
    vel_sph(:) = 0
    acc_sph(:) = 0
    loop_m_right: do m_right=sd_mmin,sd_mmax
      !$omp do
      loop_l_right: do l_right=abs(m_right),my_lmax
        if (nts%dist(l_right)>0) cycle loop_l_right ! Hook for distributed-memory parallelization
        !
        !  Precompute factors needed for dipole-velocity evaluation
        !
        if (do_term(2)) then
          prepare_spin: do ispin=1,sd_nspin
            call prepare_gradients(nr,wfn_r%wfn(:,ispin,l_right,m_right),l_right,scr,tmp,grad_m(:,ispin),grad_p(:,ispin))
          end do prepare_spin
        end if ! do_term(1)
        !
        loop_m_left: do m_left=max(sd_mmin,m_right-1),min(sd_mmax,m_right+1)
          loop_l_left: do l_left=l_right-1,l_right+1
            if (l_left<abs(m_left)) cycle loop_l_left
            if (l_left>my_lmax) cycle loop_l_left
            m_op    = m_left - m_right
            wgt_dip = 0
            wgt_vel = 0
            wgt_acc = 0
            loop_spin: do ispin=1,sd_nspin
              if (do_term(1)) then ! Dipole expectation
                wgt_dip = wgt_dip + sum(wfn_l%wfn(:nr,ispin,l_left,m_left)*wfn_r%wfn(:nr,ispin,l_right,m_right)*sd_rtab(1:nr))
              end if
              if (do_term(2)) then ! Dipole velocity
                if (l_left==l_right+1) wgt_vel = wgt_vel + sum(wfn_l%wfn(:nr,ispin,l_left,m_left)*grad_m(:nr,ispin))
                if (l_left==l_right-1) wgt_vel = wgt_vel + sum(wfn_l%wfn(:nr,ispin,l_left,m_left)*grad_p(:nr,ispin))
              end if
              if (do_term(3)) then ! Dipole acceleration
                wgt_acc = wgt_acc + sum(wfn_l%wfn(:nr,ispin,l_left,m_left)*wfn_r%wfn(:nr,ispin,l_right,m_right)*sd_dvdrtab(:nr))
              end if
            end do loop_spin
            wgt_ang       = angular_term(l_left,m_left,1_ik,m_op,l_right,m_right)
            dip_sph(m_op) = dip_sph(m_op) + electron_charge*wgt_ang*wgt_dip
            acc_sph(m_op) = acc_sph(m_op) - electron_charge*wgt_ang*wgt_acc/electron_mass
            !  Accumulating velocity terms is a bit ugly ...
            wgt_vel = wgt_vel * (0._rk,-1._rk) * electron_charge/electron_mass
            if ( (l_left==l_right+1) .and. (m_left==m_right  ) ) vel_sph( 0) = vel_sph( 0) + wgt_vel * vel_c(l_right+1, m_right  )
            if ( (l_left==l_right-1) .and. (m_left==m_right  ) ) vel_sph( 0) = vel_sph( 0) + wgt_vel * vel_c(l_right,   m_right  )
            if ( (l_left==l_right+1) .and. (m_left==m_right+1) ) vel_sph( 1) = vel_sph( 1) + wgt_vel * vel_d(l_right+1, m_right  )
            if ( (l_left==l_right-1) .and. (m_left==m_right+1) ) vel_sph( 1) = vel_sph( 1) - wgt_vel * vel_d(l_right,  -m_right-1)
            if ( (l_left==l_right+1) .and. (m_left==m_right-1) ) vel_sph(-1) = vel_sph(-1) + wgt_vel * vel_d(l_right+1,-m_right  )
            if ( (l_left==l_right-1) .and. (m_left==m_right-1) ) vel_sph(-1) = vel_sph(-1) - wgt_vel * vel_d(l_right,   m_right-1)
          end do loop_l_left
        end do loop_m_left
      end do loop_l_right
      !$omp end do nowait
    end do loop_m_right
    !$omp end parallel
    !
    !  Transform to Cartesian coordinates
    !
    call sph2cart(dip_sph,dipole(:,1))
    call cyc2cart(vel_sph,dipole(:,2))
    call sph2cart(acc_sph,dipole(:,3))
    !
    call nt_add(dipole)
    call TimerStop('WF dipole')
    !
    contains
    subroutine prepare_gradients(nr,wfn,lval,scr,tmp,grad_m,grad_p)
      integer(ik), intent(in)    :: nr         ! Radial extent of the non-zero part of the WF
      complex(rk), intent(in)    :: wfn(:)     ! Radial wavefunction
      integer(ik), intent(in)    :: lval       ! Angular momentum of the channel
      complex(rk), intent(inout) :: scr(:,:)
      complex(rk), intent(inout) :: tmp(:)
      complex(rk), intent(out)   :: grad_m(:)  !  (d/dr) Psi - ((l+1)/r)*psi
      complex(rk), intent(out)   :: grad_p(:)  !  (d/dr) Psi +     (l/r)*psi
      !
      if (lval==0) then
        call m3d_multiply(sd_d1n_l0,wfn,tmp)
        call m3d_solve(sd_m1n_l0,sd_m1nf_l0,sd_m1np_l0,tmp,grad_m,scr)
      else
        call m3d_multiply(sd_d1n_lx,wfn,tmp)
        call m3d_solve(sd_m1n_lx,sd_m1nf_lx,sd_m1np_lx,tmp,grad_m,scr)
      end if
      !
      grad_p(:nr) = grad_m(:nr) + (lval  )*sd_rminus(:nr)*wfn(:nr)
      grad_m(:nr) = grad_m(:nr) - (lval+1)*sd_rminus(:nr)*wfn(:nr)
    end subroutine prepare_gradients
    !
    real(rk) function angular_term(l1,m1,l2,m2,l3,m3)
      integer(ik), intent(in) :: l1, m1, l2, m2, l3, m3
      ! Brute-force expression using 3J symbols; Mathematica phase convenstions
      angular_term = ((-1)**m1) * sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/(4*pi)) &
                   * Math3J(l1,l2,l3,0_ik,0_ik,0_ik) * Math3J(l1,l2,l3,-m1,m2,m3)
    end function angular_term
    !
    real(rk) function vel_c(l,m) 
      integer(ik), intent(in) :: l, m
      real(rk) :: r_l, r_m
      !
      r_l   = l 
      r_m   = m
      vel_c = sqrt( (r_l-r_m)*(r_l+r_m)/( (2*r_l-1)*(2*r_l+1) ) )
    end function vel_c
    !
    real(rk) function vel_d(l,m) 
      integer(ik), intent(in) :: l, m
      real(rk) :: r_l, r_m
      !
      r_l   = l 
      r_m   = m
      vel_d = sqrt( (r_l+r_m)*(r_l+r_m+1)/( 2*(2*r_l-1)*(2*r_l+1) ) )
    end function vel_d
    !
    subroutine sph2cart(sph,cart)
      complex(rk), intent(in)  :: sph (-1:1) ! Spherical tensor
      complex(rk), intent(out) :: cart(3)    ! Equivalent Cartesian tensor
      !
      cart(1) =       sqrt((2*pi)/3) * (sph(-1)-sph(+1))
      cart(2) = (0._rk,1._rk)*sqrt((2*pi)/3) * (sph(-1)+sph(+1))
      cart(3) =       sqrt((4*pi)/3) *  sph( 0)
    end subroutine sph2cart
    !
    subroutine cyc2cart(sph,cart)
      complex(rk), intent(in)  :: sph (-1:1) ! Cyclic vector
      complex(rk), intent(out) :: cart(3)    ! Equivalent Cartesian vector
      !
      cart(1) =         (sph(-1)-sph(+1)) / sqrt(2._rk)
      cart(2) = (0._rk,1._rk) * (sph(-1)+sph(+1)) / sqrt(2._rk)
      cart(3) =          sph( 0)
    end subroutine cyc2cart
  end subroutine wt_dipole
  !
  !  Calculate the conversion factor for transformation of the right wavefunction
  !  into the real-space wavefunction.
  !
  !  The idea of the calculation is very simple: we know that the real-space
  !  wavefunction coincides with the right wavefunction up to the normalization
  !  factor. We therefore try to calculate the real-space norm of the right
  !  wavefunction as carefully as possible. Using the known norm of the left/right
  !  pair, this immediately gives us the scaling coefficient.
  !
  function wt_r2r_scale(wfn_l,wfn_r) result(scale)
    type(sd_wfn), intent(in)      :: wfn_l     ! Left wavefunction
    type(sd_wfn), intent(in)      :: wfn_r     ! Right wavefunction
    real(rk)                      :: scale     ! Normalization factor converting the right wavefunction 
                                               ! into the real-space wavefunction.
    !
    integer(ik) :: lv, mv, ispin, ipt, my_lmax
    complex(rk) :: grad(sd_nradial)   ! Gradient of the right wavefunction
    complex(rk) :: lap (sd_nradial)   ! Laplacian of the right wavefunction
    complex(rk) :: tmp (sd_nradial)   ! Temporaries for calculating laplacian and gradient
    complex(rk) :: scr (sd_nradial,m3d_sc_size)
    complex(rk) :: lr_norm            ! Wavefunction norm calculated from PsiL * PsiR
    real(rk)    :: r2_norm            ! Wavefunction norm calculated from |PsiR|**2
    real(rk)    :: rl, rh             ! Lower and upper bounds of the current volume element for the radial integral
    complex(rk) :: psir(0:2)          ! Power series for the right wavefunction
    real(rk)    :: rhos(0:4)          ! Power series for |rhor|^2
    !
    call TimerStart('WF R2R')
    !
    my_lmax = max(wfn_l%lmax,wfn_r%lmax)
    lr_norm = 0
    r2_norm = 0
    !$omp parallel default(none) &
    !$omp& reduction(+:lr_norm,r2_norm) &
    !$omp& shared(sd_mmin,sd_mmax,my_lmax,sd_nspin,sd_nradial) &
    !$omp& shared(sd_d1n_l0,sd_m1n_l0,sd_m1nf_l0,sd_m1np_l0) &
    !$omp& shared(sd_d1n_lx,sd_m1n_lx,sd_m1nf_lx,sd_m1np_lx) &
    !$omp& shared(sd_d2n_l0,sd_m2n_l0,sd_m2nf_l0,sd_m2np_l0) &
    !$omp& shared(sd_d2n_lx,sd_m2n_lx,sd_m2nf_lx,sd_m2np_lx) &
    !$omp& shared(wfn_l,wfn_r,sd_drtab,nts) &
    !$omp& private(mv,lv,ispin,tmp,scr,grad,lap,ipt,rl,rh,psir,rhos)
    lr_norm = 0
    r2_norm = 0
    scan_mv: do mv=sd_mmin,sd_mmax
      !$omp do schedule(guided,1)
      scan_lv: do lv=abs(mv),my_lmax
        if (nts%dist(lv)>0) cycle scan_lv ! Hook for shared-memory parallelization
        scan_spin: do ispin=1,sd_nspin
          !
          !  We'll need the derivatives of the right wavefunction to get the real-space
          !  integral as accurately as possible.
          !
          if (lv==0) then
            call m3d_multiply(sd_d1n_l0,wfn_r%wfn(:,ispin,lv,mv),tmp)
            call m3d_solve(sd_m1n_l0,sd_m1nf_l0,sd_m1np_l0,tmp,grad,scr)
            call m3d_multiply(sd_d2n_l0,wfn_r%wfn(:,ispin,lv,mv),tmp)
            call m3d_solve(sd_m2n_l0,sd_m2nf_l0,sd_m2np_l0,tmp,lap,scr)
          else
            call m3d_multiply(sd_d1n_lx,wfn_r%wfn(:,ispin,lv,mv),tmp)
            call m3d_solve(sd_m1n_lx,sd_m1nf_lx,sd_m1np_lx,tmp,grad,scr)
            call m3d_multiply(sd_d2n_lx,wfn_r%wfn(:,ispin,lv,mv),tmp)
            call m3d_solve(sd_m2n_lx,sd_m2nf_lx,sd_m2np_lx,tmp,lap,scr)
          end if
          !
          !  Our Hilbert-space norm is easy ...
          !
          lr_norm = lr_norm + sum(wfn_l%wfn(:,ispin,lv,mv)*wfn_r%wfn(:,ispin,lv,mv))
          !
          !  For the real-space norm, assume second-order Taylor expansion around the grid point
          !
          scan_pt: do ipt=1,sd_nradial
            !
            !  For the first point, start at the world origin (not necessarily zero); otherwise, we 
            !  start at the midpoint.
            !  For the last point, go to the end of the world; otherwise, to the midpoint
            !
            rl = -sd_drtab(ipt)
            rh =  sd_drtab(ipt+1)
            if (ipt/=1) rl = 0.5_rk * rl
            if (ipt/=sd_nradial) rh = 0.5_rk * rh
            !
            !  Construct power series for psir and |psir|^2 expanded around the grid point
            !
            psir(0) = wfn_r%wfn(ipt,ispin,lv,mv)
            psir(1) = grad(ipt)
            psir(2) = lap (ipt)
            !
            rhos(0) = abs(psir(0))**2
            rhos(1) = 2._rk * real(psir(0)*conjg(psir(1)),kind=rk)
            rhos(2) = abs(rhos(1))**2 + real(psir(0)*conjg(psir(2)),kind=rk)
            rhos(3) = real(psir(1)*conjg(psir(2)),kind=rk)
            rhos(4) = 0.25_rk * abs(psir(2))**2
            !
            !  Integrate
            !
            r2_norm = r2_norm + rhos(0) *                 (rh   -rl   ) &
                              + rhos(1) * (1._rk/2._rk) * (rh**2-rl**2) &
                              + rhos(2) * (1._rk/3._rk) * (rh**3-rl**3) &
                              + rhos(3) * (1._rk/4._rk) * (rh**4-rl**4) &
                              + rhos(4) * (1._rk/5._rk) * (rh**5-rl**5)
          end do scan_pt
        end do scan_spin
      end do scan_lv
      !$omp end do
    end do scan_mv
    !$omp end parallel
    !
    call nt_add(r2_norm)
    call nt_add(lr_norm)
    !
    ! write (out,"('WF R2R: (L,R) = ',g32.24,1x,g32.24)") lr_norm
    ! write (out,"('WF R2R: <R|R> = ',g32.24)") r2_norm
    if (abs(r2_norm)>spacing(abs(lr_norm))) then
      scale = sqrt(abs(lr_norm)/r2_norm)
    else
      !
      !  Right wavefunction vanishes; it does not matter what the scale is as long as it is not singular
      !
      scale = 1
    end if
    call TimerStop('WF R2R')
  end function wt_r2r_scale
end module wavefunction_tools
