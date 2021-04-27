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
!  Very simple biconjugate gradient solver, cribbed directly from Numerical Recipes [sec. 2.7]
!  Linear system matrix is assumed to be strongly diagonally dominant, so that very
!  few iterations are expected to be needed, and no stability problems are envisaged.
!
!  This routine does not need to be fully robust: it is used as a speed optimization.
!  If it fails, the calling routine will fall back to a slow, but robust canned solver.
!  However, it is imperative to keep the number of failures small.
!
module bicg_tools
  use accuracy
  use timer
  implicit none
  private
  public bicg_solve
  public bicg_failtrace, bicg_maxiter, bicg_epsilon
  public bicg_failure_count
  public rcsid_bicg_tools
  !
  character(len=clen) :: rcsid_bicg_tools = "$Id: bicg_tools.f90,v 1.11 2021/04/26 15:44:44 ps Exp ps $"
  !
  logical, save      :: bicg_failtrace     = .true. ! Produce a verbose report when solution fails; 
                                                    ! May increase the runtime if vectors are traced
                                                    ! (see the commented-out section below)
  integer(ik), save  :: bicg_maxiter       = 8_ik   ! Max. number of iteratios; if bi-CG iterations fail
                                                    ! to converge after pt_bicg_maxiter restarts, each with
                                                    ! pt_bicg_maxiter iterations, we'll declare failure to
                                                    ! the calling routine
  real(rk), save     :: bicg_epsilon       = 0      ! Relative accuracy required from bi-CG solver. Zero value
                                                    ! will use smaller of 1e-12 and spacing(1._rk). Note that
                                                    ! bicg_epsilon values greater than 1e-12 will likely produce
                                                    ! incorrect results.
  !
  !  Statistics
  !
  integer(ik), save  :: bicg_failure_count = 0      ! Number of times bi-CG solver failed to find a solution
  !
  contains
  !
  subroutine bicg_solve(rhs,x,iblob,cblob,matvec,fail)
    complex(rk), intent(in)  :: rhs(:)   ! Right-hand side of the equations
    complex(rk), intent(out) :: x  (:)   ! Solution vector
    integer(ik), intent(in)  :: iblob(:) ! Opaque internal state, passed along to the matrix-vector multiplication
    complex(rk), intent(in)  :: cblob(:) ! Ditto, but complex
    interface 
      subroutine matvec(op,iblob,cblob,vec,av)
        use accuracy
        character(len=*), intent(in) :: op       ! 'N' or 'T'
        integer(ik), intent(in)      :: iblob(:) ! An opaque binary blob, passed along from the cg_solve() call
        complex(rk), intent(in)      :: cblob(:) ! Ditto, complex
        complex(rk), intent(in)      :: vec(:)   ! Vector
        complex(rk), intent(out)     :: av (:)   ! Matrix-vector product: av = a^op . vec
      end subroutine matvec
    end interface 
    logical, intent(out)     :: fail             ! bi-CG iterations failed to converge
    !
    real(rk)    :: eps                           ! Desired relative norm of the residue
    complex(rk) :: rr(size(x)), rl(size(x))      ! Right and left residue
    complex(rk) :: pr(size(x)), pl(size(x))      ! Right and left search directions
    complex(rk) :: tmp(size(x))                  ! Scratch for matrix-vector product
    complex(rk) :: alpha, beta, rnorm, rnorm2, denom
    real(rk)    :: rn                            ! Cartesian norm of the residue
    real(rk)    :: rhsn                          ! Cartesian norm of the R.H.S.
    integer(ik) :: iter, irest                   ! Iteration and restart count
    !
    ! Data buffers for tracing failures; can get fairly large!
    ! complex(rk) :: trace_rr(size(x),bicg_maxiter,bicg_maxiter)
    ! complex(rk) :: trace_rl(size(x),bicg_maxiter,bicg_maxiter)
    ! complex(rk) :: trace_pr(size(x),bicg_maxiter,bicg_maxiter)
    ! complex(rk) :: trace_pl(size(x),bicg_maxiter,bicg_maxiter)
    complex(rk) :: trace_alpha(     bicg_maxiter,bicg_maxiter)
    complex(rk) :: trace_beta (     bicg_maxiter,bicg_maxiter)
    real(rk)    :: trace_rn   (     bicg_maxiter,bicg_maxiter)
    integer(ik) :: trace_niter(                  bicg_maxiter)  ! Number of iteraction before restart
    integer(ik) :: trace_code (     bicg_maxiter,bicg_maxiter)  ! Progress of last iteration before restart
    !
    !  Sanity checking
    !
    if (size(rhs)/=size(x)) stop 'bicg_tools%bicg_solve - dimension mismatch'
    !
    !  Initialization; since linear system matrix is known to be heavily
    !  diagonally-dominant, use the RHS as the initial guess.
    !
    eps  = bicg_epsilon
    if (eps<=0) eps = min(1e-12_rk,spacing(1._rk))
    rhsn = sqrt(sum(real(abs(rhs)**2,kind=rk)))
    x    = rhs
    restarts: do irest=1,bicg_maxiter
      call matvec('N',iblob,cblob,x,tmp)
      rr = rhs - tmp
      !
      rl = conjg(rr) ! The choice here is arbitrary; choosing the conjugate let's
                     ! us start with the usual Cartesian norm of the residue
      pr = rr ; pl = rl
      rnorm2 = sum(rl*rr)
      !
      !  If we guessed the solution already, this may cause iterations to fail;
      !  check whether the residue is below the threshold, and leave if it is.
      !
      if (sqrt(abs(rnorm2))<=eps*rhsn) then
        fail = .false.
        return
      end if
      iterations: do iter=1,bicg_maxiter
        if (bicg_failtrace) then
          trace_niter(       irest) = iter
          trace_code (  iter,irest) = 0
        ! trace_rl   (:,iter,irest) = rl
        ! trace_rr   (:,iter,irest) = rr
        ! trace_pr   (:,iter,irest) = pr
        ! trace_pl   (:,iter,irest) = pl
        end if
        rnorm  = rnorm2
        call matvec('N',iblob,cblob,pr,tmp)
        denom  = sum(pl*tmp)
        if (abs(denom)<=spacing(abs(rnorm))) then
          exit iterations ! Iteration failed; try to restart
        end if
        alpha  = rnorm/denom
        if (bicg_failtrace) trace_code (iter,irest) = 1
        if (bicg_failtrace) trace_alpha(iter,irest) = alpha
        x      = x + alpha*pr
        rr     = rr - alpha*tmp
        rn     = sqrt(sum(real(abs(rr)**2,kind=rk)))
        if (bicg_failtrace) trace_rn(iter,irest) = rn
        if (rn<=eps*rhsn) then
          fail = .false.
          return
        end if
        !
        call matvec('T',iblob,cblob,pl,tmp)
        rl     = rl - alpha*tmp
        rnorm2 = sum(rl*rr)
        if (abs(rnorm)<=spacing(abs(rnorm2))) then
          exit iterations ! Iteration failed; try to restart
        end if
        beta   = rnorm2/rnorm
        if (bicg_failtrace) trace_code(iter,irest) = 2
        if (bicg_failtrace) trace_beta(iter,irest) = beta
        pr     = rr + beta*pr
        pl     = rl + beta*pl
      end do iterations
    end do restarts
    !$omp atomic
    bicg_failure_count = bicg_failure_count + 1
    fail = .true.
    if (bicg_failtrace) then
      !$omp critical
      write (out,"(/'Bi-CG iterations failed to converge. Initiating post-mortem dump'/)")
      write (out,"('eps = ',g24.13,' r.h.s. norm = ',g24.13,' convergence criterion = ',g24.13)") &
             eps, rhsn, eps*rhsn
      write (out,"()")
      write (out,"(1x,a3,1x,a3,2x,a24,1x,a24,2x,a24,2x,a24,1x,a24)") &
             'Rst', 'Itr', ' Re[alpha] ', ' Im[alpha] ', ' Residue norm ', ' Re[beta] ', ' Im[beta] ', &
             '---', '---', '-----------', '-----------', '--------------', '----------', '----------'
      ft_dump_restarts: do irest=1,bicg_maxiter
        ft_dump_iterations: do iter=1,trace_niter(irest)
          write (out,"(1x,i3,1x,i3)",advance='no') irest, iter
          if (trace_code(iter,irest)<=0) then
            write (out,"('  Alpha denominator vanished')")
            cycle ft_dump_iterations
          end if
          write (out,"(2x,g24.13,1x,g24.13,2x,g24.13)",advance='no') trace_alpha(iter,irest), trace_rn(iter,irest)
          if (trace_code(iter,irest)<=1) then
            write (out,"('  Beta denominator vanished')")
            cycle ft_dump_iterations
          end if
          write (out,"(2x,g24.13,1x,g24.13)") trace_beta(iter,irest)
        end do ft_dump_iterations
      end do ft_dump_restarts
      ! write (out,"(/t5,'Initial R.H.S.'/)")
      ! write (out,"((t10,5(g14.7,1x,g14.7,2x)))") rhs
      ! call ft_dump_vectors('Right residue',         trace_rr)
      ! call ft_dump_vectors('Right search direction',trace_pr)
      ! call ft_dump_vectors('Left residue',          trace_rl)
      ! call ft_dump_vectors('Left search direction', trace_pl)
      write (out,"(/'End of Bi-CG post-mortem')")
      call flush_wrapper(out)
      !$omp end critical
    end if
    !
    contains
    subroutine ft_dump_vectors(tag,vec)
      character(len=*), intent(in) :: tag
      complex(rk), intent(in)      :: vec(:,:,:) ! Trace of a vector to dump
      !
      integer(ik) :: irest, iter
      !
      write (out,"(/t5,a/)") trim(tag)
      ftdv_restarts: do irest=1,bicg_maxiter
        ftdv_iterations: do iter=1,trace_niter(irest)
          write (out,"(1x,i3,1x,i3,(t10,5(g14.7,1x,g14.7,2x)))") irest, iter, vec(:,iter,irest)
        end do ftdv_iterations 
      end do ftdv_restarts
      write (out,"()")
    end subroutine ft_dump_vectors
  end subroutine bicg_solve
  !
end module bicg_tools
