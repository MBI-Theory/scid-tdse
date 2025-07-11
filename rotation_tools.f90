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
!  Spherical harmonics rotation.
!
!  The routines in this modules are shared-memory parallelized (OpenMP). They
!  are aware of a possible distributed-memory parallelization, but rely on the 
!  upstream routines to set up and finalize distributed-memory calculations.
!
module rotation_tools
  use accuracy
  use constants
  use math
  use spherical_data
  use tridiagonal_tools
  use node_tools
  use timer
  implicit none
  private
  public rt_max_rot, rt_blocksize, rt_sense_real, rt_nstep_max
  public rt_rotate
  public rt_rotate_bruteforce
  public rcsid_rotation_tools
  !
  character(len=clen), save :: rcsid_rotation_tools = "$Id: rotation_tools.f90,v 1.37 2024/04/23 14:32:48 ps Exp $"
  !
  real(rk), save    :: rt_max_rot   = 1e-3_rk  ! Maximum allowed rotation angle, in Radian
                                               ! The actual rotation threshold will be (rt_max_rot/(sd_lmax+1))

  integer(ik), save :: rt_blocksize = 16_ik    ! Max number of R points to process as a block. This is
                                               ! a performance optimization parameter; it does not change the results
  logical, save     :: rt_sense_real= .true.   ! If .true., try to detect real rotational propagator, and use
                                               ! the corresponding special case
  integer(ik), save :: rt_nstep_max =  100_ik  ! The number of rotation steps could potentially get very large;
                                               ! we need to avoid allocating excessive amounts of memory here.
  !
  contains
  !
  subroutine rt_prepare_rotation_steps(from,to,nstep,alp)
    real(xk), intent(in)                 :: from(:)   ! Initial (th,ph) in radian, in at least double precision
    real(xk), intent(in)                 :: to  (:)   ! Final   (th,ph)
    integer(ik), intent(out)             :: nstep     ! Number of subdivided steps
    real(xk), allocatable, intent(inout) :: alp (:,:) ! Final small-rotation angles
    !
    integer(ik) :: is, alloc
    real(xk)    :: delta(2)
    real(xk)    :: step_from(2), step_to(2), step_mid(2), step_delta(2)
    !
    if (size(from)/=2 .or. size(to)/=2) then
      stop 'rotation_tools%rt_prepare_rotation_steps - bad argument sizes'
    end if
    !
    if (allocated(alp)) deallocate(alp)
    !
    delta = to - from
    delta = delta - (2._xk*pi_xk*nint(delta/(2._xk*pi_xk)))
    !
    !  For runs involving nodes with different architectures, it is safer to choose master's result here.
    !  Because of the use of nint() intrinsic, the phase wrap may happen differently on different nodes!
    !
    call nt_broadcast(delta)
    if (all(abs(delta)<=spacing(1._xk))) then
      nstep = 0
      return
    end if
    !
    nstep = ceiling(sd_lmax*maxval(abs(delta),dim=1)/rt_max_rot)
    ! write (out,"('Need ',i0,' steps to achieve angular step of ',g24.13)") nstep, rt_max_rot/sd_lmax
    allocate (alp(3,nstep),stat=alloc)
    if (alloc/=0) then
      stop 'rotation_tools%rt_prepare_rotation_steps - allocation failed'
    end if
    break_steps: do is=1,nstep
      step_from  = from + (delta*(is-1))/nstep
      step_to    = from + (delta*(is-0))/nstep
      step_mid   = 0.5_xk * (step_from + step_to)
      step_delta = step_to - step_from
      !
      alp(1,is) =  step_delta(2) * sin(step_mid(1))
      alp(2,is) = -step_delta(1)
      alp(3,is) = -step_delta(2) * cos(step_mid(1))
    end do break_steps
  end subroutine rt_prepare_rotation_steps
  !
  !  rt_rotate can be applied to either left or right wavefunctions; the only
  !  difference is the sign of the time step and swapping of the upper/lower diagonals
  !  of the coupling matrix
  !
  subroutine rt_rotate(psi,from,to,left)
    type(sd_wfn), intent(inout) :: psi       ! Wavefunction we act upon
    real(xk), intent(in)        :: from(:)   ! Initial (th,ph) in radian
    real(xk), intent(in)        :: to  (:)   ! Final   (th,ph)
    logical, intent(in)         :: left      ! Set to .true. if rotating left wavefunction
    !
    !
    integer(ik)              :: lval, mval, nstep, istep, ispin
    integer(ik)              :: my_lmax, my_mmin, my_mmax
    integer(ik)              :: mmin, mmax            ! Smallesr and largest M for the current L
    integer(ik)              :: ir1, irn, irc
    integer(ik)              :: sstep, estep          ! First (starting) and last sub-steps
    integer(ik)              :: cstep                 ! Number of sub-steps in the current block
    integer(ik)              :: nst                   ! Smaller of nstep and rt_nstep_max
    integer(ik)              :: alloc
    complex(rk)              :: r_tmp (sd_mmin:sd_mmax,rt_blocksize)
    complex(rk)              :: r_tmp2(sd_mmin:sd_mmax,rt_blocksize)
    real(xk), allocatable    :: alp(:,:)              ! Table of rotation angles;
                                                      ! first index:  (ax,ay,az)
                                                      ! second index: rotation step
    complex(rk), allocatable :: prop_fwd   (:,:,:)    ! Forward part of the propagator
    complex(rk), allocatable :: prop_rev   (:,:,:)    ! Inverse part of the propagator
    complex(rk), allocatable :: prop_revf  (:,:,:)    ! Factored inverse
    real(rk), allocatable    :: prop_fwd_r (:,:,:)    ! Forward part of the propagator (special real case)
    real(rk), allocatable    :: prop_rev_r (:,:,:)    ! Inverse part of the propagator (special real case)
    real(rk), allocatable    :: prop_revf_r(:,:,:)    ! Factored inverse (special real case)
    logical, allocatable     :: prop_revp  (  :,:)    ! Pivot list for the inverse (both complex and real)
    complex(rk), allocatable :: scr(:,:)              ! Scratch for iterative refinement
    complex(rk), allocatable :: srm(:,:)              ! Dense rotation matrix, from small-angle expansion
    logical                  :: use_dense             ! Use dense rotation branch; it's faster
    real(rk)                 :: eps
    logical                  :: use_real              ! Propagator is real; use faster special case
    ! 
    call TimerStart('Small-angle rotate')
    !
    call rt_prepare_rotation_steps(from,to,nstep,alp)
    if (nstep==0) then
      ! Early exit; nothing to do
      call TimerStop('Small-angle rotate')
      return
    end if
    !
    my_lmax = psi%lmax
    my_mmin = max(-my_lmax,sd_mmin)
    my_mmax = min( my_lmax,sd_mmax)
    nst     = min(nstep,rt_nstep_max)
    !$omp parallel default(none) &
    !$omp& shared(my_mmin,my_mmax,my_lmax,sd_nradial,sd_nspin,nstep,alp) &
    !$omp& shared(rt_blocksize,psi,left,rt_sense_real,nst,nts) &
    !$omp& private(prop_fwd,prop_rev,prop_fwd_r,prop_rev_r) &
    !$omp& private(prop_revf,prop_revf_r,prop_revp) &
    !$omp& private(srm,alloc,lval,mmin,mmax) &
    !$omp& private(use_dense,irn,irc,r_tmp,r_tmp2,eps,use_real) &
    !$omp& private(scr,sstep,estep,cstep)
    !
    !  Parallelizing over L allows non-interleaved memory access from each
    !  thread. The alternative would be to parallelize over R, in which point
    !  we'll run a risk of false sharing and cache-line contention between
    !  the threads executing on separate CPUs.
    !
    !  We'll run in the -decreasing- order of L, since this keeps smaller
    !  work units at the end of the loop
    !
    allocate (prop_fwd  (my_mmin:my_mmax,3,nst), prop_fwd_r (my_mmin:my_mmax,          3,nst), &
              prop_rev  (my_mmin:my_mmax,3,nst), prop_revf  (my_mmin:my_mmax,m3d_dc_size,nst), &
              prop_rev_r(my_mmin:my_mmax,3,nst), prop_revf_r(my_mmin:my_mmax,m3d_dc_size,nst), &
              prop_revp (my_mmin:my_mmax,  nst), &
              scr(my_mmin:my_mmax,m3d_sc_size*rt_blocksize), &
              srm(my_mmin:my_mmax,my_mmin:my_mmax),stat=alloc)
    if (alloc/=0) then
      write (out,*) '     my_mmin = ',my_mmin
      write (out,*) '     my_mmax = ',my_mmax
      write (out,*) '       nstep = ',nstep
      write (out,*) '         nst = ',nst
      write (out,*) ' m3d_dc_size = ',m3d_dc_size
      write (out,*) ' m3d_sc_size = ',m3d_sc_size
      write (out,*) 'rt_blocksize = ',rt_blocksize
      call flush_wrapper(out)
      stop 'rotation_tools%rt_rotate - allocation failed'
    end if
    !
    ! Note that L=0 does not require rotation: it is an invariant
    !$omp do schedule(dynamic,1)
    process_l: do lval=my_lmax,1,-1
      if (nts%dist(lval)>0) cycle process_l ! Hook for distributed-memory parallel execution
      step_block: do sstep=1,nstep,nst
        estep = min(nstep,sstep+nst-1)
        cstep = estep - sstep + 1
        ! Range of allowable M values is from -L to L
        mmin = max(-lval,my_mmin)
        mmax = min( lval,my_mmax)
        call build_small_angle_propagators(lval,mmin,mmax,alp(:,sstep:estep),prop_fwd(mmin:mmax,:,:cstep), &
                      prop_rev(mmin:mmax,:,:cstep),prop_revf(mmin:mmax,:,:cstep),prop_revp(mmin:mmax,:cstep))
        !
        ! Would it be faster to collapse series of sparse propagators into a dense matrix?
        !   "sparse" version costs 26*(mmax-mmin+1)*nstep FLOP
        !   "dense" version costs 5*(mmax-mmin+1)**3 FLOP; it is also more regular and vectorizes better.
        ! The "dense" version carries an additional cost of 26*(mmax-mmin+1)**2*nstep/sd_nradial for
        ! computing the rotation matrices, which may become significant when sd_nradial is relatively
        ! small compared to the maximum angular momentum.
        !   So, if nstep is much greater than (mmax-mmin+1)**2/5 or so, we should be switching to the dense code
        ! "much" is the measure of the relative efficiency of the code implementing dense matrix multiply and 
        ! tri-diagonal matrices; it is probably OK to assume dense linear algebra is considerably more efficient.
        !
        use_dense = nstep>(mmax-mmin+1._rk)**2/200 + ((mmax-mmin+1._rk)*nstep)/sd_nradial
        use_real  = .false.  ! Only relevant if use_dense is not set
        if (use_dense) then
          call expand_small_angle_propagators(prop_fwd(mmin:mmax,:,:cstep), &
                    prop_rev(mmin:mmax,:,:cstep),prop_revf(mmin:mmax,:,:cstep),prop_revp(mmin:mmax,:cstep), &
                    srm(mmin:mmax,mmin:mmax))
        else ! .not.use_dense
          !  If the propagator is real, we can save a bit of effort here
          ! use_real = .false. ! Mover the assignment up, to suppress a spurious gfortran warning
          if (rt_sense_real) then
            eps = spacing(max(maxval(abs(prop_fwd(mmin:mmax,:,:cstep))), &
                              maxval(abs(prop_rev(mmin:mmax,:,:cstep)))))
            use_real = all(abs(aimag(prop_fwd(mmin:mmax,:,:cstep)))<eps) .and. &
                       all(abs(aimag(prop_rev(mmin:mmax,:,:cstep)))<eps)
            if (use_real) then
              prop_fwd_r (mmin:mmax,:,:cstep) = real(prop_fwd (mmin:mmax,:,:cstep),kind=rk)
              prop_rev_r (mmin:mmax,:,:cstep) = real(prop_rev (mmin:mmax,:,:cstep),kind=rk)
              prop_revf_r(mmin:mmax,:,:cstep) = real(prop_revf(mmin:mmax,:,:cstep),kind=rk)
            end if
          end if
        end if
        !
        !  Now comes the expensive part: the actual transformation, which must be done for each R point. 
        !
        process_radial: do ir1=1,sd_nradial,rt_blocksize
          irn = min(sd_nradial,ir1+rt_blocksize-1)
          irc = irn - ir1 + 1
          process_spin: do ispin=1,sd_nspin
            ! Fetch wavefunction to cache
            fetch_radial: do mval=mmin,mmax
              r_tmp(mval,:irc) = psi%wfn(ir1:irn,ispin,lval,mval)
            end do fetch_radial
            ! If we are dealing with the left wavefunction, we need to conjugate it here,
            ! then again after rotation is complete
            if (left) r_tmp(mmin:mmax,:irc) = conjg(r_tmp(mmin:mmax,:irc))
            ! Apply the propagators to each R point in turn
            if (use_dense) then
              ! Gfortran creates a (potentially large!) temporary array for the result of matmult
              ! here. Using r_tmp2 as a temporary buffer does not seem to make any difference -
              ! it seems to be a problem with the matmul intrinsic itself.
              r_tmp(mmin:mmax,1:irc) = matmul(srm(mmin:mmax,mmin:mmax),r_tmp(mmin:mmax,1:irc))
            else
              apply_rotations: do istep=1,cstep
                if (use_real) then
                  call m3d_multiply(prop_fwd_r (mmin:mmax,:,istep),r_tmp (mmin:mmax,1:irc),r_tmp2(mmin:mmax,1:irc))
                  call m3d_solve   (prop_rev_r (mmin:mmax,:,istep), &
                                    prop_revf_r(mmin:mmax,:,istep), &
                                    prop_revp  (mmin:mmax,  istep),   &
                                    r_tmp2(mmin:mmax,1:irc),r_tmp (mmin:mmax,1:irc),scr(mmin:mmax,1:m3d_sc_size*irc))
                else ! .not. use_real
                  call m3d_multiply(prop_fwd   (mmin:mmax,:,istep),r_tmp (mmin:mmax,1:irc),r_tmp2(mmin:mmax,1:irc))
                  call m3d_solve   (prop_rev   (mmin:mmax,:,istep),  &
                                    prop_revf  (mmin:mmax,:,istep),  &
                                    prop_revp  (mmin:mmax,  istep),  &
                                    r_tmp2(mmin:mmax,1:irc),r_tmp (mmin:mmax,1:irc),scr(mmin:mmax,1:m3d_sc_size*irc))
                end if
              end do apply_rotations
            end if
            ! Undo the conjugation for the left wavefunction
            if (left) r_tmp(mmin:mmax,:irc) = conjg(r_tmp(mmin:mmax,:irc))
            ! Store transformed wavefunction back
            store_radial: do mval=mmin,mmax
              psi%wfn(ir1:irn,ispin,lval,mval) = r_tmp(mval,:irc)
            end do store_radial
          end do process_spin
        end do process_radial
      end do step_block
    end do process_l
    !$omp end do
    deallocate (prop_fwd,prop_rev,prop_revf,prop_fwd_r,prop_rev_r,prop_revf_r,prop_revp,scr,srm)
    !$omp end parallel
    if (allocated(alp)) deallocate (alp)  ! Theoretically, there is no gurantee it was allocated
    call TimerStop('Small-angle rotate')
  end subroutine rt_rotate
  !
  subroutine build_small_angle_propagators(lval,mmin,mmax,alp,prop_fwd,prop_rev,prop_revf,prop_revp)
    integer(ik), intent(in)  :: lval             ! Angular momentum L
    integer(ik), intent(in)  :: mmin             ! Lower bound of the angular momentum projection
    integer(ik), intent(in)  :: mmax             ! Upper bound ..
    real(xk), intent(in)     :: alp(:,:)         ! Rotation steps
    complex(rk), intent(out) :: prop_fwd (:,:,:) ! Forward halfs of the rotation propagators
    complex(rk), intent(out) :: prop_rev (:,:,:) ! Reverse halfs of the rotation propagators
    complex(rk), intent(out) :: prop_revf(:,:,:) ! Factores prop_rev
    logical,     intent(out) :: prop_revp(  :,:) ! Pivot list for prop_rev
    !
    integer(ik) :: nstep, istep
    integer(ik) :: mval
    real(rk)    :: ax, ay, az
    real(rk)    :: d0  (mmin:mmax)   ! D0 rotation coefficients: coupling to same L
    real(rk)    :: dp  (mmin:mmax)   ! D+ rotation coefficients: coupling to higher L
    real(rk)    :: dm  (mmin:mmax)   ! D- rotation coefficients: coupling to lower L
    complex(rk) :: tmp (mmin:mmax,3)
    complex(rk) :: tmp2(mmin:mmax,3)
    !
    nstep = size(alp,dim=2)
    !
    !  First, form the L-dependent coupling coefficients
    !
    fill_d_tables: do mval=mmin,mmax
      d0(mval) = mval
      dp(mval) = 0.5_rk*sqrt(real(lval-mval,kind=rk)*real(lval+mval+1,kind=rk))
      dm(mval) = 0.5_rk*sqrt(real(lval+mval,kind=rk)*real(lval-mval+1,kind=rk))
    end do fill_d_tables
    !
    !  Now, form propagator matrices at each time step
    !
    build_propagators: do istep=1,nstep
      ax = -0.5_rk * real(alp(1,istep),kind=rk)
      ay = -0.5_rk * real(alp(2,istep),kind=rk)
      az = -0.5_rk * real(alp(3,istep),kind=rk)
      tmp(mmin:mmax  ,1) = d0(mmin  :mmax  ) * (-az          )
      tmp(mmin:mmax-1,2) = dp(mmin  :mmax-1) * (-ax +(0._rk,1._rk)*ay)
      tmp(mmin:mmax-1,3) = dm(mmin+1:mmax  ) * (-ax -(0._rk,1._rk)*ay)
      ! Construct direct and inverse propagators at this time step
      ! [1 - i H dt/2]
      prop_fwd(:,:,istep) = -(0._rk,1._rk)*tmp(:,:)
      prop_fwd(:,1,istep) = 1._rk + prop_fwd(:,1,istep)
      ! [1 + i H dt/2]^-1
      tmp2    (:,:)       = +(0._rk,1._rk)*tmp(:,:)
      tmp2    (:,1)       = 1._rk + tmp2(:,1)
      !
      prop_rev(:,:,istep) = tmp2
      call m3d_decompose(prop_rev (:,:,istep), &
                         prop_revf(:,:,istep), &
                         prop_revp(:,  istep))
    end do build_propagators
  end subroutine build_small_angle_propagators
  !
  subroutine expand_small_angle_propagators(prop_fwd,prop_rev,prop_revf,prop_revp,rm)
    complex(rk), intent(in)  :: prop_fwd (:,:,:) ! Forward halfs of the rotation propagators
    complex(rk), intent(in)  :: prop_rev (:,:,:) ! Reverse halfs of the rotation propagators
    complex(rk), intent(in)  :: prop_revf(:,:,:) ! Factored reverse
    logical,     intent(in)  :: prop_revp(  :,:) ! Pivot list for the factored reverse
    complex(rk), intent(out) :: rm(:,:)          ! Explicit, dense rotation matrix
    !
    integer(ik) :: nstep, istep
    integer(ik) :: sz, icol
    complex(rk) :: scr(size(rm,dim=1),m3d_sc_size*size(rm,dim=1))
    complex(rk) :: rhs(size(rm,dim=1),size(rm,dim=1))
    !
    if (any(ubound(prop_fwd)/=ubound(prop_rev)) .or. &
        size(rm,dim=1)/=size(rm,dim=2) .or. size(rm,dim=1)/=size(prop_fwd,dim=1)) then
      stop 'rotation_tools%expand_small_angle_propagators - bad argument sizes'
    end if
    !
    nstep = size(prop_fwd,dim=3)
    sz    = size(rm,dim=1)
    !
    rhs(:,:) = 0
    fill_columns: do icol=1,sz
      rhs(icol,icol) = 1._rk
    end do fill_columns
    apply_terms: do istep=1,nstep
      call m3d_multiply(prop_fwd(:,:,istep),rhs,rm)
      call m3d_solve(prop_rev(:,:,istep),prop_revf(:,:,istep),prop_revp(:,istep),rm,rhs,scr)
    end do apply_terms
    rm = rhs
  end subroutine expand_small_angle_propagators
  !
  subroutine rt_finite_rotation_matrix(from,to,mult,rm)
    real(xk), intent(in)     :: from(:)  ! Initial theta, phi in radian
    real(xk), intent(in)     :: to  (:)  ! Final theta, phi in radian
    integer(ik), intent(in)  :: mult     ! Multiplicity, 2*j+1 (half-integer j is OK)
    complex(rk), intent(out) :: rm(:,:)  ! Finite-rotation matrix, transforming spherical
                                         ! harmonics at [from] to spherical harmonics at [to]
    !
    complex(rk) :: rm_from(mult,mult), rm_to(mult,mult)
    integer(ik) :: im1, im2
    real(rk)    :: e_from(3), e_to(3)
    !
    !  Conceptually, we need two rotations:
    !  1) From the initial point [from] to the lab system
    !  2) From the lab system to the final point [to]
    !
    !  It looks like there is an error in L&L expression for the Wigner functions,
    !  so that we have to flip the sign of the beta Euler angle to get the correct
    !  rotation matrix. Oh well.
    !
    e_from(1) = -from(2) ; e_from(2) = -from(1) ; e_from(3) = 0._rk
    e_to  (1) = -to  (2) ; e_to  (2) = -to  (1) ; e_to  (3) = 0._rk
    call MathYJMRotationMatrix(euler_angles=e_from,mult=mult,mat=rm_from)
    call MathYJMRotationMatrix(euler_angles=e_to,  mult=mult,mat=rm_to)
    rm_to = conjg(rm_to)
    rm = matmul(rm_to,transpose(rm_from))
    !
    !  MathYJMRotationMatrix assumes harmonics use L&L III phase conventions.
    !  Our propagator uses Mathematica phase conventions, which differ by
    !  an extra factor (-I)**(L+2M). Let's fix up the matrix now.
    !
    !  The factor we need is:  (-I)**(L+2*M1)/(-I)**(L+2*M2) = (-1)**(M1-M2),
    !  so it's sufficient to change the sign of the rotation matrix whenever 
    !  m1 and m2 have different parities.
    !
    fix_m2: do im2=1,mult,2
      fix_m1_even: do im1=2,mult,2
        rm(im1,im2) = -rm(im1,im2)
      end do fix_m1_even
      if (im2+1>mult) cycle fix_m2
      fix_m1_odd: do im1=1,mult,2
        rm(im1,im2+1) = -rm(im1,im2+1)
      end do fix_m1_odd
    end do fix_m2
  end subroutine rt_finite_rotation_matrix
  !
  subroutine rt_rotate_bruteforce(psi,from,to,left)
    type(sd_wfn), intent(inout) :: psi       ! Wavefunction we act upon
    real(xk), intent(in)        :: from(:)   ! Initial (th,ph) in radian
    real(xk), intent(in)        :: to  (:)   ! Final   (th,ph)
    logical, intent(in)         :: left      ! Set to .true. if rotating left wavefunction
    !
    integer(ik) :: lval, mval, ispin, my_lmax
    integer(ik) :: mmin, mmax
    integer(ik) :: ir1, irn, irc, irx
    complex(rk) :: rm(-sd_lmax:sd_lmax,-sd_lmax:sd_lmax)
    complex(rk) :: r_tmp (sd_mmin:sd_mmax,rt_blocksize)
    complex(rk) :: r_tmp2(sd_mmin:sd_mmax,rt_blocksize)
    !
    if (size(from)/=2 .or. size(to)/=2) stop 'rotation_tools%rt_rotate_bruteforce - bad dimensions'
    !
    !  Detect a no-op
    !
    if (all(abs(from-to)<=spacing(1._xk))) return
    !
    my_lmax = psi%lmax
    call TimerStart('Rotate brute force')
    !$omp parallel default(none) &
    !$omp& shared(my_lmax,sd_mmin,sd_mmax,sd_nradial,sd_nspin,rt_blocksize) &
    !$omp& shared(from,to,psi,left,nts) &
    !$omp& private(lval,rm,mmin,mmax,ir1,irn,irc,irx,r_tmp,r_tmp2)
    !$omp do
    process_l: do lval=my_lmax,0,-1
      if (nts%dist(lval)>0) cycle process_l ! Hook for distributed-memory parallelization
      call rt_finite_rotation_matrix(from,to,2_ik*lval+1_ik,rm(-lval:lval,-lval:lval))
      mmin = max(sd_mmin,-lval)
      mmax = min(sd_mmax, lval)
      process_radial: do ir1=1,sd_nradial,rt_blocksize
        irn  = min(sd_nradial,ir1+rt_blocksize-1)
        irc  = irn - ir1 + 1
        process_spin: do ispin=1,sd_nspin
          ! Fetch wavefunction to cache
          fetch_radial: do mval=mmin,mmax
            r_tmp(mval,:irc) = psi%wfn(ir1:irn,ispin,lval,mval)
          end do fetch_radial
          ! If we are dealing with the left wavefunction, we need to conjugate it here,
          ! then again after rotation is complete
          if (left) r_tmp(mmin:mmax,:irc) = conjg(r_tmp(mmin:mmax,:irc))
          ! Apply the rotation to each R point in turn
          process_radial_cache: do irx=1,irc
            r_tmp2(mmin:mmax,irx) = matmul(rm(mmin:mmax,mmin:mmax),r_tmp(mmin:mmax,irx))
          end do process_radial_cache
          ! Undo the conjugation for the left wavefunction
          if (left) r_tmp2(mmin:mmax,:irc) = conjg(r_tmp2(mmin:mmax,:irc))
          ! Store transformed wavefunction back
          store_radial: do mval=mmin,mmax
            psi%wfn(ir1:irn,ispin,lval,mval) = r_tmp2(mval,:irc)
          end do store_radial
        end do process_spin
      end do process_radial
    end do process_l
    !$omp end do
    !$omp end parallel
    call TimerStop('Rotate brute force')
  end subroutine rt_rotate_bruteforce
  !
end module rotation_tools
