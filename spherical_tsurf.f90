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
!  Implementation of the "t-SURF" method for calculation of photoelectron
!  spectrum in spherical coordinates. We support multiple t-SURF calculations
!  running at the same time, provided that they use the same k grid and the
!  same matching-sphere radius. 
!
!  Additionally, we can also calculate the amplitudes of the Coulomb partial
!  waves present in the simulation volume at the end of the simulation. These
!  amplitudes can be converted to the photoelectron spectrum, using Rutherford's
!  asymptotic Coulomb plane waves.
!
!  This module is prallelized for shared-memory (OpenMP) execution. Parts of the
!  time step step and the final analysis/corrections are also parallelized for 
!  distributed-memory runs. Enough of the work is repeated on distributed-memory
!  nodes, so that scaling with the number of nodes is unlikely to be great.
!
module spherical_tsurf
  use accuracy
  use constants
  use coulomb_functions
  use math
  use node_tools
  use spherical_bessel
  use spherical_data
  use spherical_tsurf_data
  use timer
  use tridiagonal_tools
  use wavefunction_tools
  implicit none
  private
  public sts_initialize_global
  public sts_initialize_instance, sts_destroy_instance
  public sts_timestep
  public sts_atend_prepare, sts_atend_terms
  public sts_atend_direct
  public sts_report
  public rcsid_spherical_tsurf
  !
  character(len=clen), save :: rcsid_spherical_tsurf = "$Id: spherical_tsurf.f90,v 1.39 2024/02/13 14:22:14 ps Exp $"
  !
  !  Private data, not visible outside the module
  !
  logical, parameter          :: detail_timer      = .false.
  !
  real(rk), allocatable, save :: sts_ktab(:)      ! Grid of radial k values; 1:sts_kgrid_count
  real(rk), allocatable, save :: sts_dtab(:,:)    ! Grid of k directions (laboratory frame). 
                                                  ! First index:  1:3 (Cartesian components of the direction)
                                                  ! Second index: 1:sts_dgrid_count 
  real(rk), allocatable, save :: sts_btab(:,:)    ! Table of spherical Bessel function values
                                                  ! First index: -1:sd_lmax+1 (order of the Bessel function)
                                                  ! Second index: 1:sts_kgrid_count (sts_rmatch * radial grid position)
  !
  contains
  !
  subroutine sts_initialize_global(dt,amax)
    real(xk), intent(in) :: dt   ! Time step in TDSE propagation; t-SURF time step may be longer
    real(xk), intent(in) :: amax ! Estimated peak of the vector potential
    !
    integer(ik) :: total_kvec
    real(rk)    :: nyquist_pimax ! Maximum kinetic momentum supported by the time step
    !
    !  First, decide what kind of calculation was asked for
    !
    sts_active       = sts_volkov
    !
    select case (sts_atend_mode)
      case default
        write (out,"('iSURF: Final projection algorithm is not implemented: ',a)") trim(sts_atend_mode)
        stop 'spherical_tsurf%sts_initialize_global - Bad sts_atend_mode'
      case ('direct')
        sts_active_atend  = .false.
        sts_active_direct = sts_volkov_atend .or. sts_coulomb_atend
      case ('spectral')
        sts_active_atend  = sts_volkov_atend .or. sts_coulomb_atend
        sts_active_direct = .false.
    end select
    !
    if (.not.(sts_active.or.sts_active_atend.or.sts_active_direct)) return
    !
    call TimerStart('t-SURF: Global initialize')
    if (sts_verbose>=0) then
      write (out,"(/'Preparing for photoelectron spectrum calculation'/)")
      write (out,"( 'Project on Volkov states during the TDSE (t-SURF) = ',L6)") sts_volkov
      write (out,"( '          Project on Volkov states after the TDSE = ',L6)") sts_volkov_atend
      write (out,"( '         Project on Coulomb states after the TDSE = ',L6)") sts_coulomb_atend
      write (out,"( '          Algorithm used for the final projection = ',a/)") trim(sts_atend_mode)
    end if
    !
    if (.not.sd_capped) then
      write (out,"('Photoelectron spectrum calculation is meaningless in the absence of the absorbing boundary')")
      stop 'spherical_tsurf%sts_initialize_global - t-SURF requires an absorbing boundary'
    end if
    !
    call initialize_global_rmatch
    !
    call initialize_global_kgrid
    !
    call initialize_global_dgrid
    !
    nyquist_pimax = sqrt(2._rk*pi/(real(dt,kind=rk)*sts_step_frequency))
    if (sts_volkov .and. ((maxval(sts_ktab)+real(amax,kind=rk)>0.1_rk*nyquist_pimax) .or. sts_verbose>=0)) then
      write (out,"(/'       Effective time step in t-SURF: ',g24.12,' jiffies')") dt*sts_step_frequency
      write (out,"( '  Nyquist limit for kinetic momentum: ',g24.12,' Bohr/jiffy')") nyquist_pimax
      write (out,"( 'Maximum canonical momentum requested: ',g24.12,' Bohr/jiffy')") maxval(sts_ktab)
      write (out,"( '  Estimated maximum kinetic momentum: ',g24.12,' Bohr/jiffy'/)") maxval(sts_ktab) + real(amax,kind=rk)
    end if
    !
    if (sts_verbose>=0) then
      total_kvec = sts_dgrid_count * sts_kgrid_count
      write (out,"()")
      write (out,"(' There are altogether ',i10,' distinct K vectors')") total_kvec
      write (out,"('       Global t-SURF data requires ',f0.3,' Mbytes')") &
           (size(sts_ktab)+size(sts_dtab)+size(sts_btab))*rk_bytes()/1024._rk**2
      write (out,"(' Per-instance t-SURF data requires ',f0.3,' Mbytes')") &
           5._rk*total_kvec*rk_bytes()/1024._rk**2
      write (out,"()")
    end if
    !
    call TimerStop('t-SURF: Global initialize')
  end subroutine sts_initialize_global
  !
  subroutine initialize_global_rmatch
    !
    !  Figure out radial position of the matching sphere
    !
    if (sts_rmatch<=0) then
      sts_rmatch = sd_rtab(sd_cap_start) - abs(sts_rmatch)
    end if
    sts_ipt_match = minloc(abs(sd_rtab-sts_rmatch),dim=1)
    !
    !  If processes happen to execute on nodes with different architecture, there is no guarantee that
    !  all nodes will yield the same result for sts_ipt_match. Which will be A VERY BAD THING. It is
    !  much safer to broadcast whatever the master computes.
    !
    call nt_broadcast(sts_ipt_match)
    if (sts_verbose>=0) then
      write (out,"('      matching at radial grid point ',i8)") sts_ipt_match
      write (out,"('  matching distance from the origin ',f24.15,' Bohr (requested ',f14.5,')')") &
             sd_rtab(sts_ipt_match), sts_rmatch
    end if
    sts_rmatch    = sd_rtab(sts_ipt_match)
    if (sts_ipt_match<=1 .or. sts_ipt_match>=sd_nradial) then
      stop 'spherical_tsurf%sts_initialize_global - can''t match at the grid edge'
    end if
    if (sts_ipt_match>=sd_cap_start) then
      write (out,"(/'WARNING: t-SURF matching point is inside the absorbing boundary'/)")
    end if
  end subroutine initialize_global_rmatch
  !
  subroutine initialize_global_kgrid
    integer(ik) :: alloc, ikp, lv, iu_temp
    !
    !  Initialize radial grid for k
    !
    if (sts_verbose>=0) then
      write (out,"(/'  initializing ',a,' radial grid for the final momentum k')") trim(sts_kgrid)
    end if
    select case (sts_kgrid)
      case default
        write (out,"('spherical_tsurf: Radial k grid ',a,' is not recognized')") trim(sts_kgrid)
        stop 'spherical_tsurf%sts_initialize_global - bad std_kgrid'
      case ('uniform')
        if (sts_verbose>=0) then
          write (out,"('  radial k grid starts at ',f14.6,' Bohr/jiffy (exclusive)')") sts_kgrid_min 
          write (out,"('    radial k grid ends at ',f14.6,' Bohr/jiffy (inclusive)')") sts_kgrid_max 
          write (out,"('   radial k grid contains ',i8,' points')") sts_kgrid_count
        end if
        if (sts_kgrid_min<0 .or. sts_kgrid_min>=sts_kgrid_max .or. sts_kgrid_count<=0) then
          stop 'spherical_tsurf%sts_initialize_global - bad std_kgrid_*'
        end if
      case ('read')
        if (sts_verbose>=0) then
          write (out,"('   radial k grid will be read from ',a)") trim(sts_kgrid_file)
          write (out,"('   radial k grid contains ',i8,' points')") sts_kgrid_count
        end if
    end select
    allocate (sts_ktab(sts_kgrid_count),sts_btab(-1:sd_lmax+1,sts_kgrid_count),stat=alloc)
    if (alloc/=0) then
      stop 'spherical_tsurf%sts_initialize_global - allocate sts_ktab and sts_btab'
    end if
    select case (sts_kgrid)
      case default ; stop 'spherical_tsurf%sts_initialize_global - kgrid logic error'
      case ('uniform')
        ! sts_kgrid_min is not a part of the grid: we can't handle k=0
        fill_uniform_kgrid: do ikp=1,sts_kgrid_count
          sts_ktab(ikp) = sts_kgrid_min + (ikp*(sts_kgrid_max-sts_kgrid_min))/sts_kgrid_count
        end do fill_uniform_kgrid
      case ('read')
        open (newunit=iu_temp,form='formatted',action='read',position='rewind',status='old',file=trim(sts_kgrid_file))
        read (iu_temp,*) sts_ktab
        close (iu_temp)
        if (any(sts_ktab<=0)) then
          write (out,"(/'Encountered a non-positive k-vector magnitude'/)")
          write (out,*) sts_ktab
          stop 'spherical_tsurf%sts_initialize_global - k-vector magnitudes must be positive'
        end if
    end select
    !$omp parallel do default(none) shared(sts_kgrid_count,sd_lmax,sts_rmatch,sts_ktab,sts_btab) private(ikp)
    fill_bessel_table: do ikp=1,sts_kgrid_count
      call besselj_table(sd_lmax+1,sts_rmatch*sts_ktab(ikp),sts_btab(:,ikp))
    end do fill_bessel_table
    !$omp end parallel do
    !
    if (sts_verbose>=1) then
      write (out,"(/t5,'Radial k-vector grid'/)")
      write (out,"((1x,a5,1x,a20,t27,2x,a5,1x,a20,1x,a3))") &
             ' IKM ', ' |K|, Bohr/jiffy ', ' L ', '  SphericalBessel(L,Rmatch*|K|) ', '...', &
             '-----', '-----------------', '---', '--------------------------------' 
      report_radial_kgrid: do ikp=1,sts_kgrid_count
        write (out,"(1x,i5,1x,g20.12,(t27,5(2x,i5,1x,g20.12)))") &
               ikp, sts_ktab(ikp), (lv, sts_btab(lv,ikp), lv=-1,sd_lmax+1)
      end do report_radial_kgrid
      write (out,"()")
    end if
  end subroutine initialize_global_kgrid
  !
  subroutine initialize_global_dgrid
    integer(ik) :: alloc, ith, iph, idp, iu_temp
    real(rk)    :: th, ph
    !
    !  Fill direction grid for the final K vector
    !
    if (sts_verbose>=0) then
      write (out,"(/'  initializing ',a,' direction grid for the final momentum k')") trim(sts_dgrid)
    end if
    select case (sts_dgrid)
      case default
        write (out,"('spherical_tsurf: Angular k grid ',a,' is not recognized')") trim(sts_dgrid)
        stop 'spherical_tsurf%sts_initialize_global - bad std_dgrid'
      case ('product')
        if (sts_verbose>=0) then
          write (out,"('   using ',i8,' points along spherical theta angle')") sts_dgrid_ntheta
          write (out,"('   using ',i8,' points along spherical phi angle')") sts_dgrid_nphi
        end if
        if (sts_dgrid_ntheta<=0 .or. sts_dgrid_nphi<=0) then
          stop 'spherical_tsurf%sts_initialize_global - bad std_dgrid_*'
        end if
        sts_dgrid_count = sts_dgrid_ntheta * sts_dgrid_nphi
      case ('read')
        if (sts_verbose>=0) then
          write (out,"('   reading direction grid from ',a)") trim(sts_dgrid_file)
        end if
    end select
    if (sts_verbose>=0) then
      write (out,"('   direction grid containts ',i8,' points overall')") sts_dgrid_count
    end if
    if (sts_dgrid_count<=0) then
      stop 'spherical_tsurf%sts_initialize_global - sts_dgrid_count must be positive'
    end if
    allocate (sts_dtab(3,sts_dgrid_count),stat=alloc)
    if (alloc/=0) then
      stop 'spherical_tsurf%sts_initialize_global - allocate sts_dtab'
    end if
    select case (sts_dgrid)
      case default ; stop 'spherical_tsurf%sts_initialize_global - dgrid logic error'
      case ('product')
        !$omp parallel do default(none) shared(sts_dgrid_ntheta,sts_dgrid_nphi,sts_dtab) &
        !$omp&            private(ith,th,iph,ph)
        fill_direction_theta: do ith=1,sts_dgrid_ntheta
          th = (ith-0.5_rk)*pi/sts_dgrid_ntheta
          fill_direction_phi: do iph=1,sts_dgrid_nphi
            ph = (iph-0.5_rk)*2._rk*pi/sts_dgrid_nphi
            ! gfortran really generates bad code for array assignments ...
            sts_dtab(1,iph+(ith-1)*sts_dgrid_nphi) = sin(th)*cos(ph)
            sts_dtab(2,iph+(ith-1)*sts_dgrid_nphi) = sin(th)*sin(ph)
            sts_dtab(3,iph+(ith-1)*sts_dgrid_nphi) = cos(th)
          end do fill_direction_phi
        end do fill_direction_theta
        !$omp end parallel do
      case ('read')
        open (newunit=iu_temp,form='formatted',action='read',position='rewind',status='old',file=trim(sts_dgrid_file))
        read (iu_temp,*) sts_dtab
        close (iu_temp)
        normalize_direction_grid: do idp=1,sts_dgrid_count
          sts_dtab(:,idp) = sts_dtab(:,idp) / sqrt(sum(sts_dtab(:,idp)**2))
        end do normalize_direction_grid
    end select
    !
    if (sts_verbose>=1) then
      write (out,"(/t5,'Direction k-vector grid'/)")
      write (out,"((1x,a6,1x,a15,1x,a15,1x,3(1x,a17)))") &
             ' IPT ', ' Theta, degree ', ' Phi, degree ', ' Kx/|K| ', ' Ky/|K| ', ' Kz/|K| ', &
             '-----', '---------------', '-------------', '--------', '--------', '--------'
      report_direction_grid: do idp=1,sts_dgrid_count
        th = acos(sts_dtab(3,idp))*180._rk/pi
        ph = atan2(sts_dtab(2,idp),sts_dtab(1,idp))*180._rk/pi
        write (out,"(1x,i6,1x,f15.10,1x,f15.10,1x,3(1x,g17.10))") idp, th, ph, sts_dtab(:,idp)
      end do report_direction_grid
      write (out,"()")
    end if
  end subroutine initialize_global_dgrid
  !
  subroutine sts_initialize_instance(sdt,wfn_l,wfn_r,volkov_opendx,volkov_table, &
                                     coulomb_waves,coulomb_table,coulomb_opendx)
    type(sts_data), intent(inout)          :: sdt
    type(sd_wfn), intent(in)               :: wfn_l          ! Left wavefunction
    type(sd_wfn), intent(in)               :: wfn_r          ! Right wavefunction
    character(len=*), intent(in), optional :: volkov_opendx  ! See corresponding sts_* global variables
    character(len=*), intent(in), optional :: volkov_table   ! 
    character(len=*), intent(in), optional :: coulomb_waves  ! 
    character(len=*), intent(in), optional :: coulomb_table  ! 
    character(len=*), intent(in), optional :: coulomb_opendx ! 
    !
    integer(ik) :: alloc1, alloc2
    !
    alloc1 = 0 ; alloc2 = 0
    !
    if (sts_volkov .or. sts_volkov_atend) then
      allocate (sdt%vphase   (sts_kgrid_count,sts_dgrid_count), &
                sdt%pref     (sts_kgrid_count,sts_dgrid_count), &
                sdt%amplitude(sts_kgrid_count,sts_dgrid_count),stat=alloc1)
    end if
    !
    if (sts_coulomb_atend) then
      allocate (sdt%fcgr (2,sts_kgrid_count,0:sd_lmax), &
                sdt%fcamp(sts_kgrid_count,0:sd_lmax,sd_mmin:sd_mmax), &
                sdt%coulamp(sts_kgrid_count,sts_dgrid_count),stat=alloc2)
    end if
    !
    if (alloc1/=0 .or. alloc2/=0) then
      stop 'spherical_tsurf%sts_initialize_instance - allocation failure'
    end if
    !
    if (sts_volkov .or. sts_volkov_atend) then
      sdt%vphase    = 0
      sdt%amplitude = 0
    end if
    !
    if (sts_r2r_scale<=0._rk) then
      sdt%wfn_scale = wt_r2r_scale(wfn_l,wfn_r)
    else
      sdt%wfn_scale = sts_r2r_scale
    end if
    if (sts_verbose>=0) then
      write (out,"(/' Right to real-space wavefunction conversion factor: ',g32.24/)") sdt%wfn_scale
    end if
    sdt%volkov_opendx  = sts_volkov_opendx  ; if (present(volkov_opendx))  sdt%volkov_opendx  = volkov_opendx
    sdt%volkov_table   = sts_volkov_table   ; if (present(volkov_table ))  sdt%volkov_table   = volkov_table
    sdt%coulomb_waves  = sts_coulomb_waves  ; if (present(coulomb_waves))  sdt%coulomb_waves  = coulomb_waves
    sdt%coulomb_table  = sts_coulomb_table  ; if (present(coulomb_table))  sdt%coulomb_table  = coulomb_table
    sdt%coulomb_opendx = sts_coulomb_opendx ; if (present(coulomb_opendx)) sdt%coulomb_opendx = coulomb_opendx
  end subroutine sts_initialize_instance
  !
  subroutine sts_destroy_instance(sdt)
    type(sts_data), intent(inout) :: sdt
    !
    if (.not.(sts_active .or. sts_active_atend .or. sts_active_direct)) return
    !
    if (allocated(sdt%vphase)) deallocate (sdt%vphase)
    if (allocated(sdt%pref)) deallocate (sdt%pref)
    if (allocated(sdt%amplitude)) deallocate (sdt%amplitude)
    if (allocated(sdt%fcgr)) deallocate (sdt%fcgr)
    if (allocated(sdt%fcamp)) deallocate (sdt%fcamp)
    if (allocated(sdt%coulamp)) deallocate (sdt%coulamp)
  end subroutine sts_destroy_instance
  !
  !  The interface below is a bit eclectic, but it uses the quantities which are readily
  !  available during the time propagation - so I guess there is no harm done ...
  !
  subroutine sts_timestep(wfn_l,wfn_r,sdt,dt1,dt2,vp,efield)
    type(sd_wfn), intent(in)      :: wfn_l     ! Left wavefunction; not used
    type(sd_wfn), intent(in)      :: wfn_r     ! Right wavefunction
    type(sts_data), intent(inout) :: sdt       ! t-SURF data to update
    real(xk), intent(in)          :: dt1       ! Time step from the beginning of the interval to the time
                                               ! where vector-potential and electric field are specified
    real(xk), intent(in)          :: dt2       ! Time step from the VP/E expansion point to the end of the 
                                               ! interval. The time steps do not need to coincide with the
                                               ! propagation time step, or to remain constant throughout the
                                               ! simulation.
    real(xk), intent(in)          :: vp(:)     ! Vector-potential in laboratory spherical coordinates: (a,theta,phi)
                                               ! Note that (a) does not have to be positive!
    real(xk), intent(in)          :: efield(:) ! Electric field (aka -d A/d t) in laboratory Cartesian 
                                               ! coordinates: (ex,ey,ez)
    !
    integer(ik)              :: ikm, ikd, lv, mv, my_lmax, my_mmin, my_mmax
    integer(ik)              :: alloc     
    real(rk)                 :: th, ph       ! Theta and phi angles for the vector-potential in the lab c.s.
    real(rk)                 :: euler(3)     ! Euler angles for MathRotationMatrix - stupic gfortran ...
    real(rk)                 :: az           ! Vector-potential along local Z
    real(rk)                 :: krm(3,3)     ! Rotation matrix for the lab->local coordinate system transformation
    real(rk)                 :: kd_loc(3)    ! K-vector direction in the local coordinate system
    complex(rk), allocatable :: psilm(:,:,:) ! Wavefunction and radial derivative at the matching point; global
                                             ! First index: 1 = wf; 2 = derivative; second index: L; third index: M
    complex(rk), allocatable :: ylm(:,:)     ! Spherical harmonics for the k-vector direction; per-thread
                                             ! First index: M; Second index: L
    complex(rk), allocatable :: sylm  (:,:)  ! 1 = (-I)**(L+1)  L/(4L+2)      Sum Y_{LM} Psi_{LM}(R)
                                             ! 2 = (-I)**(L+1)  -(L+1)/(4L+2) Sum Y_{LM} Psi_{LM}(R)
                                             ! 3 = (-I)**(L+1)  e Az          Sum C_{LM} Y_{L-1,M} Psi_{LM}(R)
                                             ! 4 = (-I)**(L+1)  -e Az         Sum C_{L+1,M} Y_{L+1,M} Y_{LM} Psi_{LM}(R)
                                             ! 5 = (-I)**(L+1)  0.5           Sum Y_{LM} d Psi_{LM}(R)/d R
                                             ! Second index: L, 0 to sd_lmax
    complex(rk)              :: accu
    complex(rk)              :: cs    
    !
    if (.not.sts_active) return
    !
    if (size(vp)/=3) stop 'spherical_tsurf%sts_timestep - bad vp parameter'
    !
    call TimerStart('t-SURF: Time step')
    my_lmax = wfn_r%lmax ! We never use the left wavefunction here
    my_mmin = max(-my_lmax,sd_mmin)
    my_mmax = min( my_lmax,sd_mmax)
    !
    !  Step 1: Calculate wavefunction and its radial derivative at the matching point for each channel
    !
    allocate (psilm(2,0:sd_lmax,sd_mmin:sd_mmax),stat=alloc)
    if (alloc/=0) then
      stop 'spherical_tsurf%sts_timestep - error allocating global t-SURF buffers'
    end if
    call fill_wavefunction_at_matching_point(wfn_r,sdt%wfn_scale,psilm)
    !
    !  Step 2: Update Volkov phase for all k-vectors we care about
    !
    call update_volkov_phases(sdt,dt1,dt2,vp,efield)
    !
    !  Step 3: Prepare rotation matrix for the lab to local coordinate transformation
    !
    az = real(vp(1),kind=rk)
    th = real(vp(2),kind=rk)
    ph = real(vp(3),kind=rk)
    euler(1) = ph ; euler(2) = th ; euler(3) = 0._rk
    call MathRotationMatrix(euler,krm)
    if (sts_verbose>=2) then
      write (out,"('       Time step to expansion point = ',g25.17)") dt1
      write (out,"('          Time step to interval end = ',g25.17)") dt2
      write (out,"('Vector-potential along local Z axis = ',g25.17)") az
      write (out,"('     Lab theta for the local Z axis = ',g25.17)") th
      write (out,"('       Lab phi for the local Z axis = ',g25.17)") ph
      write (out,"('             Electric field (X,Y,Z) = ',3g25.17)") efield
      write (out,"('Lab to local rotation matrix:')")
      write (out,"((5x,3(1x,f25.17)))") transpose(krm)
      write (out,"()")
    end if
    !
    !  Steps 4 through 6 are repeated for each direction of the K grid
    !
    !$omp parallel default(none) &
    !$omp& shared(my_lmax,my_mmin,my_mmax,sts_dgrid_count,krm,psilm,az,sts_kgrid_count) &
    !$omp& shared(sts_ktab,sts_dtab,sts_btab,sdt,sts_rmatch,sts_verbose,nts) &
    !$omp& private(ylm,sylm,ikd,kd_loc,ikm,accu,lv,cs,alloc)
    allocate (ylm(my_mmin:my_mmax,0:my_lmax+1),sylm(5,0:my_lmax),stat=alloc)
    if (alloc/=0) then
      stop 'spherical_tsurf%sts_timestep - error allocating per-thread t-SURF buffers'
    end if
    !$omp do schedule(dynamic,1)
    kvec_direction: do ikd=1,sts_dgrid_count
      !
      !  Step 4: Precalculate spherical harmonics for this direction.
      !
      kd_loc(:) = matmul(krm,sts_dtab(:,ikd))
      if (sts_verbose>=2) then
        write (out,"(/'Direction ',i6,' Lab: ',3(1x,f12.9),' Local: ',3(1x,f12.9))") ikd, sts_dtab(:,ikd), kd_loc
      end if
      if (detail_timer) call TimerStart('t-SURF: Time step: YLM')
      call MathAllYLM2(my_lmax+1,my_mmin,my_mmax,kd_loc,ylm(my_mmin:my_mmax,0:my_lmax+1),phase='Arfken')
      if (detail_timer) call TimerStop('t-SURF: Time step: YLM')
      if (sts_verbose>=3) then
        write (out,"(/t5,'Spherical harmonics'/)")
        write (out,"((1x,a5,1x,a6,3(1x,a12),1x,2(1x,a20)))") &
               ' L ', ' M ', ' Kx ', ' Ky ', ' Kz ', ' Re[YLM] ', ' Im[YLM] ', &
               '---', '---', '----', '----', '----', '---------', '---------'
        print_ylm_l: do lv=0,my_lmax+1
          print_ylm_m: do mv=my_mmin,my_mmax
            write (out,"(1x,i5,1x,i6,3(1x,f12.9),1x,2(1x,f20.15))") lv, mv, kd_loc, ylm(mv,lv)
          end do print_ylm_m
        end do print_ylm_l
        write (out,"()")
      end if
      !
      !  Step 5: Contract spherical harmonics with wavefunction and derivatives
      !
      if (detail_timer) call TimerStart('t-SURF: Time step: Partial-M')
      call calculate_partial_sums_over_m(my_lmax,my_mmin,my_mmax,ylm,psilm,az,sylm)
      if (detail_timer) call TimerStop('t-SURF: Time step: Partial-M')
      !
      !  Step 6: Contract with all possible values of K magnitude in this direction
      !
      if (detail_timer) call TimerStart('t-SURF: Time step: |K|-accumulate')
      kvec_magnitude: do ikm=1,sts_kgrid_count
        accu = 0
        sum_lv: do lv=0,my_lmax
          if (nts%dist(lv)>0) cycle sum_lv ! Hook for distributed-memory parallization
          accu = accu + sts_ktab(ikm) * sts_btab(lv+1,ikm) * sylm(1,lv)
          accu = accu + sts_ktab(ikm) * sts_btab(lv-1,ikm) * sylm(2,lv)
          accu = accu +                 sts_btab(lv-1,ikm) * sylm(3,lv)
          accu = accu +                 sts_btab(lv+1,ikm) * sylm(4,lv)
          accu = accu +                 sts_btab(lv  ,ikm) * sylm(5,lv)
        end do sum_lv
        !
        !  Add the prefactor, and we are done here. Note that sdt%pref
        !  already includes the time step.
        !
        cs = sqrt(2._rk/pi) * sdt%pref(ikm,ikd) * sts_rmatch / electron_mass
        sdt%amplitude(ikm,ikd) = sdt%amplitude(ikm,ikd) + cs * accu
      end do kvec_magnitude
      if (detail_timer) call TimerStop('t-SURF: Time step: |K|-accumulate')
    end do kvec_direction
    !$omp end do
    !
    !  Free per-thread buffers
    !
    deallocate (ylm,sylm)
    !$omp end parallel
    !
    !  Free global buffers
    !
    deallocate (psilm)
    !
    call TimerStop('t-SURF: Time step')
  end subroutine sts_timestep
  !
  subroutine fill_wavefunction_at_matching_point(wfn_r,scale,psilm)
    type(sd_wfn), intent(in) :: wfn_r                              ! Right wavefunction
    real(rk), intent(in)     :: scale                              ! Scaling factor to convert right wavefunction -> real-space w.f.
    complex(rk), intent(out) :: psilm(2,0:sd_lmax,sd_mmin:sd_mmax) ! Wavefunction (1) and derivative (2) at matching point
    !
    complex(rk) :: tmp(sd_nradial), grad(sd_nradial) ! Per-thread temporary arrays
    complex(rk) :: scr(sd_nradial,m3d_sc_size)       ! ditto, scratch for iterative linear solver
    integer(ik) :: lv, mv, mmin, mmax, my_lmax
    !
    call TimerStart('t-SURF: Time step: Matching point')
    if (wfn_r%nradial<sts_ipt_match) then
      !
      !  Wavefunction has not reached the sensing surface; the result is zero
      !
      psilm = 0
      call TimerStop('t-SURF: Time step: Matching point')
      return
    end if
    !
    if (sd_nspin/=1) then
      stop 'spherical_tsurf%fill_wavefunction_at_matching_point - sd_nspin/=1 not implemented'
    end if
    if (sts_ipt_match<=1 .or. sts_ipt_match>=sd_nradial) then
      stop 'spherical_tsurf%fill_wavefunction_at_matching_point - bad sts_ipt_match'
    end if
    my_lmax = wfn_r%lmax
    !$omp parallel do default(none) &
    !$omp& shared(my_lmax,sd_mmin,sd_mmax,nts) &
    !$omp& shared(sd_d1n_l0,sd_m1n_l0,sd_m1nf_l0,sd_m1np_l0) &
    !$omp& shared(sd_d1n_lx,sd_m1n_lx,sd_m1nf_lx,sd_m1np_lx) &
    !$omp& shared(wfn_r,psilm,sts_ipt_match,scale) &
    !$omp& private(lv,mmin,mmax,mv,tmp,grad,scr)
    wavefunction_l: do lv=0,my_lmax
      if (nts%dist(lv)>0) cycle wavefunction_l ! Hook for distributed-memory parallization
      mmin = max(-lv,sd_mmin)
      mmax = min( lv,sd_mmax)
      wavefunction_m: do mv=mmin,mmax
        !
        !  Evaluate gradient; even though we need just one point, we have to evaluate the gradient
        !  at all radial points.
        !
        if (lv==0) then
          call m3d_multiply(sd_d1n_l0,wfn_r%wfn(:,1,lv,mv),tmp)
          call m3d_solve(sd_m1n_l0,sd_m1nf_l0,sd_m1np_l0,tmp,grad,scr)
        else
          call m3d_multiply(sd_d1n_lx,wfn_r%wfn(:,1,lv,mv),tmp)
          call m3d_solve(sd_m1n_lx,sd_m1nf_lx,sd_m1np_lx,tmp,grad,scr)
        end if
        psilm(1,lv,mv) = scale * wfn_r%wfn(sts_ipt_match,1,lv,mv)
        psilm(2,lv,mv) = scale * grad     (sts_ipt_match)
      end do wavefunction_m
    end do wavefunction_l
    !$omp end parallel do
    if (sts_verbose>=2) then
      write (out,"(/t5,'Wavefunction and gradient at the matching point'/)")
      write (out,"((1x,a5,1x,a6,2(2x,a24,1x,a24)))") &
             ' L ', ' M ', ' Re[Psi] ', ' Im[Psi] ', ' Re[d Psi/d R] ', ' Im[d Psi/d R] ', &
             '---', '---', '---------', '---------', '---------------', '---------------'
      print_wavefunction_l: do lv=0,my_lmax
        if (nts%dist(lv)>0) cycle print_wavefunction_l ! Hook for distributed-memory parallization
        mmin = max(-lv,sd_mmin)
        mmax = min( lv,sd_mmax)
        print_wavefunction_m: do mv=mmin,mmax
          write (out,"(1x,i5,1x,i6,2(2x,g24.16,1x,g24.16))") lv, mv, psilm(:,lv,mv)
        end do print_wavefunction_m
      end do print_wavefunction_l
      write (out,"()")
    end if
    call TimerStop('t-SURF: Time step: Matching point')
  end subroutine fill_wavefunction_at_matching_point
  !
  !  update_volkov_phases is parallelized with OpenMP, but does not make sense to
  !  parallelize across nodes (we need all phases everywhere!). It therefore limits
  !  the achievable parallel scaling for distributed-memory runs.
  !
  subroutine update_volkov_phases(sdt,dt1,dt2,vp,efield)
    type(sts_data), intent(inout) :: sdt        ! t-SURF data to updata
    real(xk), intent(in)          :: dt1        ! Time step from beginning of the interval to A(t) reference point
    real(xk), intent(in)          :: dt2        ! Time step from A(t) reference point to the end of the interval
    real(xk), intent(in)          :: vp(:)      ! Vector-potential in laboratory spherical coordinates: (a,theta,phi)
    real(xk), intent(in)          :: efield(:)  ! Electric field (=-dA/dt) in laboratory Cartesian coordinates
    !
    real(rk)    :: vp_xyz(3), kvec(3) ! Cartesian representation of the vector-potential and final electron momentum
    real(rk)    :: pvec(3)            ! Cartesian final momentum at the A(t) reference point
    real(rk)    :: ekin, rdt1, rdt2, f2
    real(rk)    :: pdotf              ! Scalar product of kinetic momentum and electric field
    complex(rk) :: eff_dt
    integer(ik) :: ik, id
    !
    call TimerStart('t-SURF: Time step: Volkov')
    vp_xyz(1) = real(vp(1)*sin(vp(2))*cos(vp(3)),kind=rk)
    vp_xyz(2) = real(vp(1)*sin(vp(2))*sin(vp(3)),kind=rk)
    vp_xyz(3) = real(vp(1)*cos(vp(2)),kind=rk)
    rdt1      = real(dt1,kind=rk)
    rdt2      = real(dt2,kind=rk)
    f2        = real(sum(efield**2),kind=rk)
    !
    if (sts_verbose>=2) then
      write (out,"(/t6,'Updating Volkov phases')")
      write (out,"('      Time steps = ',1x,g24.16,1x,g24.16,' jiffies')") dt1, dt2
      write (out,"('Vector-potential = ',3(1x,g24.16))") vp_xyz
      write (out,"('  Electric field = ',3(1x,g24.16))") efield
      write (out,"('             F^2 = ',1(1x,g24.16))") f2
    end if
    !$omp parallel do default(none) &
    !$omp& shared(sts_dgrid_count,sts_kgrid_count,sts_ktab,sts_dtab) &
    !$omp& shared(sts_verbose,vp_xyz,efield,rdt1,rdt2,sdt,f2) &
    !$omp& private(id,ik,kvec,pvec,ekin,pdotf,eff_dt) &
    !$omp& schedule(dynamic,1)
    kvec_direction: do id=1,sts_dgrid_count
      kvec_magnitude: do ik=1,sts_kgrid_count
        kvec(:) = sts_ktab(ik)*sts_dtab(:,id)
        pvec(:) = kvec(:) - electron_charge*vp_xyz(:)  ! Kinetic momentum at the expansion point
        ekin    = 0.5_rk * sum(pvec**2) / electron_mass
        pdotf   = sum(pvec*efield)
        ! 
        !  Calculate integral:  exp(I*phase) Int( Exp(I*Ekin*tau) d tau ) from 0 to dt
        !
        eff_dt = effective_dt(ekin,pdotf,rdt1,rdt2)
        !
        if (sts_verbose>=3) then
          write (out,"('Direction ',i6,' magnitude ',i6)") id, ik
          write (out,"('     kvec = ',3(1x,f24.16))") kvec
          write (out,"('     pvec = ',3(1x,f24.16))") pvec
          write (out,"('     ekin = ',1(1x,f24.16))") ekin
          write (out,"('    pdotf = ',1(1x,f24.16))") pdotf
          write (out,"('   eff_dt = ',2(1x,f24.16))") eff_dt
          write (out,"('old phase = ',1(1x,f24.16))") sdt%vphase(ik,id)
        end if
        !
        ! For some reason, Intel Fortran can't produce good code for the line below
        ! sdt%pref(ik,id) = exp((0._rk,1._rk)*sdt%vphase(ik,id)) * eff_dt
        sdt%pref(ik,id) = eff_dt * cmplx(cos(sdt%vphase(ik,id)),sin(sdt%vphase(ik,id)),kind=rk)
        !
        !  Advance the phase to the end of the time interval. We are being a bit inconsistent here,
        !  and integrate the phase to the second order in time, even though we only included the
        !  first order when doing the prefactor.
        !
        sdt%vphase(ik,id) = sdt%vphase(ik,id) + ekin * (rdt1+rdt2) &
                          + 0.5_rk * (electron_charge/electron_mass) * pdotf * (rdt2**2-rdt1**2) &
                          + (electron_charge**2/(6._rk*electron_mass)) * f2 * (rdt2**3 + rdt1**3)
        if (sdt%vphase(ik,id)>=2*pi) sdt%vphase(ik,id) = sdt%vphase(ik,id) - 2*pi
        if (sdt%vphase(ik,id)< 0   ) sdt%vphase(ik,id) = sdt%vphase(ik,id) + 2*pi
        !
        if (sts_verbose>=3) then
          write (out,"('new phase = ',1(1x,f24.16))") sdt%vphase(ik,id)
          write (out,"('prefactor = ',2(1x,f24.16))") sdt%pref  (ik,id)
        end if
      end do kvec_magnitude
    end do kvec_direction
    !$omp end parallel do
    if (sts_verbose>=3) then
      write (out,"(/t5,'Volkov phases and prefactors')")
      write (out,"((1x,a4,1x,a6,3(1x,a24),2x,a24,2x,a24,2x,a24))") &
             ' IK ', ' ID ', ' Kx ', ' Ky ', ' Kz ', ' Phase ', ' Re[pref] ', ' Im[pref] ', &
             '----', '----', '----', '----', '----', '-------', '----------', '----------'
      print_kvec_magnitude: do ik=1,sts_kgrid_count
        print_kvec_direction: do id=1,sts_dgrid_count
          kvec(:) = sts_ktab(ik)*sts_dtab(:,id)
          write (out,"(1x,i4,1x,i6,3(1x,g24.16),2x,g24.16,1x,2(g24.16))") id, ik, kvec, sdt%vphase(ik,id), sdt%pref(ik,id)
        end do print_kvec_direction
      end do print_kvec_magnitude
      write (out,"()")
    end if
    call TimerStop('t-SURF: Time step: Volkov')
    !
    contains
    !
    !  Calculate: 
    !
    !    (I ekin)**-1 (exp(I ekin (dt1+dt2))-1) [ 1 - (I e/(6 m)) pdotf (2*t1**2-t2**2+t1*t2) ]
    !
    !  taking care not to loose accuracy for small dt
    !
    function effective_dt(ekin,pdotf,dt1,dt2) result(edt)
      real(rk), intent(in) :: ekin, pdotf, dt1, dt2
      complex(rk)          :: edt
      !
      real(rk) :: dt
      real(rk) :: arg
      !
      dt  = dt1 + dt2
      arg = ekin*dt
      if (abs(arg**4)<=spacing(120._rk)) then
        !
        !  Order-3 Series expansion is exact, so use it.
        !
        edt = dt * cmplx(1._rk-arg**2/6._rk,arg/2._rk-arg**3/24._rk,kind=rk)
      else
        !
        !  Argument is too large; go with the exponent. We may lose at most 4 digits in double precision here
        !
        !  For some reason, Intel Fortran has trouble producing efficient code for the commented line below
        ! edt = -(0._rk,1._rk)*(exp((0._rk,1._rk)*arg)-1._rk)/ekin
        edt = cmplx(sin(arg),1._rk-cos(arg),kind=rk)/ekin
      end if
      edt = edt * cmplx(1._rk,(-electron_charge/(6._rk*electron_mass))*pdotf*(2*dt1**2-dt2**2+dt1*dt2),kind=rk)
    end function effective_dt
  end subroutine update_volkov_phases
  !
  subroutine calculate_partial_sums_over_m(my_lmax,my_mmin,my_mmax,ylm,psilm,az,sylm)
    integer(ik), intent(in)  :: my_lmax                            ! Largest angular momentum to consider
    integer(ik), intent(in)  :: my_mmin, my_mmax                   ! Smallest and largest angular momentum projections to consider
    complex(rk), intent(in)  :: ylm(my_mmin:my_mmax,0:my_lmax+1)   ! Spherical harmonics for the K direction
    complex(rk), intent(in)  :: psilm(2,0:sd_lmax,sd_mmin:sd_mmax) ! Per-channel wavefunction and derivative at matching point
    real(rk), intent(in)     :: az                                 ! Vector-potential along local Z axis
    complex(rk), intent(out) :: sylm(5,0:my_lmax)                  ! Partial sums over the M quantum number; see comment in sts_timestep
    !
    integer(ik) :: lv, mv, mmin, mmax
    real(rk)    :: rlv
    complex(rk) :: accu
    !
    loop_l: do lv=0,my_lmax
      if (nts%dist(lv)>0) cycle loop_l ! Hook for distributed-memory parallelization
      !
      rlv  = real(lv,kind=rk)
      mmin = max(sd_mmin,-lv)
      mmax = min(sd_mmax, lv)
      !
      !  1 = L/(4L+2)      Sum Y_{LM} Psi_{LM}(R)
      !  2 = -(L+1)/(4L+2) Sum Y_{LM} Psi_{LM}(R)
      !  5 = 0.5           Sum Y_{LM} d Psi_{LM}(R)/d R
      !
      accu = sum(ylm(mmin:mmax,lv)*psilm(1,lv,mmin:mmax))
      sylm(1,lv) =          (rlv/(4._rk*rlv+2._rk)) * accu
      sylm(2,lv) = -((rlv+1._rk)/(4._rk*rlv+2._rk)) * accu
      sylm(5,lv) = 0.5_rk * sum(ylm(mmin:mmax,lv)*psilm(2,lv,mmin:mmax))
      !
      !  3 = e Az Sum C_{LM} Y_{L-1,M} Psi_{LM}(R)
      !
      accu = 0
      minus_1: do mv=max(mmin,-(lv-1)),min(mmax,lv-1)
        accu = accu + clm(lv,mv)*ylm(mv,lv-1)*psilm(1,lv,mv) ! Spurious gfortran warning here
      end do minus_1
      sylm(3,lv) = electron_charge * az * accu
      !
      !  4 = -e Az Sum C_{L+1,M} Y_{L+1,M} Y_{LM} Psi_{LM}(R)
      !
      accu = 0
      plus_1: do mv=mmin,mmax
        accu = accu + clm(lv+1,mv)*ylm(mv,lv+1)*psilm(1,lv,mv)
      end do plus_1
      sylm(4,lv) = -electron_charge * az * accu
      !
      !  Add the overall factor of I**(L+1)
      !
      ! For some strange reason, Intel Fortran has trouble generating good code for 
      ! the commented line below
      ! sylm(:,lv) = ((0,-1)**(lv+1)) * sylm(:,lv)
      sylm(:,lv) = mipow(modulo(lv+1,4)) * sylm(:,lv)
    end do loop_l
    !
    if (sts_verbose>=2) then
      write (out,"(/t5,'Partial contractions with spherical harmonics over M'/)")
      write (out,"((1x,a5,5(1x,2(1x,a15))))") &
             ' L ', ' YL (1) ', ' ... ', ' YL (2) ', ' ... ', ' YL-1 ', ' ... ', ' YL+1 ', ' ... ', ' YL G ', ' ... ', &
             '---', '--------', '-----', '--------', '-----', '------', '-----', '------', '-----', '------', '-----'
      print_contractions_1: do lv=0,sd_lmax
        if (nts%dist(lv)>0) cycle print_contractions_1 ! Hook for distributed-memory parallelization
        write (out,"(1x,i5,5(1x,2(1x,g15.7)))") lv, sylm(:,lv)
      end do print_contractions_1
      write (out,"()")
    end if
    !
    contains 
    real(rk) function clm(l,m)
      integer(ik), intent(in) :: l, m
      !
      real(rk) :: lr, mr
      !
      lr  = real(l,kind=rk)
      mr  = real(m,kind=rk)
      clm = sqrt((lr**2-mr**2)/(4._rk*lr**2-1._rk))
    end function clm
  end subroutine calculate_partial_sums_over_m
  !
  !  Precompute common factors needed in evaluating field-free contributions to t-SURF
  !
  !  Volkov phases. It is trivial, but we need a lot of trig operations, and this
  !  works better when we call them all at once.
  !
  subroutine sts_atend_prepare_volkov(sdt)
    type(sts_data), intent(inout) :: sdt    ! t-SURF data structure. 
    !
    integer(ik) :: ikm, ikd
    !
    if (.not.sts_volkov_atend) return
    !
    call TimerStart('t-SURF: At end: Init: Volkov')
    !
    !  Update phases of the Volkov solutions until the end of the simulation
    !
    !$omp parallel do default(none) private(ikm,ikd) shared(sdt,sts_dgrid_count,sts_kgrid_count)
    kvec_direction: do ikd=1,sts_dgrid_count
      kvec_magnitude: do ikm=1,sts_kgrid_count
        ! Intel Fortran generates suboptimal code for the commented line below!
        ! sdt%pref(ikm,ikd) = exp((0._rk,1._rk)*sdt%vphase(ikm,ikd))
        sdt%pref(ikm,ikd) = cmplx(cos(sdt%vphase(ikm,ikd)),sin(sdt%vphase(ikm,ikd)),kind=rk)
      end do kvec_magnitude
    end do kvec_direction
    !$omp end parallel do
    call TimerStop('t-SURF: At end: Init: Volkov')
  end subroutine sts_atend_prepare_volkov
  !
  !  Precompute Coulomb scattering waves in each channel at the scattering point.
  !  The calculation is much more efficient of all L,M values are done at once.
  !
  subroutine sts_atend_prepare_coulomb(sdt)
    type(sts_data), intent(inout) :: sdt    ! t-SURF data structure. 
    !
    integer(ik) :: ikm
    real(rk)    :: c_f(0:sd_lmax), c_g(0:sd_lmax)
    !
    if (.not.sts_coulomb_atend) return
    !
    !  Precalculate Coulomb scattering waves at the matching point
    !
    call TimerStart('t-SURF: At end: Init: Coulomb')
    !
    !$omp parallel do default(none) private(ikm,c_f,c_g) &
    !$omp& shared(sts_kgrid_count,sd_lmax,sts_asymptotic_q,sts_ktab,sts_rmatch,sdt)
    coulomb_kvec_magnitude: do ikm=1,sts_kgrid_count
      !
      !  coulombF expects Z>0 to mean attraction (ie it assumes that electron_charge must be -1).
      !  Hence we use Z = -sts_asymptotic_q*electron_charge
      !
      call coulombF(0_ik,sd_lmax,-sts_asymptotic_q*electron_charge,sts_ktab(ikm),sts_rmatch,c_f,c_g)
      sdt%fcgr(1,ikm,0:sd_lmax) = c_f(:)
      sdt%fcgr(2,ikm,0:sd_lmax) = c_g(:)
    end do coulomb_kvec_magnitude
    !$omp end parallel do
    call TimerStop('t-SURF: At end: Init: Coulomb')
  end subroutine sts_atend_prepare_coulomb
  !
  subroutine sts_atend_prepare(sdt,wfn_l,wfn_r)
    type(sts_data), intent(inout) :: sdt    ! t-SURF data structure. 
    type(sd_wfn), intent(in)      :: wfn_l  ! Left wavefunction
    type(sd_wfn), intent(in)      :: wfn_r  ! Right wavefunction
    !
    if (.not.sts_active_atend) return
    !
    call TimerStart('t-SURF: At end: Init')
    !
    call sts_atend_prepare_volkov(sdt) 
    !
    call sts_atend_prepare_coulomb(sdt)
    !
    call TimerStop('t-SURF: At end: Init')
  end subroutine sts_atend_prepare
  !
  !  Solve auxiliary linear systems for each K-vector magnitude.
  !  Remember the wavefunctions and their derivatives.
  ! 
  !  We need to compute an auxiliary quantity Z:
  !
  !    Z = (H - eps)^{-1} Psi
  !
  !  As usual, this is done by solving an auxiliary tridiagonal system of equations:
  !
  !    [-0.5 D + M (V-eps) ] Z = M Psi
  !
  subroutine sts_atend_direct_auxiliary(sdt,wfn_r,aux_rd)
    type(sts_data), intent(inout) :: sdt             ! t-SURF data structure. 
    type(sd_wfn), intent(in)      :: wfn_r           ! Right wavefunction
    complex(rk), intent(out)      :: aux_rd(2,0:sd_lmax,sd_mmin:sd_mmax,sts_kgrid_count)
                                                     ! Solutions of auxiliary linear systems and derivatives
                                                     ! at the matching point, for each L, M, and K magnitude
                                                     ! Index 1: 1 = function, 2 = derivative
                                                     ! Index 2: 0:sd_lmax
                                                     ! Index 3: sd_mmin:sd_mmax
                                                     ! Index 4: 1 .. sts_kgrid_count
    integer(ik)                   :: ikm, lv, mv                  ! Counters: |K|, L, and M
    integer(ik)                   :: mmin, mmax
    real(rk)                      :: eps                          ! Continuum energy = k^2/(2*m)
    complex(rk)                   :: vminuseps(sd_nradial)        ! potential - continuum energy
    complex(rk)                   :: leq (sd_nradial,3)           ! Auxiliary linear system: -0.5*D + M 
    complex(rk)                   :: leqf(sd_nradial,m3d_dc_size) ! Factorized linear system
    logical                       :: leqp(sd_nradial)             ! Pivot list
    complex(rk)                   :: scr (sd_nradial,m3d_sc_size) ! Scratch for linear solver
    complex(rk)                   :: rhs (sd_nradial)             ! Right-hand side for linear solver
    complex(rk)                   :: aux (sd_nradial)             ! Z
    complex(rk)                   :: auxg(sd_nradial)             ! Gradient of Z
    !
    !$omp parallel do collapse(2) &
    !$omp& default(none) &
    !$omp& shared(sts_kgrid_count,sts_ktab,sd_lmax,sd_pottab,sd_capped,sd_captab) &
    !$omp& shared(sd_m2n_l0,sd_m2n_lx,sd_d2n_l0,sd_d2n_lx,sd_mmin,sd_mmax) &
    !$omp& shared(sd_d1n_l0,sd_m1n_l0,sd_m1nf_l0,sd_m1np_l0) &
    !$omp& shared(sd_d1n_lx,sd_m1n_lx,sd_m1nf_lx,sd_m1np_lx) &
    !$omp& shared(aux_rd,sd_cap_start,wfn_r,sts_ipt_match,sdt,nts) &
    !$omp& private(ikm,eps,lv,vminuseps,leq,leqf,leqp,mv,rhs,aux,scr,auxg)
    kvec_magnitude: do ikm=1,sts_kgrid_count
      lval_loop: do lv=0,sd_lmax
        if (nts%dist(lv)>0) cycle lval_loop  ! Hook for distributed-memory parallelization
        !
        eps = (0.5_rk/electron_mass) * sts_ktab(ikm)**2
        ! Multiplicative part of the (H-eps)
        vminuseps(:) = sd_pottab(:,1,lv) - eps
        if (sd_capped) then
          vminuseps(sd_cap_start:) = vminuseps(sd_cap_start:) + sd_captab(sd_cap_start:)
        end if
        !
        !  Build auxiliary linear system and factorize it. We need different matrices
        !  for different angular channels.
        !
        if (lv==0) then
          call m3d_right_scale(sd_m2n_l0,vminuseps,leq)
          leq = leq - (0.5_rk/electron_mass) * sd_d2n_l0
        else
          call m3d_right_scale(sd_m2n_lx,vminuseps,leq)
          leq = leq - (0.5_rk/electron_mass) * sd_d2n_lx
        end if
        call m3d_decompose(leq,leqf,leqp)
        !
        !  The linear system is the same for all M channels; only the right-hand side changes
        !
        mval_loop: do mv=sd_mmin,sd_mmax
          if (lv==0) then
            call m3d_multiply(sd_m2n_l0,wfn_r%wfn(:,1,lv,mv),rhs)         ! rhs = M . psi
            call m3d_solve(leq,leqf,leqp,rhs,aux,scr)                     ! Z
            call m3d_multiply(sd_d1n_l0,aux,rhs)                          ! Evaluating the gradient
            call m3d_solve(sd_m1n_l0,sd_m1nf_l0,sd_m1np_l0,rhs,auxg,scr)
          else
            call m3d_multiply(sd_m2n_lx,wfn_r%wfn(:,1,lv,mv),rhs)         ! rhs = M . psi
            call m3d_solve(leq,leqf,leqp,rhs,aux,scr)                     ! Z
            call m3d_multiply(sd_d1n_lx,aux,rhs)                          ! Evaluating the gradient
            call m3d_solve(sd_m1n_lx,sd_m1nf_lx,sd_m1np_lx,rhs,auxg,scr)
          end if
          ! Remember auxiliary and gradient at the matching point; we do not need the rest.
          ! Don't forget to apply the scale factor to the right wavefunction!
          aux_rd(1,lv,mv,ikm) = sdt%wfn_scale * aux (sts_ipt_match)
          aux_rd(2,lv,mv,ikm) = sdt%wfn_scale * auxg(sts_ipt_match)
        end do mval_loop
      end do lval_loop
    end do kvec_magnitude
    !$omp end parallel do
    if (sts_verbose>=2) then
      write (out,"(/t5,'Auxiliary solution and gradient at the matching point'/)")
      write (out,"((1x,a6,1x,a5,1x,a6,1x,a18, 2(2x,a24,1x,a24)))") &
             ' IK ', ' L ', ' M ', '  |K|  ', '  Re[Z]  ', '  Im[Z]  ', '  Re[d Z/d R]  ', '  Im[d Z/d R]  ', &
             '----', '---', '---', '-------', '---------', '---------', '---------------', '---------------'
      print_aux_k: do ikm=1,sts_kgrid_count
        print_aux_l: do lv=0,sd_lmax
          if (nts%dist(lv)>0) cycle print_aux_l ! Hook for distributed-memory parallelization
          !
          mmin = max(-lv,sd_mmin)
          mmax = min( lv,sd_mmax)
          print_aux_m: do mv=mmin,mmax
            write (out,"(1x,i6,1x,i5,1x,i6,g16.8,1x,2(2x,g24.16,1x,g24.16))") ikm, lv, mv, sts_ktab(ikm), aux_rd(:,lv,mv,ikm)
          end do print_aux_m
        end do print_aux_l
      end do print_aux_k
      write (out,"()")
      call flush_wrapper(out)
    end if
  end subroutine sts_atend_direct_auxiliary
  !
  subroutine sts_atend_direct_coulomb(sdt,aux_rd)
    type(sts_data), intent(inout) :: sdt             ! t-SURF data structure.
    complex(rk), intent(in)       :: aux_rd(2,0:sd_lmax,sd_mmin:sd_mmax,sts_kgrid_count)
                                                     ! Solutions of auxiliary linear systems and derivatives
                                                     ! at the matching point, for each L, M, and K magnitude
                                                     ! Index 1: 1 = function, 2 = derivative
                                                     ! Index 2: 0:sd_lmax
                                                     ! Index 3: sd_mmin:sd_mmax
                                                     ! Index 4: 1 .. sts_kgrid_count
    !
    integer(ik)              :: ikm
    integer(ik)              :: mmin, mmax
    integer(ik)              :: lv, mv
    complex(rk)              :: gc
    !
    call TimerStart('t-SURF: Direct: Coulomb')
    !
    !$omp parallel do collapse(2) default(none) &
    !$omp& private(ikm,lv,mmin,mmax,mv,gc) &
    !$omp& shared(sts_kgrid_count,sd_lmax,sd_mmin,sd_mmax,sdt,aux_rd,nts)
    coulomb_kvec_magnitude: do ikm=1,sts_kgrid_count
      lval_loop: do lv=0,sd_lmax
        sdt%fcamp(ikm,lv,:) = 0
        if (nts%dist(lv)>0) cycle lval_loop ! Hook for distributed-memory parallization
        mmin = max(sd_mmin,-lv)
        mmax = min(sd_mmax, lv)
        mval_loop: do mv=mmin,mmax
          !
          !  We assume that vector potential is zero!
          !
          !  gc = (hbar**2)/(2 Me) * [ (d chi*/d r) Z - chi* (d Z/d r) ]
          !
          !  where chi is the Coulomb scattering solution (real, so conjugation is immaterial)
          !  and Z is the auxiliary linear system solution
          !
          gc = sdt%fcgr(2,ikm,lv)*aux_rd(1,lv,mv,ikm) - sdt%fcgr(1,ikm,lv)*aux_rd(2,lv,mv,ikm)
          !
          sdt%fcamp(ikm,lv,mv) = (0.5_rk/electron_mass) * gc
        end do mval_loop
      end do lval_loop
    end do coulomb_kvec_magnitude
    !$omp end parallel do
    !
    call TimerStop('t-SURF: Direct: Coulomb')
  end subroutine sts_atend_direct_coulomb
  !
  subroutine sts_atend_direct_volkov(sdt,aux_rd)
    type(sts_data), intent(inout) :: sdt             ! t-SURF data structure.
    complex(rk), intent(in)       :: aux_rd(2,0:sd_lmax,sd_mmin:sd_mmax,sts_kgrid_count)
                                                     ! Solutions of auxiliary linear systems and derivatives
                                                     ! at the matching point, for each L, M, and K magnitude
                                                     ! Index 1: 1 = function, 2 = derivative
                                                     ! Index 2: 0:sd_lmax
                                                     ! Index 3: sd_mmin:sd_mmax
                                                     ! Index 4: 1 .. sts_kgrid_count
    !
    integer(ik)              :: alloc
    integer(ik)              :: ikd, ikm
    integer(ik)              :: lv, mv, mmin, mmax
    real(rk)                 :: gfacr, gfacd
    complex(rk)              :: gfac, gc
    complex(rk)              :: milp        ! -(-I)**(L)
    complex(rk), allocatable :: ylm (:,:)   ! Spherical harmonics; First index: M; Second index: L
    !
    call TimerStart('t-SURF: Direct: Volkov')
    !
    !$omp parallel default(none) &
    !$omp& private(ylm,alloc,ikd,lv,milp,mmin,mmax,ikm,gfacr,gfacd,gc,mv,gfac) &
    !$omp& shared(sts_dgrid_count,sd_lmax,sd_mmin,sd_mmax,sts_dtab,sts_kgrid_count) &
    !$omp& shared(sts_rmatch,sts_ktab,sts_btab,sdt,aux_rd,nts)
    allocate (ylm(sd_mmin:sd_mmax,0:sd_lmax),stat=alloc)
    if (alloc/=0) then
      stop 'spherical_tsurf%sts_atend_direct_volkov - error allocating per-thread t-SURF buffers'
    end if
    !
    !  We evaluate spherical harmonics too many times in the parallel loop below. Strictly speaking,
    !  it is enough to call MathAllYLM2 once per direction. However, MathAllYLM2 is quite efficient,
    !  so that we gain more by having more parallelism than we lose by repeating evaluation of the
    !  harmonics.
    !
    !$omp do collapse(2)
    kvec_directions: do ikd=1,sts_dgrid_count
      lval_loop: do lv=0,sd_lmax
        if (nts%dist(lv)>0) cycle lval_loop ! Hook for distributed-memory parallelization
        mmin = max(sd_mmin,-lv)
        mmax = min(sd_mmax, lv)
        if (mmin>mmax) cycle lval_loop
        !
        milp = -((0._rk,-1._rk)**(lv))
        !
        !  Evaluate spherical harmonics
        !
        call MathAllYLM2(lv,mmin,mmax,sts_dtab(:,ikd),ylm(mmin:mmax,0:lv),phase='Arfken')
        !
        !  Go over all possible values of K magnitude
        !
        kvec_magnitude: do ikm=1,sts_kgrid_count
          !
          !  Collect all terms which do not depend of the M quantum number or the state identity
          !  
          !   gfacr is the prefactor for the wavefunction at the boundary;
          !   gfacd is the prefactor for the derivative
          !
          gfacr = (sts_rmatch/electron_mass) * sts_ktab(ikm) * (lv*sts_btab(lv+1,ikm) - (lv+1)*sts_btab(lv-1,ikm))/(4*lv+2)
          gfacd = (sts_rmatch/electron_mass) * sts_btab(lv,ikm) * 0.5_rk
          !
          !  We finally have all the ingredients for calculating the amplitude. Hurray!
          !
          gc    = 0
          mval_loop: do mv=mmin,mmax
            gfac = sqrt(2._rk/pi) * sdt%pref(ikm,ikd) * milp * ylm(mv,lv) * (gfacr*aux_rd(1,lv,mv,ikm) + gfacd*aux_rd(2,lv,mv,ikm))
            gc   = gc + gfac
          end do mval_loop
          !$omp atomic
          sdt%amplitude(ikm,ikd) = sdt%amplitude(ikm,ikd) + gc
        end do kvec_magnitude
      end do lval_loop
    end do kvec_directions
    !$omp end do
    deallocate (ylm)
    !$omp end parallel
    !
    call TimerStop('t-SURF: Direct: Volkov')
  end subroutine sts_atend_direct_volkov
  !
  !  Direct implementation of the at-end term follows Eq. 17 of the paper, avoiding
  !  evaluation of the spectral form of the Hamiltonian. To achieve this, we must
  !  solve an auxiliary linear system for each final wavevector magnitude and L,M 
  !  channel. Note that this part does not require the left wavefunction!
  !
  subroutine sts_atend_direct(sdt,wfn_r)
    type(sts_data), intent(inout) :: sdt    ! t-SURF data structure. 
    type(sd_wfn), intent(in)      :: wfn_r  ! Right wavefunction
    !
    integer(ik)              :: alloc
    complex(rk), allocatable :: aux_rd(:,:,:,:) ! Solutions of auxiliary linear systems and derivatives
                                                ! at the matching point, for each L, M, and K magnitude
                                                ! Index 1: 1 = function, 2 = derivative
                                                ! Index 2: 0:sd_lmax
                                                ! Index 3: sd_mmin:sd_mmax
                                                ! Index 4: 1 .. sts_kgrid_count
    !
    if (.not.sts_active_direct) return
    !
    if (sd_nspin/=1) stop 'spherical_tsurf%sts_atend_direct - sd_nspin must be 1'
    !
    call TimerStart('t-SURF: At end: direct')
    !
    !  Precompute data we'll need frequently later on
    !
    call sts_atend_prepare_volkov(sdt)
    !
    call sts_atend_prepare_coulomb(sdt)
    !
    allocate (aux_rd(2,0:sd_lmax,sd_mmin:sd_mmax,sts_kgrid_count),stat=alloc)
    if (alloc/=0) then
      stop 'spherical_tsurf%sts_atend_direct - allocation failed'
    end if
    !
    call sts_atend_direct_auxiliary(sdt,wfn_r,aux_rd)
    !
    if (sts_coulomb_atend) then
      call sts_atend_direct_coulomb(sdt,aux_rd)
    end if
    !
    if (sts_volkov_atend) then
      call sts_atend_direct_volkov(sdt,aux_rd)
    end if
    !
    deallocate (aux_rd)
    call TimerStop('t-SURF: At end: direct')
  end subroutine sts_atend_direct
  !
  !  WARNING: sts_atend_terms() will be called from inside an OpenMP
  !  WARNING: parallel region. It therefore MUST protect all updates to sdt%amplitude
  !  WARNING: with a lock! (It is not necessary to lock fcgr() - is is accessed by L)
  !
  !  It may still be beneficial to run it in parallel: the number of calling
  !  threads nay be limited by the amount of memory available.
  !
  subroutine sts_atend_terms(sdt,lval,en,evec,amp)
    type(sts_data), intent(inout)  :: sdt                              ! t-SURF data structure. LOCK THE UPDATES!
    integer(ik), intent(in)        :: lval                             ! Angular momentum L
    complex(rk), intent(in)        :: en  (sd_nradial)                 ! Eigenvalues of the field-free states
    complex(rk), intent(in)        :: evec(sd_nradial,sd_nradial)      ! Right eigenvectors of the field-free states
    complex(rk), intent(in)        :: amp (sd_nradial,sd_mmin:sd_mmax) ! Amplitudes of the field-free states in the 
                                                                       ! TDSE solutions at the time of the analysis
    !
    complex(rk) :: evrd(2,sd_nradial)   ! Stationary solutions and their gradients at the matching point
                                        ! First index: 1=function; 2=gradient; second index = solution ordinal
    !
    if (.not.sts_active_atend) return
    !
    if (sd_nspin/=1) stop 'spherical_tsurf%sts_atend_terms - sd_nspin must be 1'
    !
    call TimerStart('t-SURF: At end')
    !
    !  Step numbering below matches sts_timestep; however, we end up skipping a lot
    !  of steps because of: a) absence of the field; and b) pure L, M values.
    !
    !  Step 1: Extract wavefunctions and radial derivatives at the matching point
    !          for all stationary solutions.
    !
    call fill_stationary_solutions_at_matching_point(sdt%wfn_scale,lval,evec,evrd)
    !
    if (sts_volkov_atend) then
      call volkov_atend_terms(sdt,lval,en,amp,evrd)
    end if
    !
    !  Calculate amplitudes for the Coulomb scattering waves; this should cost almost nothing 
    !  compared to the Volkov solutions above.
    !
    if (sts_coulomb_atend) then
      call coulomb_atend_terms(sdt,lval,en,amp,evrd)
    end if
    !
    call TimerStop('t-SURF: At end')
  end subroutine sts_atend_terms
  !
  subroutine fill_stationary_solutions_at_matching_point(scale,lval,evec,evrd)
    real(rk), intent(in)     :: scale                       ! Scaling factor for converting right wavefunction into the real-space wavefunction
    integer(ik), intent(in)  :: lval
    complex(rk), intent(in)  :: evec(sd_nradial,sd_nradial) ! Stationary solutions. (Right eigenvectors)
    complex(rk), intent(out) :: evrd(2,sd_nradial)          ! Solution and gradient at the matching point
    !
    complex(rk)              :: tmp (sd_nradial,sts_atend_block)
    complex(rk)              :: grad(sd_nradial,sts_atend_block)  
    complex(rk)              :: scr (sd_nradial,sts_atend_block*m3d_sc_size) 
    integer(ik)              :: iv, ie, ns
    !
    call TimerStart('t-SURF: At end: Matching point')
    if (sts_ipt_match<=1 .or. sts_ipt_match>=sd_nradial) then
      stop 'spherical_tsurf%fill_stationary_solutions_at_matching_point - bad sts_ipt_match'
    end if
    if (sts_atend_block<=0) stop 'spherical_tsurf%fill_stationary_solutions_at_matching_point - bad sts_atend_block'
    !
    !$omp parallel do default(none) &
    !$omp& shared(sd_nradial,sts_atend_block,evec,evrd,sts_ipt_match,lval) &
    !$omp& shared(sd_d1n_l0,sd_m1n_l0,sd_m1nf_l0,sd_m1np_l0) &
    !$omp& shared(sd_d1n_lx,sd_m1n_lx,sd_m1nf_lx,sd_m1np_lx) &
    !$omp& shared(scale) &
    !$omp& private(iv,ie,ns,tmp,grad,scr)
    solution_blocks: do iv=1,sd_nradial,sts_atend_block
      ie = min(sd_nradial,iv+sts_atend_block-1)
      ns = ie - iv + 1
      !
      if (lval==0) then
        call m3d_multiply(sd_d1n_l0,evec(:,iv:ie),tmp(:,1:ns))
        call m3d_solve(sd_m1n_l0,sd_m1nf_l0,sd_m1np_l0,tmp(:,1:ns),grad(:,1:ns),scr(:,1:m3d_sc_size*ns))
      else
        call m3d_multiply(sd_d1n_lx,evec(:,iv:ie),tmp(:,1:ns))
        call m3d_solve(sd_m1n_lx,sd_m1nf_lx,sd_m1np_lx,tmp(:,1:ns),grad(:,1:ns),scr(:,1:m3d_sc_size*ns))
      end if
      evrd(1,iv:ie) = scale * evec(sts_ipt_match,iv:ie)
      evrd(2,iv:ie) = scale * grad(sts_ipt_match,1:ns)
    end do solution_blocks
    !$omp end parallel do
    !
    if (sts_verbose>=2) then
      write (out,"(/t5,'Stationary solution and gradient at the matching point'/)")
      write (out,"((1x,a3,1x,a4,2(2x,a24,1x,a24)))") &
             ' L ', ' I ', ' Re[Psi] ', ' Im[Psi] ', ' Re[d Psi/d R] ', ' Im[d Psi/d R] ', &
             '---', '---', '---------', '---------', '---------------', '---------------'
      print_state: do iv=1,sd_nradial
        write (out,"(1x,i3,1x,i4,2(2x,g24.16,1x,g24.16))") lval, iv, evrd(:,iv)
      end do print_state
      write (out,"()")
    end if
    call TimerStop('t-SURF: At end: Matching point')
  end subroutine fill_stationary_solutions_at_matching_point
  !
  subroutine volkov_atend_terms(sdt,lval,en,amp,evrd)
    type(sts_data), intent(inout)  :: sdt                              ! t-SURF data structure. LOCK THE UPDATES!
    integer(ik), intent(in)        :: lval                             ! Angular momentum L
    complex(rk), intent(in)        :: en  (sd_nradial)                 ! Eigenvalues of the field-free states
    complex(rk), intent(in)        :: amp (sd_nradial,sd_mmin:sd_mmax) ! Amplitudes of the field-free states in the 
                                                                       ! TDSE solutions at the time of the analysis
    complex(rk), intent(in)        :: evrd(2,sd_nradial)               ! Stationary solutions and their gradients at the matching point
                                                                       ! First index: 1=function; 2=gradient; second index = solution ordinal
    !
    integer(ik)              :: alloc
    integer(ik)              :: ikd, ikm
    integer(ik)              :: mmin, mmax
    integer(ik)              :: is
    real(rk)                 :: gfacr, gfacd
    real(rk)                 :: ekin
    complex(rk)              :: gfac, gc
    complex(rk)              :: ilp1        ! (-I)**(L+1)
    complex(rk), allocatable :: ylm (:,:)   ! Spherical harmonics; First index: M; Second index: L
    complex(rk), allocatable :: accu(:,:)   ! Amplitude accumulator, to avoid locking (sdt) more than once
                                            ! First index: k magnitude; Second index: k direction (same as sdt%amplitude)
    !
    mmin = max(sd_mmin,-lval)
    mmax = min(sd_mmax, lval)
    if (mmin>mmax) return ! Nothing to do
    !
    call TimerStart('t-SURF: At end: Volkov')
    !
    allocate (accu(sts_kgrid_count,sts_dgrid_count),stat=alloc)
    if (alloc/=0) then
      stop 'spherical_tsurf%volkov_atend_terms - error allocating global t-SURF buffers'
    end if
    !
    !  Step 4 through 6 are repeated for all k vector directions
    !
  ! ilp1 = (0._rk,-1._rk)**(lval+1)
    ilp1 = mipow(modulo(lval+1,4))
    accu = 0
    !$omp parallel default(none) &
    !$omp& shared(mmin,mmax,lval,sts_dgrid_count,sts_dtab,sts_kgrid_count) &
    !$omp& shared(sts_rmatch,sts_ktab,sts_btab,sd_nradial,sdt,evrd) &
    !$omp& shared(ilp1,accu,amp,en) &
    !$omp& private(ylm,alloc,ikd,ikm,gfacr,gfacd,ekin,is,gfac,gc)
    allocate (ylm(mmin:mmax,0:lval),stat=alloc)
    if (alloc/=0) then
      stop 'spherical_tsurf%volkov_atend_terms - error allocating per-thread t-SURF buffers'
    end if
    !$omp do
    kvec_directions: do ikd=1,sts_dgrid_count
      !
      !  Step 4: evaluate spherical harmonics
      !
      call MathAllYLM2(lval,mmin,mmax,sts_dtab(:,ikd),ylm,phase='Arfken')
      !
      !  Step 6: go over all possible values of K magnitude
      !
      kvec_magnitude: do ikm=1,sts_kgrid_count
        !
        !  Collect all terms which do not depend of the M quantum number or the state identity
        !  
        !   gfacr is the prefactor for the wavefunction at the boundary;
        !   gfacd is the prefactor for the derivative
        !
        gfacr = (sts_rmatch/electron_mass) * sts_ktab(ikm) * (lval*sts_btab(lval+1,ikm) - (lval+1)*sts_btab(lval-1,ikm))/(4*lval+2)
        gfacd = (sts_rmatch/electron_mass) * sts_btab(lval,ikm) * 0.5_rk
        ekin  = sts_ktab(ikm)**2/(2._rk*electron_mass)
        !
        !  We finally have all the ingredients for calculating the amplitude. Hurray!
        !
        state_loop: do is=1,sd_nradial
          !
          !  Assemble all factors independent of the angular momentum projection M
          !
          gfac = (0._rk,1._rk) * sqrt(2._rk/pi) * sdt%pref(ikm,ikd) * ilp1 * (gfacr*evrd(1,is) + gfacd*evrd(2,is))
          !
        ! mloop: do mv=mmin,mmax
        !   gc = gfac * ylm(mv,lval) * amp(is,mv)
        !   accu(ikm,ikd) = accu(ikm,ikd) + gc / ( ekin - en(is) )
        ! end do mloop
          gc = gfac * sum(ylm(mmin:mmax,lval) * amp(is,mmin:mmax))
          accu(ikm,ikd) = accu(ikm,ikd) + gc / ( ekin - en(is) )
        end do state_loop
      end do kvec_magnitude
    end do kvec_directions
    !$omp end do
    deallocate (ylm)
    !$omp end parallel
    !
    !$omp critical
    sdt%amplitude = sdt%amplitude + accu
    !$omp end critical
    deallocate (accu)
    call TimerStop('t-SURF: At end: Volkov')
  end subroutine volkov_atend_terms
  !
  subroutine coulomb_atend_terms(sdt,lval,en,amp,evrd)
    type(sts_data), intent(inout)  :: sdt                              ! t-SURF data structure. LOCK THE UPDATES!
    integer(ik), intent(in)        :: lval                             ! Angular momentum L
    complex(rk), intent(in)        :: en  (sd_nradial)                 ! Eigenvalues of the field-free states
    complex(rk), intent(in)        :: amp (sd_nradial,sd_mmin:sd_mmax) ! Amplitudes of the field-free states in the 
                                                                       ! TDSE solutions at the time of the analysis
    complex(rk), intent(in)        :: evrd(2,sd_nradial)               ! Stationary solutions and their gradients at the matching point
                                                                       ! First index: 1=function; 2=gradient; second index = solution ordinal
    !
    integer(ik)              :: ikm
    integer(ik)              :: mmin, mmax
    integer(ik)              :: mv, is
    real(rk)                 :: ekin
    complex(rk)              :: gfac, gc
    !
    mmin = max(sd_mmin,-lval)
    mmax = min(sd_mmax, lval)
    if (mmin>mmax) return ! Nothing to do
    !
    call TimerStart('t-SURF: At end: Coulomb')
    !
    !  Calculate amplitudes for the Coulomb scattering waves; this should cost almost nothing 
    !  compared to the Volkov solutions above.
    !
    ! sdt%fcamp(:,lval,:) = 0
    !$omp parallel do default(none) &
    !$omp& shared(sts_kgrid_count,sdt,sts_ktab,evrd,amp,mmin,mmax) &
    !$omp& shared(en,sd_nradial,lval) &
    !$omp& private(ikm,is,ekin,gfac,gc)
    coulomb_kvec_magnitude: do ikm=1,sts_kgrid_count
      sdt%fcamp(ikm,lval,:) = 0
      ekin  = sts_ktab(ikm)**2/(2._rk*electron_mass)
      coulomb_state_loop: do is=1,sd_nradial
        !
        !  gfac = (I hbar)/(2 Me) * [ (d chi*/d r) psi - chi* (d psi/d r) ]
        !
        !  where chi is the Coulomb scattering solution (real, so conjugation is immaterial)
        !  and psi is one of the stationary solutions.
        !
        gfac = sdt%fcgr(2,ikm,lval)*evrd(1,is) - sdt%fcgr(1,ikm,lval)*evrd(2,is)
        gfac = ((0._rk,0.5_rk)/electron_mass) * gfac
        gfac = ((0._rk,1.0_rk)/(ekin-en(is))) * gfac
        !
        coulomb_mloop: do mv=mmin,mmax
          gc = gfac * amp(is,mv)
          sdt%fcamp(ikm,lval,mv) = sdt%fcamp(ikm,lval,mv) + gc
        end do coulomb_mloop
      end do coulomb_state_loop
    end do coulomb_kvec_magnitude
    !$omp end parallel do
    call TimerStop('t-SURF: At end: Coulomb')
  end subroutine coulomb_atend_terms
  !
  subroutine sts_report(sdt,mode)
    type(sts_data), intent(inout) :: sdt
    character(len=*), intent(in)  :: mode ! Either 'In TDSE' or 'At end'
    !
    !
    call TimerStart('t-SURF: Report')
    !
    call nt_merge_all(sdt)  ! Reporting phase runs on the master node only
    if (nts%this_node/=1) then
      call TimerStop('t-SURF: Report')
      return
    end if
    !
    select case (mode)
      case default
        stop 'spherical_tsurf%sts_report - bad mode parameter'
      case ('In TDSE')
        if (sts_volkov .and. .not. sts_volkov_atend) then
          write (out,"(/t5,'Projection on Volkov states is calculated from boundary flux up to now (t-SURF)')")
          write (out,"( t5,'Photoelectrons still within the simulation volume will not be included in the spectrum'/)")
          call table_report(sdt%volkov_table,sdt%amplitude,sdt%vphase)
          call opendx_report(sdt%volkov_opendx,sdt%amplitude)
        end if
      case ('At end')
        if (sts_volkov_atend) then
          if (sts_volkov) then
            write (out,"(/t5,'Projection on Volkov states is calculated from boundary flux until infinite time'/)")
          else
            write (out,"(/t5,'Projection on Volkov states is calculated from boundary flux from now until infinite time')")
            write (out,"( t5,'Photoelectrons reaching the boundary up to now will not be included in the spectrum'/)")
          end if
          call table_report(sdt%volkov_table,sdt%amplitude,sdt%vphase)
          call opendx_report(sdt%volkov_opendx,sdt%amplitude)
        end if
        if (sts_coulomb_atend) then
          write (out,"(/t5,'Projection on Coulomb states is calculated from boundary flux from now until infinite time')")
          write (out,"( t5,'Photoelectrons reaching the boundary up to now will not be included in the spectrum'/)")
          call phase_report(sdt%coulomb_waves,sdt%fcamp)
          call convert_coulomb_phase_to_spectrum(sdt)
          call table_report(sdt%coulomb_table,sdt%coulamp)
          call opendx_report(sdt%coulomb_opendx,sdt%coulamp)
        end if
    end select
    !
    call TimerStop('t-SURF: Report')
  end subroutine sts_report
  !
  subroutine convert_coulomb_phase_to_spectrum(sdt)
    type(sts_data), intent(inout) :: sdt
    !
    integer(ik)              :: ikm, ikd, lv, mv, mmin, mmax
    integer(ik)              :: alloc
    real(rk)                 :: cp, kv
    complex(rk)              :: accu
    complex(rk), allocatable :: ylm(:,:)      ! Spherical harmonics for each direction
    complex(rk), allocatable :: cphase(:,:)   ! Phase prefactor, containing:
                                              ! sqrt(8/pi) (1/k) (-I)**L Exp(-I cphase_L)
    !
    call TimerStart('t-SURF: Coulomb spectrum')
    allocate (cphase(0:sd_lmax,sts_kgrid_count),stat=alloc)
    if (alloc/=0) then
      stop 'spherical_tsurf%convert_coulomb_phase_to_spectrum - allocation failed (1)'
    end if
    !
    !$omp parallel do default(none) &
    !$omp& shared(sts_kgrid_count,sd_lmax,sts_ktab,sts_asymptotic_q,cphase) &
    !$omp& private(ikm,lv,kv,cp)
    fill_coulomb_phases_k: do ikm=1,sts_kgrid_count
      fill_coulomb_phases_l: do lv=0,sd_lmax
        !
        !  Arg(Gamma(L+1-I Z /k)). We need to be consistent with the invocation of coulombF here
        !
        kv = sts_ktab(ikm)
        !
        !  A bit of a problem here: Our derivation indicates sts_asymptotic_q*electron_charge/kv
        !  for the imaginary part of the gamma; however, the results are internally consistent
        !  with the Volkov function projection when it is -sts_asymptotic_q*electron_charge/kv.
        !
        !  Oops.
        !
        ! cp = aimag(MathLogGamma(cmplx(lv+1,sts_asymptotic_q*electron_charge/kv,kind=rk)))
        cp = aimag(MathLogGamma(cmplx(lv+1,-sts_asymptotic_q*electron_charge/kv,kind=rk)))
        !
        cphase(lv,ikm) = (1._rk/sqrt(2._rk*pi)) * (1._rk/kv) * (0._rk,-1._rk)**lv * exp((0._rk,-1._rk)*cp)
      end do fill_coulomb_phases_l
    end do fill_coulomb_phases_k
    !$omp end parallel do
    !
    !$omp parallel default(none) &
    !$omp& shared(sts_dgrid_count,sd_lmax,sts_dtab,sts_kgrid_count,sd_mmin,sd_mmax,cphase,sdt) &
    !$omp& private(alloc,ikd,ylm,ikm,accu,lv,mmin,mmax,mv)
    allocate (ylm(sd_mmin:sd_mmax,0:sd_lmax),stat=alloc)
    if (alloc/=0) then
      stop 'spherical_tsurf%convert_coulomb_phase_to_spectrum - allocation failed (2)'
    end if
    !$omp do schedule(guided,1)
    process_directions: do ikd=1,sts_dgrid_count
      !
      !  Spherical harmonics
      !
      call MathAllYLM2(sd_lmax,sd_mmin,sd_mmax,sts_dtab(:,ikd),ylm,phase='Arfken')
      !
      process_magnitudes: do ikm=1,sts_kgrid_count
        accu = 0
        sum_lv: do lv=0,sd_lmax
          mmin = max(sd_mmin,-lv)
          mmax = min(sd_mmax, lv)
          sum_mv: do mv=mmin,mmax
            accu = accu + cphase(lv,ikm) * ylm(mv,lv) * sdt%fcamp(ikm,lv,mv)
          end do sum_mv
        end do sum_lv
        sdt%coulamp(ikm,ikd) = accu
      end do process_magnitudes
    end do process_directions
    !$omp end do
    deallocate (ylm)
    !$omp end parallel
    !
    deallocate (cphase)
    !
    call TimerStop('t-SURF: Coulomb spectrum')
  end subroutine convert_coulomb_phase_to_spectrum
  !
  subroutine table_report(file,amp,phase)
    character(len=*), intent(in)   :: file        ! Destigation; use ' ' for standard output
    complex(rk), intent(in)        :: amp   (:,:) ! Amplitudes
    real(rk), intent(in), optional :: phase (:,:) ! Volkov phases at the time of the report
                                                  ! Default is zero
    !
    integer     :: iu_out
    integer(ik) :: ikm, ikd
    real(rk)    :: kvec(3), p_phase, th, ph
    !
    if (file==' ') then
      iu_out = out
    else
      open (newunit=iu_out,form='formatted',status='replace',file=trim(file))
    end if
    !
    write (iu_out,"('#')")
    write (iu_out,"('# The velocity-gauge Volkov states are taken as:')")
    write (iu_out,"('#')")
    write (iu_out,"('#  chi(k) = (2*Pi)**(-1.5) * Exp(I*(K.R) + Phi0)')")
    write (iu_out,"('#')")
    write (iu_out,"(('#',6(1x,a15),1x,3(1x,a25)))") &
      ' Kx ', ' Ky ', ' Kz ', ' |K| ', ' Theta ', ' Phi ', '  Re[A(K)]  ', '  Im[A(K)]  ', ' Phi0 ', &
      '----', '----', '----', '-----', '-------', '-----', '------------', '------------', '------'
    report_directions: do ikd=1,sts_dgrid_count
      th = acos(sts_dtab(3,ikd))*180._rk/pi
      ph = atan2(sts_dtab(2,ikd),sts_dtab(1,ikd))*180._rk/pi
      report_magnitudes: do ikm=1,sts_kgrid_count
        kvec(:) = sts_ktab(ikm) * sts_dtab(:,ikd)
        p_phase = 0
        if (present(phase)) p_phase = phase(ikm,ikd)
        write (iu_out,"(6(1x,g15.8),1x,3(1x,g25.14e3))") kvec, sts_ktab(ikm), th, ph, amp(ikm,ikd), p_phase
      end do report_magnitudes
      write (iu_out,"()")
    end do report_directions
    write (iu_out,"()")
    !
    if (file/=' ') then
      close (iu_out)
      write (out,"(/t5,'Photoelectron spectrum written to ',a)") trim(file)
    end if
  end subroutine table_report
  !
  !  Generate an OpenDX file containing the spectrum.
  !
  subroutine opendx_report(file,amp)
    character(len=*), intent(in) :: file     ! Output file to use; ' ' to supress output
    complex(rk), intent(in)      :: amp(:,:) ! Amplitudes
    !
    integer(ik) :: ikm, ikd, ith, iph, ith_eff, iph_eff
    integer(ik) :: icol, iu_odx
    !
    if (file==' ') return
    if (sts_dgrid/='product') then
      write (out,"(/'OpenDX output requires sts_dgrid=product'/)")
      return
    end if
    !
    call TimerStart('t-SURF: Report: OpenDX')
    !
    write (out,"(t5,'Generating OpenDX file for photoelectron spectrum: ',a)") trim(file)
    open (newunit=iu_odx,form='formatted',status='replace',file=trim(file))
    !
    !  Prepare OpenDX header describing our grid.
    !
    !  Object 1: regular product grid for the angular part: (theta,phi)
    !            Add two extra points for both theta and phi, one on each end (wrap-around)
    !            Adding extra theta points yields artifacts in OpenDX; don't do it.
    !
    write (iu_odx,"('object 1 class gridpositions counts ',i0,1x,i0,1x,i0)") sts_dgrid_ntheta, sts_dgrid_nphi + 2, 1
    write (iu_odx,"(' origin   ',3(1x,f14.8))") 0.5_rk*pi/sts_dgrid_ntheta, -0.5_rk*2._rk*pi/sts_dgrid_nphi, 0.0_rk
    write (iu_odx,"(' delta    ',3(1x,f14.8))") pi/sts_dgrid_ntheta, 0._rk, 0._rk
    write (iu_odx,"(' delta    ',3(1x,f14.8))") 0._rk, 2._rk*pi/sts_dgrid_nphi, 0._rk
    !
    !  Object 2: an irregular array for the radial part: (k)
    !
    write (iu_odx,"('object 2 class array type float rank 1 shape 3 items ',i0,' data follows')") sts_kgrid_count
    write (iu_odx,"((' 0.0 0.0 ',g18.9e3))") sts_ktab(:)
    !
    !  Object 3: product grid: regular theta, phi; irregular k
    !
    write (iu_odx,"('object 3 class productarray')")
    write (iu_odx,"(' term 1')")
    write (iu_odx,"(' term 2')")
    !
    !  Object 4: grid connections
    !
    write (iu_odx,"('object 4 class gridconnections counts ',3(1x,i0))") sts_dgrid_ntheta, sts_dgrid_nphi+2, sts_kgrid_count
    !
    !  Object 5: Photoelectron spectrum
    !
    write (iu_odx,"('object 5 class array type float category complex rank 0 items ',i0,' data 0')") &
           (sts_dgrid_ntheta)*(sts_dgrid_nphi+2)*sts_kgrid_count
    write (iu_odx,"('attribute ""dep"" string ""positions""')")
    !
    !  Overall data field
    !
    write (iu_odx,"('object ""field0"" class field')")
    write (iu_odx,"('component ""data"" value 5')")
    write (iu_odx,"('component ""positions"" value 3')")
    write (iu_odx,"('component ""connections"" value 4')")
    write (iu_odx,"('attribute ""name"" string ""field0""')")
    write (iu_odx,"('end')")
    !
    !  Done with the header; write the data
    !
    icol = 0
    directions_th: do ith=1,sts_dgrid_ntheta
      ! Extra code for ith=0 and its=(sts_dgrid_ntheta+1) is to let positions to wrap-around if we want to do it
      if (ith==0) then
        ith_eff = sts_dgrid_ntheta
      else if (ith==sts_dgrid_ntheta+1) then
        ith_eff = 1
      else
        ith_eff = ith
      end if
      directions_ph: do iph=0,sts_dgrid_nphi+1
        if (iph==0) then
          iph_eff = sts_dgrid_nphi
        else if (iph==sts_dgrid_nphi+1) then
          iph_eff = 1
        else
          iph_eff = iph
        end if
        ikd = iph_eff+(ith_eff-1)*sts_dgrid_nphi
        magnitudes: do ikm=1,sts_kgrid_count
          call write_output(real(amp(ikm,ikd),kind=rk))
          call write_output(aimag(amp(ikm,ikd)))
        end do magnitudes
      end do directions_ph
    end do directions_th
    if (icol/=0) write (iu_odx,"()")
    close (iu_odx)
    !
    call TimerStop('t-SURF: Report: OpenDX')
    !
    contains
      subroutine write_output(val)
        real(rk), intent(in) :: val
        !
        write (iu_odx,"(g18.9e3,1x)",advance='no') val
        icol = icol + 1
        if (icol==6) then
          icol = 0
          write (iu_odx,"()")
        end if
     end subroutine write_output
  end subroutine opendx_report
  !
  subroutine phase_report(file,fcamp)
    character (len=*), intent(in) :: file         ! File to use for the report; ' ' for standart output
    complex(rk), intent(in)       :: fcamp(sts_kgrid_count,0:sd_lmax,sd_mmin:sd_mmax) 
                                                  ! Amplitides of the Coulomb functions
    !
    integer     :: iu_out
    integer(ik) :: ikm, lv, mv, mmin, mmax
    !
    if (file==' ') then
      iu_out = out
    else
      open (newunit=iu_out,form='formatted',status='replace',file=trim(file))
    end if
    !
    write (iu_out,"('#')")
    write (iu_out,"('# Amplitudes of the spherical Coulomb continuum states (Z=',f12.6,')')") sts_asymptotic_q
    write (iu_out,"('# Asymptotic form:')")
    write (iu_out,"('# ')")
    write (iu_out,"('#    (2/r) Sin( k r + (Z/k)*log(2 k r) - l pi/2 + phi)')")
    write (iu_out,"('# ')")
    write (iu_out,"('#    phi = Arg( Gamma(l + 1 - I Z/k) )')")
    write (iu_out,"('#')")
    write (iu_out,"(('#',1x,a5,1x,a6,3(1x,a25)))") &
           ' L ', ' M ', ' |K| ', ' Re[amp] ', ' Im[amp] ', &
           '---', '---', '-----', '---------', '---------'
    report_l: do lv=0,sd_lmax
      mmin = max(sd_mmin,-lv)
      mmax = min(sd_mmax, lv)
      report_m: do mv=mmin,mmax
        report_magnitudes: do ikm=1,sts_kgrid_count
          write (iu_out,"(1x,i5,1x,i6,3(1x,g25.14e3))") lv, mv, sts_ktab(ikm), fcamp(ikm,lv,mv)
        end do report_magnitudes
        write (iu_out,"()")
      end do report_m
      write (iu_out,"()")
    end do report_l
    write (iu_out,"()")
    !
    if (file/=' ') then
      close (iu_out)
      write (out,"(/t5,'Coulomb-wave phases are written to ',a)") trim(file)
    end if
  end subroutine phase_report
  !
end module spherical_tsurf
