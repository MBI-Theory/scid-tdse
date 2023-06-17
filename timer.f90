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
module timer
!
!  Timing routines. Because of the possibility for parallel
!  execution, potentially on two scales at the same time 
!  (OpenMP and MPI), all timings are done using real time.
!
!  For OpenMP programs, all calls to timer routines from 
!  secondary threads are silently ignored.
!
!  The externally visible routines are:
!
!    TimerStart  - Begins a timing region
!    TimerStop   - Ends a timing region
!    TimerReport - Produces timing report
!
!  Internally, timers are implemented using a hash table and
!  should have constant cost regardless of the number of
!  timers. Hash implementation follows Algorthim L (linear 
!  insertion with linear search) of chapter 6.4 of Knuth's 
!  vol 3. Hash function is taken from Aho, Sethi, and Ullman, 
!  pp. 434-438. Because timers are persistent, there is no
!  deletion.
!
  use accuracy
  implicit none
  private
  public TimerStart, TimerStop, TimerReport
  public timer_disable
  public rcsid_timer
!
  character(len=clen), save :: rcsid_timer = "$Id: timer.f90,v 1.6 2023/06/09 14:10:24 ps Exp $"
!
  integer, parameter     :: trk        = selected_real_kind(14)  ! Must be of default integer kind
  integer(ik), parameter :: table_size = 1000 ! Max number of entries to track
  integer(ik), parameter :: name_len   =   40 ! Max length of timer name
!
  type tim
    logical                 :: used       ! Slot used?
    logical                 :: active     ! Currently active?
    character(len=name_len) :: name       ! Timer name
    real(trk)               :: calls      ! Number of times the timer was invoked
                                          !
                                          ! All times below are in seconds
                                          !
    real(trk)               :: real_time  ! Total real time on this timer
    real(trk)               :: cpu_time   ! ditto for CPU time
    real(trk)               :: real_kids  ! Real time spent in nested timers
    real(trk)               :: cpu_kids   ! ditto for CPU time
    real(trk)               :: real_start ! For active timers, time of activation
    real(trk)               :: cpu_start  ! ditto for CPI time
    integer(ik)             :: stack_p    ! For active timers, position in the stack
  end type tim
!
!  Timers and sundry data
!
  type(tim), target :: t_table (table_size) ! All of our timers
  integer(ik)       :: t_nested(table_size) ! Stack of currently active timers
  integer(ik)       :: t_appear(table_size) ! Appearance order for the timers
  integer(ik)       :: t_count  = 0         ! Number of defined timers
  integer(ik)       :: t_active = 0         ! Number of currently active timers
  real(trk)         :: prog_start           ! Timebase for the real time
  real(trk)         :: count_backward = 0   ! Number of backward real-time steps
  real(trk)         :: count_overflow = 0   ! Number of real-time overflows
  real(trk)         :: count_bad_real = 0   ! Number of times time-ordering was violated for real time
  real(trk)         :: count_bad_cpu  = 0   ! Number of times time-ordering was violated for CPU time
  logical, save     :: timer_disable  = .false.
!
  contains
  !
  !  Public interfaces
  !
    !
    !  Start timing region. 
    !
    subroutine TimerStart(region)
      character(len=*), intent(in) :: region  ! Timer name
      !
      integer(ik)          :: pos  ! Timer position
      type(tim), pointer   :: t    ! Current timer (for convenience)
      !
      if (timer_disable) return
      if (omp_secondary()) return ! Timer routines are not thread-safe otherwise!
      !
      !  One-time initialization
      !
      if (t_count==0) then
        prog_start      = get_real_time()
        t_table(:)%used = .false.
      end if
      !
      pos =  insert_item(region)
      t   => t_table(pos)
      !
      if (t%active) then
        write (out,"('TimerStart: timer ',a,' is already active')") trim(region)
        stop 'TimerStart - nested timer'
      end if
      !
      !  Push the new timer to the timer stack
      !
      t_active = t_active + 1
      t_nested(t_active) = pos
      !
      t%active     = .true.
      t%stack_p    = t_active
      t%calls      = t%calls + 1
      t%real_start = get_real_time()
      t%cpu_start  = get_cpu_time ()
    end subroutine TimerStart
    !
    !  End timing region. 
    !
    subroutine TimerStop(region)
      character(len=*), intent(in) :: region  ! Timer name
      !
      integer(ik)        :: pos        ! Timer position
      type(tim), pointer :: t          ! Current timer (for convenience)
      type(tim), pointer :: pt         ! Parent timer
      real(trk)          :: real_time  ! Real time for this invocation
      real(trk)          :: cpu_time   ! ditto CPU
      !
      if (timer_disable) return
      if (omp_secondary()) return ! Timer routines are not thread-safe otherwise!
      !
      pos =  insert_item(region)
      t   => t_table(pos)
      !
      if (.not.t%active) then
        write (out,"('TimerStop: timer ',a,' is not running')") trim(region)
        stop 'TimerStop - inactive timer'
      end if
      !
      !  Get timings for this invocation
      !
      real_time = get_real_time() - t%real_start
      cpu_time  = get_cpu_time()  - t%cpu_start
      !
      if (real_time<0 .or. cpu_time<0) then
        write (out,"(' TIMER BUG DETECTED:')")
        write (out,"(' region     = ',a)") trim(region)
        write (out,"(' real_time  = ',g24.16)") real_time
        write (out,"(' cpu_time   = ',g24.16)") cpu_time
        write (out,"(' real_start = ',g24.16)") t%real_start
        write (out,"(' cpu_start  = ',g24.16)") t%cpu_start
      end if
      !
      !  Update time counts for the leaving timer
      !
      t%real_time = t%real_time + real_time
      t%cpu_time  = t%cpu_time + cpu_time
      !
      !  Mark as inactive
      !
      t%active    = .false.
      !
      !  Pop the timer from stack, and update counts for the parent
      !
      t_active = t_active - 1
      if (t_active>0) then
        pt => t_table(t_nested(t_active))
        pt%real_kids = pt%real_kids + real_time
        pt%cpu_kids  = pt%cpu_kids  + cpu_time
      end if
      !
    end subroutine TimerStop
    !
    !  Produce timing report
    !
    subroutine TimerReport
      real(trk)          :: real_now
      real(trk)          :: cpu_now 
      real(trk)          :: real_time, cpu_time 
      real(trk)          :: real_kids, cpu_kids
      real(trk)          :: real_threshold
      real(trk)          :: cpu_threshold
      character(len=19)  :: timestamp
      integer(ik)        :: ord
      integer(ik)        :: pos, kid_pos
      type(tim), pointer :: t, k
      character(len=1)   :: active
      integer(ik)        :: omitted
      !
      if (timer_disable) return
      if (omp_secondary()) return ! Timer routines are not thread-safe otherwise!
      !
      real_now  = get_real_time()
      cpu_now   = get_cpu_time()
      timestamp = get_timestamp()
      !
      real_threshold = 0.01_trk * (real_now - prog_start)
      cpu_threshold  = 0.01_trk * cpu_now
      !
      if (real_threshold<0 .or. cpu_threshold<0) then
        write (out,"('TIMER BUG DETECTED')")
        write (out,"(' real_threshold = ',g24.16)") real_threshold
        write (out,"(' cpu_threshold  = ',g24.16)") cpu_threshold
        write (out,"(' real_now       = ',g24.16)") real_now
        write (out,"(' prog_start     = ',g24.16)") prog_start
        write (out,"(' cpu_now        = ',g24.16)") cpu_now
      end if
      !
      write (out,"(/t15,'Timing data at ',a/)") timestamp
      write (out,"(t2,'     ',t42,'     ',t56,'Total time (seconds)',t87,'Self time (seconds)')")
      write (out,"(t2,'Timer',t42,'Calls',t56,'--------------------',t87,'-------------------')")
      write (out,"(t2,'-----',t42,'-----',t58,'Real',t73,'CPU',t88,'Real',t103,'CPU')")
      write (out,"()")
      !
      omitted = 0
      scan: do ord=1,t_count
        pos = t_appear(ord)
        t => t_table(pos)
        if (.not.t%used) then
          write (out,"('Timer ',i4,' in slot ',i5,' is defined but unused?!')") ord, pos
          stop 'TimerReport - logic error'
        end if
        !
        ! Calculate active-timer corrections
        ! 
        real_time = 0 ; real_kids = 0 ;
        cpu_time  = 0 ; cpu_kids  = 0 ;
        active     = ' '
        if (t%active) then
          real_time = real_now - t%real_start
          cpu_time  = cpu_now  - t%cpu_start
          if (t_active/=t%stack_p) then 
            !
            ! If we are not at the top of the stack, adjust
            ! cumulative children time.
            !
            kid_pos   = t_nested(t%stack_p+1)
            k         => t_table(kid_pos)
            real_kids = real_now - k%real_start
            cpu_kids  = cpu_now  - k%cpu_start
          end if
          active     = '*'
        end if
        !
        real_time = real_time + t%real_time
        cpu_time  = cpu_time  + t%cpu_time
        real_kids = real_kids + t%real_kids
        cpu_kids  = cpu_kids  + t%cpu_kids
        !
        !  Output needed?
        !
        if (real_time<real_threshold .and. cpu_time<cpu_threshold) then
          omitted = omitted + 1
          cycle scan
        end if
        !
        !  Output
        !
        write (out,"(t2,a30,t33,a1,t35,f12.0,t49,2(f13.1,1x,f13.1,3x))") &
               t%name, active, t%calls, real_time, cpu_time, &
               real_time - real_kids, cpu_time - cpu_kids
      end do scan
      if (omitted>0) then
        write (out,"(/' (',i3,' timers contributing less than 1% are not shown)')") &
               omitted
      else
      end if
      if (count_backward>0 .or. count_overflow>0 .or. count_bad_real>0 .or. count_bad_cpu>0) then
        write (out,"()")
        if (count_backward>0) write (out,"(t5,'Number of real-time decrements = ',f16.0)") count_backward
        if (count_overflow>0) write (out,"(t5,'including real-time roll-overs = ',f16.0)") count_overflow
        if (count_bad_real>0) write (out,"(t5,' Real-time ordering violations = ',f16.0)") count_bad_real
        if (count_bad_cpu >0) write (out,"(t5,'  CPU-time ordering violations = ',f16.0)") count_bad_cpu
      end if
      write (out,"()")
      call flush_wrapper(out)
    end subroutine TimerReport
  !
  !  Support routines
  !
    logical function omp_secondary()
    !$ use OMP_LIB
      integer :: this_thread  ! Must be of default integer kind
      !
      this_thread = 0
      !$ this_thread = omp_get_thread_num()
      omp_secondary = (this_thread /= 0)
    end function omp_secondary
    !
    !  Get real time in seconds. SYSTEM_CLOCK is allowed to roll over,
    !  so this code gets a little messy to handle the roll-over.
    !
    !  Furthermore, we have to make sure that get_real_time() is
    !  time-ordered; anything else will play havoc with timing 
    !  routines.
    !
    function get_real_time() result(t)
      real(trk) :: t
      !
      integer         :: count, count_rate, count_max ! Must be of default integer kind
      real(trk), save :: overflow   =  0
      integer, save   :: last_count = -1   ! Must be of default integer kind
      real(trk), save :: last_time  = -1   ! Initialize to an impossible small value
      !
      call system_clock(count,count_rate,count_max)
      !
      ! Try to detect a rollover. This is a little trickier than one
      ! might think: on systems with ntp-controlled time, there is a
      ! possibility for time resets. We do not want to see those as
      ! roll-overs!
      !
      if (count<last_count) then
        !? write (out,"(' # get_real_time: detected a backward step')")
        !? write (out,"(' # count = ',i20,' count_rate = ',i20,' count_max = ',i20)") &
        !?            count, count_rate, count_max
        !? write (out,"(' # oveflow = ',g20.12,' last_count = ',i20)") &
        !?            overflow, last_count
        count_backward = count_backward + 1
        if ( (last_count-count) > (count_max/2) ) then
          !? write (out,"(' # this is a roll-over, adjust counters')")
          count_overflow = count_overflow + 1
          overflow = overflow + real(count_max,kind=rk)
        end if
      end if
      last_count = count
      !
      ! Convert to seconds
      !
      t = (overflow+real(count,kind=rk))/real(count_rate,kind=rk)
      !
      if (t<last_time) then
        !
        !  Oops. The result violates time-ordering property; substitute the last time.
        !
        t = last_time
        count_bad_real = count_bad_real + 1
      else
        last_time = t
      end if
    end function get_real_time
    !
    !  Get CPU time, whatever this means
    !
    function get_cpu_time() result(t)
      real(trk)       :: t
      real(trk), save :: last_time = -1  ! An impossible small value
      !
      call cpu_time(t)
      !
      if (t<last_time) then
        !
        !  Oops. The result violates time-ordering property; substitute the last time.
        !
        t = last_time
        count_bad_cpu = count_bad_cpu + 1
      else
        last_time = t
      end if
    end function get_cpu_time
    !
    !  Return a (reasonably) pretty time stamp
    !
    function get_timestamp() result(s)
      character(len=19) :: s
      !
      character(len=8)  :: date
      character(len=10) :: time
      !
      call date_and_time(date=date,time=time)
      s = 'YYYY/MM/DD HH:MM:SS'
      s( 1: 4) = date(1:4)
      s( 6: 7) = date(5:6)
      s( 9:10) = date(7:8)
      s(12:13) = time(1:2)
      s(15:16) = time(3:4)
      s(18:19) = time(5:6)
    end function get_timestamp
    !
    !  Linear insertion with linear search (Algorithm L)
    !
    function insert_item(name) result(pos)
      character(len=*), intent(in) :: name
      integer(ik)                  :: pos
      !
      pos = string_hash(name)
      search: do
        if (.not.t_table(pos)%used) then
          !
          ! This is a new key, insert it
          !
          t_count = t_count + 1
          if (t_count>=table_size/5) then
            write (out,"('Too many timers. Increase table_size in "// &
                       "timer.f90 to at least ',i5)") t_count*5
            stop 'timer%insert_item'
          end if
          t_appear(t_count)      = pos
          t_table(pos)%used      = .true.
          t_table(pos)%active    = .false.
          t_table(pos)%name      = name
          t_table(pos)%calls     = 0
          t_table(pos)%real_time = 0
          t_table(pos)%cpu_time  = 0
          t_table(pos)%real_kids = 0
          t_table(pos)%cpu_kids  = 0
          exit search
        end if
        if (t_table(pos)%name==name) then
          !
          ! This is an existing key, simply return the location
          !
          exit search
        end if
        pos = 1 + modulo(pos-2,table_size)
      end do search
      !
    end function insert_item
    !
    !  Function below must use at least 29-bit integers. We'll stick
    !  to the default kind.
    !
    integer function string_hash(str) result(h)
      character(len=*), intent(in) :: str
!
      integer :: i, chr, g
      integer :: mask 
      data mask/Z"1FFFFFF"/
!
!    This hash assumes at least 29-bit integers. It is supposedly
!    documented in Aho, Sethi, and Ullman, pp. 434-438
!
      h = 0
      do i=1,len_trim(str)
        chr = ichar(str(i:i))
        h   = ishft(h,  4) + chr
        g   = ishft(h,-24)
        h   = iand(ieor(h,g),mask)
      end do
      h = 1 + modulo(h,int(table_size))
    end function string_hash
  !
end module timer
