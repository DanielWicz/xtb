! This file is part of xtb.
!
! Helper utilities to sample the current resident set size (RSS) from
! /proc/self/status. Designed for debugging memory growth during
! geometry optimizations; guarded by the XTB_MEMLOG environment
! variable to avoid overhead in normal production runs.
module xtb_mctc_meminfo
   use iso_fortran_env, only : int64
   implicit none
   private

   public :: rss_kb, memlog_enabled, log_memory_usage, log_memory_usage_delta

contains

   !> Returns .true. once XTB_MEMLOG is set in the environment.
   logical function memlog_enabled()
      logical, save :: initialized = .false.
      logical, save :: enabled = .false.
      integer :: len, stat

      if (.not.initialized) then
         call get_environment_variable('XTB_MEMLOG', length=len, status=stat)
         enabled = (stat == 0 .and. len > 0)
         initialized = .true.
      end if

      memlog_enabled = enabled
   end function memlog_enabled

   !> Read VmRSS from /proc/self/status. Returns -1 on failure.
   integer(int64) function rss_kb()
      integer :: unit, ios, pos
      character(len=256) :: line
      character(len=32) :: suffix

      rss_kb = -1_int64
      open(newunit=unit, file='/proc/self/status', status='old', action='read', &
         & iostat=ios)
      if (ios /= 0) return

      do
         read(unit,'(A)', iostat=ios) line
         if (ios /= 0) exit
         if (line(1:5) == 'VmRSS') then
            pos = index(line, ':')
            if (pos > 0) then
               read(line(pos+1:), *, iostat=ios) rss_kb, suffix
            end if
            exit
         end if
      end do

      close(unit)
      if (ios /= 0) rss_kb = -1_int64
   end function rss_kb

   !> Conditional RSS logger (enabled via XTB_MEMLOG).
   subroutine log_memory_usage(unit, label)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: label
      integer(int64) :: rss

      if (.not.memlog_enabled()) return
      rss = rss_kb()
      if (rss >= 0_int64) then
         write(unit,'(1x,"[mem]",1x,a,1x,i0," kB")') trim(label), rss
      end if
   end subroutine log_memory_usage

   !> Same as log_memory_usage, but also reports delta to previous sample.
   subroutine log_memory_usage_delta(unit, label, last_rss)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: label
      integer(int64), intent(inout) :: last_rss
      integer(int64) :: rss, delta

      if (.not.memlog_enabled()) return
      rss = rss_kb()
      if (rss < 0_int64) return
      if (last_rss < 0_int64) then
         write(unit,'(1x,"[mem]",1x,a,1x,i0," kB")') trim(label), rss
      else
         delta = rss - last_rss
         write(unit,'(1x,"[mem]",1x,a,1x,i0," kB",1x,"(Δ",i0," kB)")') &
            trim(label), rss, delta
      end if
      last_rss = rss
   end subroutine log_memory_usage_delta

end module xtb_mctc_meminfo
