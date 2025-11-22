! This file is part of xtb.
!
! Helper utilities to keep Intel MKL's Fast Memory Manager from
! accumulating cached buffers between geometry iterations.  These calls
! are no-ops when MKL is not available.
module xtb_mkl_memory
   use iso_c_binding, only : c_int
   implicit none
   private

   public :: release_mkl_buffers
   public :: disable_mkl_fast_mm
   public :: init_mkl_memory

contains

#ifdef WITH_MKL
   interface
      subroutine mkl_free_buffers() bind(C, name='mkl_free_buffers')
      end subroutine mkl_free_buffers

      subroutine mkl_thread_free_buffers() bind(C, name='mkl_thread_free_buffers')
      end subroutine mkl_thread_free_buffers

      integer(c_int) function mkl_disable_fast_mm() bind(C, name='mkl_disable_fast_mm')
      end function mkl_disable_fast_mm
   end interface

   !> Release MKL per-thread caches as well as the global fast memory pool.
   subroutine release_mkl_buffers()
#ifdef _OPENMP
      !$omp parallel
         call mkl_thread_free_buffers()
      !$omp end parallel
#endif
      call mkl_free_buffers()
   end subroutine release_mkl_buffers

   !> Programmatic equivalent of setting MKL_DISABLE_FAST_MM=1.
   subroutine disable_mkl_fast_mm()
      integer(c_int) :: ierr
      ierr = mkl_disable_fast_mm()
   end subroutine disable_mkl_fast_mm

   !> One-time initialization for MKL memory handling.
   !! On Intel oneAPI builds the MKL fast memory manager is prone to retain
   !! large buffers across geometry steps when OpenMP is used.  To prevent the
   !! observed monotonic RSS growth in GFN2 optimizations, disable the fast
   !! memory manager by default on Intel compilers unless the user explicitly
   !! re-enables it via `XTB_ENABLE_MKL_FAST_MM`.
   subroutine init_mkl_memory()
      character(len=32) :: env
      integer :: len, stat
      logical :: disable

      disable = .false.

#if defined(__INTEL_LLVM_COMPILER) || defined(__INTEL_COMPILER)
      ! Default to disabling on Intel oneAPI (ifx) or classic (ifort) builds
      ! where MKL FastMM tends to retain per-thread buffers across geometry
      ! steps under OpenMP.
      disable = .true.
#endif

      ! Allow opt-out: if the user sets XTB_ENABLE_MKL_FAST_MM to a true-ish
      ! value, keep the fast memory manager enabled.
      call get_environment_variable('XTB_ENABLE_MKL_FAST_MM', value=env, length=len, status=stat)
      if (stat == 0 .and. len > 0) then
         if (env(1:1) == '1' .or. env(1:1) == 'T' .or. env(1:1) == 't' .or. env(1:1) == 'Y' .or. env(1:1) == 'y') then
            disable = .false.
         end if
      end if

      ! Allow explicit disable toggle as well
      call get_environment_variable('XTB_DISABLE_MKL_FAST_MM', value=env, length=len, status=stat)
      if (stat == 0 .and. len > 0) then
         if (env(1:1) == '1' .or. env(1:1) == 'T' .or. env(1:1) == 't' .or. env(1:1) == 'Y' .or. env(1:1) == 'y') then
            disable = .true.
         end if
      end if

      if (disable) call disable_mkl_fast_mm()
   end subroutine init_mkl_memory
#else
   subroutine release_mkl_buffers()
   end subroutine release_mkl_buffers

   subroutine disable_mkl_fast_mm()
   end subroutine disable_mkl_fast_mm

   subroutine init_mkl_memory()
   end subroutine init_mkl_memory
#endif

end module xtb_mkl_memory
