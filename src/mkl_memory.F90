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
#else
   subroutine release_mkl_buffers()
   end subroutine release_mkl_buffers

   subroutine disable_mkl_fast_mm()
   end subroutine disable_mkl_fast_mm
#endif

end module xtb_mkl_memory
