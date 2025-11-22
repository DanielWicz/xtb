! This file is part of xtb.
!
! Runtime helpers to harmonize BLAS threading behaviour across vendors
! and to release vendor-managed buffers on shutdown.  OpenBLAS ignores
! changes to OMP_NUM_THREADS/omp_set_num_threads after initialization,
! so we call openblas_set_num_threads explicitly.  MKL caches buffers
! per thread and globally; we provide wrappers to flush them at exit.

module xtb_blas_runtime
   use, intrinsic :: iso_c_binding, only : c_int, c_size_t
   implicit none
   private

   public :: set_blas_threads
   public :: blas_thread_cleanup
   public :: blas_global_cleanup
   public :: blas_flush_buffers
   public :: blas_trim_os
   public :: blas_last_trim_result
   public :: blas_tune_malloc

#ifdef WITH_OPENBLAS
   interface
      subroutine openblas_set_num_threads(num_threads) bind(C, name='openblas_set_num_threads')
         import :: c_int
         integer(c_int), value :: num_threads
      end subroutine openblas_set_num_threads
   end interface
#endif

#ifdef WITH_MKL
   interface
      subroutine mkl_set_num_threads(num_threads) bind(C, name='mkl_set_num_threads')
         import :: c_int
         integer(c_int), value :: num_threads
      end subroutine mkl_set_num_threads

      subroutine mkl_free_buffers() bind(C, name='mkl_free_buffers')
      end subroutine mkl_free_buffers

      subroutine mkl_thread_free_buffers() bind(C, name='mkl_thread_free_buffers')
      end subroutine mkl_thread_free_buffers
   end interface
#endif

   interface
      function malloc_trim(pad) bind(C, name='malloc_trim')
         import :: c_int, c_size_t
         integer(c_int) :: malloc_trim
         integer(c_size_t), value :: pad
      end function malloc_trim

      function mallopt(param, value) bind(C, name='mallopt')
         import :: c_int
         integer(c_int) :: mallopt
         integer(c_int), value :: param, value
      end function mallopt
   end interface

   integer(c_int) :: last_trim = -1
   logical :: tuned = .false.

   integer(c_int), parameter :: M_TRIM_THRESHOLD = -1
   integer(c_int), parameter :: M_MMAP_THRESHOLD = -3
   integer(c_int), parameter :: M_ARENA_MAX      = -8

contains

   subroutine set_blas_threads(num_threads)
      integer, intent(in) :: num_threads
#if defined(WITH_OPENBLAS)
      call openblas_set_num_threads(num_threads)
#endif
#if defined(WITH_MKL)
      call mkl_set_num_threads(num_threads)
#endif
#if !defined(WITH_OPENBLAS) && !defined(WITH_MKL)
      ! Prevent unused variable warnings when neither backend exposes runtime control
      if (num_threads /= 0) continue
#endif
   end subroutine set_blas_threads


   subroutine blas_thread_cleanup()
#ifdef WITH_MKL
      call mkl_thread_free_buffers()
#endif
   end subroutine blas_thread_cleanup


   subroutine blas_global_cleanup()
#ifdef WITH_MKL
      call mkl_free_buffers()
#endif
   end subroutine blas_global_cleanup


   subroutine blas_flush_buffers()
#ifdef WITH_MKL
!$omp parallel
      call blas_thread_cleanup()
!$omp end parallel
      call blas_global_cleanup()
#endif
      call blas_trim_os()
   end subroutine blas_flush_buffers


   subroutine blas_trim_os()
      integer(c_int) :: ret
      if (.not.tuned) call blas_tune_malloc()
      ! Request libc to return free arenas to the OS (pad=0 => trim all)
      ret = malloc_trim(0_c_size_t)
      last_trim = ret
      if (ret == 0) continue ! avoid unused warning when assertions off
   end subroutine blas_trim_os


   function blas_last_trim_result() result(val)
      integer(c_int) :: val
      val = last_trim
   end function blas_last_trim_result


   subroutine blas_tune_malloc()
      integer(c_int) :: rc
      if (tuned) return
      rc = mallopt(M_ARENA_MAX, 2_c_int)
      rc = mallopt(M_TRIM_THRESHOLD, 0_c_int)
      rc = mallopt(M_MMAP_THRESHOLD, 128*1024_c_int)
      tuned = .true.
   end subroutine blas_tune_malloc

end module xtb_blas_runtime
