! This file is part of xtb.
!
! Runtime helpers to harmonize BLAS threading behaviour across vendors
! and to release vendor-managed buffers on shutdown.  OpenBLAS ignores
! changes to OMP_NUM_THREADS/omp_set_num_threads after initialization,
! so we call openblas_set_num_threads explicitly.  MKL caches buffers
! per thread and globally; we provide wrappers to flush them at exit.

module xtb_blas_runtime
   use, intrinsic :: iso_c_binding, only : c_int
   implicit none
   private

   public :: set_blas_threads
   public :: blas_thread_cleanup
   public :: blas_global_cleanup

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

end module xtb_blas_runtime
