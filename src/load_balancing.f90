! This file is part of xtb.
!
! Copyright (C) 2025
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

module xtb_load_balancing
   use omp_lib
   implicit none
   private

   public :: get_load_balance

contains

   !> Heuristic to determine optimal task and thread counts
   !> n_tasks: Number of concurrent singlepoint calculations (outer loop threads)
   !> n_threads_per_task: Number of threads for each singlepoint calculation (inner loop)
   subroutine get_load_balance(natoms, n_tasks, n_threads_per_task)
      integer, intent(in) :: natoms
      integer, intent(out) :: n_tasks
      integer, intent(out) :: n_threads_per_task

      integer :: total_threads
      integer :: optimal_threads_per_scf

      ! Get total available threads from OpenMP environment
      ! Use parallel region to get the max threads available for this team
      total_threads = omp_get_max_threads()

      ! Heuristic:
      ! Large systems consume significant memory per thread in SCF.
      ! To avoid OOM, we must restrict the number of concurrent SCFs.
      ! 400 atoms is the empirical threshold where OOM starts appearing with default settings.
      if (natoms > 400) then
         ! Memory-bound regime
         ! Run strictly serial displacement loop (1 task)
         ! Give all resources to the single SCF task
         n_tasks = 1
         n_threads_per_task = total_threads
      else
         ! Throughput-bound regime
         ! SCF efficiency usually saturates. 
         ! For small/medium systems, 2-4 threads is often optimal efficiency.
         ! Let's say 4 threads is a good target for efficiency.
         
         if (total_threads >= 4) then
            optimal_threads_per_scf = 4
         else
            optimal_threads_per_scf = 1
         end if

         n_tasks = max(1, total_threads / optimal_threads_per_scf)
         n_threads_per_task = max(1, total_threads / n_tasks)
         
         ! Ensure we don't oversubscribe
         if (n_tasks * n_threads_per_task > total_threads) then
             n_tasks = max(1, total_threads / n_threads_per_task)
         end if
      end if

   end subroutine get_load_balance

end module xtb_load_balancing
