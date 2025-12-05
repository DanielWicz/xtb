# MPI support sources
set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/numhess_mpi.F90"
)

set(srcs ${srcs} PARENT_SCOPE)
