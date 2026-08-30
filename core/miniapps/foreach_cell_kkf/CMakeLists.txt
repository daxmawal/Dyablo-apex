add_executable(foreach_cell_kkf_replay foreach_cell_kkf_replay.cpp)
target_include_directories(foreach_cell_kkf_replay PRIVATE "${PROJECT_SOURCE_DIR}/src" "${PROJECT_SOURCE_DIR}/src/utils/mpi")
target_link_libraries(foreach_cell_kkf_replay PRIVATE Kokkos::kokkos krepe::kernel_replayer $<$<BOOL:${DYABLO_USE_MPI}>:mpiUtils>)
