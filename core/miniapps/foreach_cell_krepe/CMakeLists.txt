foreach(replay_target IN ITEMS foreach_cell_krepe_replay foreach_cell_krepe_parallel_for_replay)
  add_executable(${replay_target} ${replay_target}.cpp)
  set_target_properties(${replay_target} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
  target_include_directories(${replay_target} PRIVATE "${PROJECT_SOURCE_DIR}/src/utils/mpi")
  target_compile_definitions(${replay_target} PRIVATE DYABLO_USE_KREPE_REPLAYER)
  target_link_libraries(${replay_target} PRIVATE dyablo_core krepe::kernel_replayer)
endforeach()
