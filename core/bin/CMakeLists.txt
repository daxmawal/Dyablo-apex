set(EXE_NAME "dyablo")
add_executable(${EXE_NAME}
  ${CMAKE_SOURCE_DIR}/src/dyablo_main.cpp
  ${BACKWARD_ENABLE}
  )

target_include_directories(${EXE_NAME}
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_BINARY_DIR}/src
  )

target_link_libraries(${EXE_NAME}
  dyablo_core
  dyablo_plugins
  )

if (DYABLO_USE_MPI)  
  target_link_libraries(${EXE_NAME} mpiUtils)
endif()

add_backward(${EXE_NAME})

# Files for gitlab-ci testcases
configure_file(blast.ini blast.ini COPYONLY)
# Blast 3D block based 
configure_file(test_blast_3D_block.ini test_blast_3D_block.ini COPYONLY)
configure_file(visu_blast_3D_block.pvsm visu_blast_3D_block.pvsm COPYONLY)
# Blast 3D cell based 
configure_file(test_blast_3D.ini test_blast_3D.ini COPYONLY)
configure_file(visu_blast_3D_cell.pvsm visu_blast_3D_cell.pvsm COPYONLY)
# Blast 2D block based
configure_file(test_blast_2D_block.ini test_blast_2D_block.ini COPYONLY)
configure_file(visu_blast_2D_block.pvsm visu_blast_2D_block.pvsm COPYONLY)
# Blast 2D cell based 
configure_file(test_blast_2D.ini test_blast_2D.ini COPYONLY)
configure_file(visu_blast_2D_cell.pvsm visu_blast_2D_cell.pvsm COPYONLY)
# Riemann 2D cell based 
configure_file(test_riemann_2D.ini test_riemann_2D.ini COPYONLY)
configure_file(visu_riemann_2D_cell.pvsm visu_riemann_2D_cell.pvsm COPYONLY)
# Riemann 2D block based 
configure_file(test_riemann_2D_block.ini test_riemann_2D_block.ini COPYONLY)
configure_file(visu_riemann_2D_block.pvsm visu_riemann_2D_block.pvsm COPYONLY)
# Sod test
configure_file(test_sod_3D.ini  test_sod_3D.ini )
configure_file(visu_sod_3D.pvsm  visu_sod_3D.pvsm )

configure_file(test_sod_2D.ini  test_sod_2D.ini )
configure_file(validate_sod.py validate_sod.py )

# Brio-Wu MHD test
configure_file(test_brio_wu.ini test_brio_wu.ini COPYONLY)
configure_file(validate_brio_wu.py validate_brio_wu.py COPYONLY)
configure_file(brio_wu_reference.npy brio_wu_reference.npy COPYONLY)

# 2D Blast MHD test
configure_file(test_mhd_blast_2D_block.ini test_mhd_blast_2D_block.ini COPYONLY)
configure_file(visu_mhd_blast_2D_block.pvsm visu_mhd_blast_2D_block.pvsm COPYONLY)

# Orszag-Tang Vortex test
configure_file(test_orszag_tang_vortex.ini test_orszag_tang_vortex.ini COPYONLY)


# Gravity 3D test
configure_file(test_gravity_spheres_3D.ini test_gravity_spheres_3D.ini COPYONLY)
configure_file(visu_gravity_spheres_3D.pvsm visu_gravity_spheres_3D.pvsm COPYONLY)

# Zeldovitch 3D test
configure_file(../../settings/cosmo/zeldovitch_dyablo.ini test_zeldovitch_dyablo.ini COPYONLY)
configure_file(../../settings/cosmo/zeldovitch_particles_dyablo.ini test_zeldovitch_particles_dyablo.ini COPYONLY)
configure_file(../../settings/cosmo/visu_zeldovitch.pvsm visu_zeldovitch.pvsm COPYONLY)

# Zeldovitch grafic test
configure_file(../../settings/cosmo/zeldovitch_grafic.ini test_zeldovitch_grafic.ini COPYONLY)

# Kelvin-Helmholtz cell-based
configure_file(test_kelvin_helmholtz.ini test_kelvin_helmholtz.ini COPYONLY)
configure_file(visu_kelvin_helmholtz.pvsm visu_kelvin_helmholtz.pvsm COPYONLY)

# Kelvin-Helmholtz block-based
configure_file(test_kelvin_helmholtz_block.ini test_kelvin_helmholtz_block.ini COPYONLY)
configure_file(visu_kelvin_helmholtz_block.pvsm visu_kelvin_helmholtz_block.pvsm COPYONLY)
# RT test
configure_file(../../settings/cosmo/beam.ini beam.ini COPYONLY)

# Heat conduction
configure_file(test_heat_conduction.ini test_heat_conduction.ini COPYONLY)
configure_file(validate_heat_conduction.py validate_heat_conduction.py COPYONLY)
