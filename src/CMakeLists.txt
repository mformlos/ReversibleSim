
include_directories ("${PROJECT_SOURCE_DIR}/lib/Particles")
include_directories ("${PROJECT_SOURCE_DIR}/lib/MolecularDynamics")

add_executable(moltest mol_test.cpp)
target_link_libraries(moltest Random Molecule Particle System)

add_executable(LangevinTest Langevin_test.cpp)
target_link_libraries(LangevinTest Random Molecule Particle System)

add_executable(ReversibleSim_auto ReversibleSim_auto.cpp)
target_link_libraries(ReversibleSim_auto Random Molecule Particle System)

add_executable(Widom Widom.cpp)
target_link_libraries(Widom Random System)

add_executable(Widom_parallel Widom_parallel.cpp)
target_link_libraries(Widom_parallel Random Molecule)

add_executable(Radial_dist radial_distribution.cpp)
target_link_libraries(Radial_dist Particle)


add_executable(Radial_dist_all radial_distribution_all.cpp)
target_link_libraries(Radial_dist_all Particle)

add_executable(Structure_factor structure_factor.cpp)
target_link_libraries(Structure_factor Molecule)

add_executable(Form_factor form_factor.cpp)
target_link_libraries(Form_factor System)

add_executable(Create_scattering_vectors create_scattering_vectors.cpp)

add_executable(cluster_calc cluster_calc.cpp)
