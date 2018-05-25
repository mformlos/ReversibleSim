
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