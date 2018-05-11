add_library(Particle Particle.cpp)

add_library(Molecule Molecule.cpp)
target_link_libraries(Molecule Particle)
