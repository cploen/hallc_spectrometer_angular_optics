cmake_minimum_required(VERSION 3.16)
project(hms_sieve_scattering LANGUAGES CXX)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(hms_sieve_scattering src/hms_sieve_scattering.cc)
target_compile_features(hms_sieve_scattering PRIVATE cxx_std_17)
target_link_libraries(hms_sieve_scattering PRIVATE ${Geant4_LIBRARIES})
