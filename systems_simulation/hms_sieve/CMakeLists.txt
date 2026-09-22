cmake_minimum_required(VERSION 3.16)
project(hms_coordinate_viewer LANGUAGES CXX)

# Disable on an iFarm installation without interactive UI support. File-based
# visualization (e.g. VRML2FILE) still requires Geant4 visualization libraries.
option(HMS_WITH_INTERACTIVE "Build the interactive Geant4 UI" ON)
if(HMS_WITH_INTERACTIVE)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED vis_all)
endif()
include(${Geant4_USE_FILE})

add_executable(hms_geometry_viewer src/hms_geometry_viewer.cc)
target_compile_features(hms_geometry_viewer PRIVATE cxx_std_17)
target_link_libraries(hms_geometry_viewer PRIVATE ${Geant4_LIBRARIES})
target_compile_definitions(hms_geometry_viewer PRIVATE
  HMS_WITH_INTERACTIVE=$<BOOL:${HMS_WITH_INTERACTIVE}>)
foreach(macro view.mac export_vrml.mac)
  configure_file(macros/${macro} ${CMAKE_CURRENT_BINARY_DIR}/${macro} COPYONLY)
endforeach()
