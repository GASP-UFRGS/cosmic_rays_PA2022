
#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(cosmic_rays_PA2022)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
set(PROJECT_SRC

     src/B1ActionInitialization.cc
     src/B1DetectorConstruction.cc
     src/B1EventAction.cc
     src/B1PrimaryGeneratorAction.cc
     src/B1RunAction.cc
     src/B1SteppingAction.cc
     src/B1TrackingAction.cc
     src/B1SD.cc
     src/B1Hits.cc
)

set(PROJECT_HEADER
include/B1ActionInitialization.hh
include/B1DetectorConstruction.hh
include/B1EventAction.hh
include/B1PrimaryGeneratorAction.hh
include/B1RunAction.hh
include/B1SteppingAction.hh
include/B1Analysis.hh
include/B1TrackingAction.hh
include/B1SD.hh
include/B1Hits.hh
)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(cosmic_rays cosmic_rays.cc ${PROJECT_SRC} ${PROJECT_HEADER})
target_link_libraries(cosmic_rays ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEB1_SCRIPTS
  exampleB1.in
  exampleB1.out
  init_vis.mac
  run1.mac
  run2.mac
  vis.mac
  comandos_shower
  )

foreach(_script ${EXAMPLEB1_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()


#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS cosmic_rays DESTINATION bin)
