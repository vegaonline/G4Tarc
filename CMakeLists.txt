#-----------------------------------------------------------------
# $ID: CMakeLists.txt
# Code for TARC
# Abhijit Bhattacharyya
# Nov 2017 : CERN
#----------------------------------------------------------------
cmake_minimum_required(VERSION 2.5 FATAL_ERROR)

project(TARC)
cmake_policy(SET CMP0003 NEW)

option(WITH_GEANT4_UIVIS "Build with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED gdml qt ui_all vis_all)
else()
  find_packages(Geant4 REQUIRED gdml qt)
endif()

#--------------------------------------------------------------
# Find package ROOT (required for histogram)
#--------------------------------------------------------------
find_package(ROOT REQUIRED)
#SET(MYCXXFLAGS "-g -O0 -pg")
#SET(MYCXXLINKFLAGS "-pg")
#SET( CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS} ${MYCXXFLAGS}" )
#SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MYCXXLINKFLAGS}" )

#---------------------------------------------------------------
# Setup include and link directory
#---------------------------------------------------------------
include(${Geant4_USE_FILE})
include_directories( ${PROJECT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR} ${ROOT_INCLUDE_DIR})

#---------------------------------------------------------------
# Locate Headers and sources
#---------------------------------------------------------------
file (GLOB headers
  ${PROJECT_SOURCE_DIR}/include/*.hh
)

file (GLOB sources
  ${PROJECT_SOURCE_DIR}/src/*.cxx
)

#---------------------------------------------------------------
# Add the execuatble and link to libraries
#---------------------------------------------------------------
add_executable(tarcVega tarcVega.cc ${sources} ${headers})
target_link_libraries(tarcVega ${Geant4_LIBRARIES} ${HBOOK_LIBRARIES} ${ROOT_LIBRARIES}) 

#---------------------------------------------------------------
# Copy scripts to build directory
#---------------------------------------------------------------
set(G4TARC_SCRIPTS
  Data/TARC_EXPT_DATA/TARC_EXPTL_DATA.txt
  Data/GeomData/tarc_geometry.gdml 
  Data/GeomData/schema/gdml_core.xsd 
  Data/GeomData/schema/gdml_define.xsd 
  Data/GeomData/schema/gdml_extensions.xsd
  Data/GeomData/schema/gdml_materials.xsd 
  Data/GeomData/schema/gdml_parameterised.xsd 
  Data/GeomData/schema/gdml_replicas.xsd
  Data/GeomData/schema/gdml_solids.xsd 
  Data/GeomData/schema/gdml.xsd
  run.mac 
  vis.mac
  TARCVEGA.run
)
foreach(_script ${G4TARC_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
  )
endforeach()

#---------------------------------------------------------------
# Install the executable to bin directory
#             under CMAKE_INSTALL_PREFIX
#---------------------------------------------------------------
install(TARGETS tarcVega DESTINATION bin)
