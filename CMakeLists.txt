#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.6 FATAL_ERROR)
project(MOS_test C CXX Fortran)
enable_language(Fortran)

#execute_process(COMMAND "cd ${PROJECT_SOURCE_DIR}/src/CRY/; make" RESULT_VARIABLE rv)
#message("${PROJECT_SOURCE_DIR}/src/CRY/")
#message("rv='${rv}'")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if (WITH_GEANT4_UIVIS)
    find_package(Geant4 REQUIRED gdml ui_all vis_all)
else ()
    find_package(Geant4 REQUIRED gdml)
endif ()

find_package(PkgConfig REQUIRED)
find_package(Threads REQUIRED)
pkg_search_module(UUID REQUIRED uuid)
#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/src/CRY/src/)
include_directories(${PROJECT_SOURCE_DIR}/build/)
include_directories(${PROJECT_SOURCE_DIR}/src/include
        ${Geant4_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src/src/C_functions/uuid/ ${PROJECT_SOURCE_DIR}/src/src/msis/C/)


file(GLOB sources ${PROJECT_SOURCE_DIR}/src/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/src/include/*.hh)
file(GLOB sourcesf ${PROJECT_SOURCE_DIR}/src/src/*.f)
#file(GLOB sourcesf2 ${PROJECT_SOURCE_DIR}/src/src/msis/*.for)
file(GLOB sourcesMSISC ${PROJECT_SOURCE_DIR}/src/src/msis/C/*.c)
file(GLOB sourcesMSISH ${PROJECT_SOURCE_DIR}/src/src/msis/C/*.h)
file(GLOB sourcesf90 ${PROJECT_SOURCE_DIR}/src/src/parma/*.f90)
file(GLOB sourcesf902 ${PROJECT_SOURCE_DIR}/src/src/unique_id/*.f90)
file(GLOB sourcesc ${PROJECT_SOURCE_DIR}/src/src/*.c)
file(GLOB sourcesc2 ${PROJECT_SOURCE_DIR}/src/src/C_functions/*.c)
file(GLOB sourcesh ${PROJECT_SOURCE_DIR}/src/include/*.h)
file(GLOB main ${PROJECT_SOURCE_DIR}/src/mos_test.cc)

#${hdf}
#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(mos_test ${main} ${sourcesf} ${sourcesMSISH} ${sourcesMSISC} ${sourcesf90} ${sourcesf902} ${sourcesc} ${sourcesc2} ${sourcesh} ${sources} ${headers})
target_link_libraries(mos_test ${Geant4_LIBRARIES})

target_include_directories(mos_test PUBLIC ${UUID_INCLUDE_DIRS})
target_link_libraries(mos_test ${UUID_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -mfpmath=sse -msse2 -ffp-contract=off")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-512")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # using Intel C++
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fp-model=extended")
endif ()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS mos_test DESTINATION bin)

