# TODO: cleanup build scripts
# TODO: move fortran to CMake build

# Enable LTO-capable GNU wrappers
if(GCC_LTO)
    message(STATUS "GCC_LTO: ${GCC_LTO}")
    set(CMAKE_AR "gcc-ar" CACHE STRING "AR" FORCE)
    set(CMAKE_RANLIB "gcc-ranlib" CACHE STRING "RANLIB" FORCE)
    set(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> rcs <TARGET> <LINK_FLAGS> <OBJECTS>" CACHE STRING "AR TEMPLATE" FORCE)
endif()

project(OpenST)
cmake_minimum_required(VERSION 3.0)
set(CMAKE_VERBOSE_MAKEFILE OFF)
message(STATUS "CMAKE_VERSION: ${CMAKE_VERSION}")

LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(SelectImplementation)
include(SetStdC99)

SetStdC99()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release
        CACHE STRING "Build type: Debug | Release" FORCE)
endif()

set(OPENST_VERSION "0.2.0")
set(OPENST_SOVERSION "0.2")
message(STATUS "OPENST_SOVERSION " ${OPENST_SOVERSION})
message(STATUS "OPENST_VERSION " ${OPENST_VERSION})
set(CMAKE_SKIP_RPATH TRUE)

if(NOT MSVC)
set(OPENST_DEFAULT_FSM3D "fsm3d_bfsm_openmp_v1.c")
set(OPENST_DEFAULT_LSM3D "lsm3d_blsm_openmp_v1.c")
set(OPENST_DEFAULT_FSM3D_BLOCKSERIAL "fsm3d_blockserial_v2.c")
set(OPENST_DEFAULT_LSM3D_BLOCKSERIAL "lsm3d_blockserial_v2.c")
elseif(MSVC)
set(OPENST_DEFAULT_FSM3D "fsm3d_bfsm_openmp_v2.c")
set(OPENST_DEFAULT_LSM3D "lsm3d_blsm_openmp_v2.c")
set(OPENST_DEFAULT_FSM3D_BLOCKSERIAL "fsm3d_blockserial_v2.c")
set(OPENST_DEFAULT_LSM3D_BLOCKSERIAL "lsm3d_blockserial_v2.c")
endif()

if(NOT DEFINED FSM3D_IMP)
    set(FSM3D_IMP "${OPENST_DEFAULT_FSM3D}")
endif()

if(NOT DEFINED LSM3D_IMP)
    set(LSM3D_IMP "${OPENST_DEFAULT_LSM3D}")
endif()

if(NOT DEFINED FSM3D_BLOCKSERIAL_IMP)
    set(FSM3D_BLOCKSERIAL_IMP "${OPENST_DEFAULT_FSM3D_BLOCKSERIAL}")
endif()

if(NOT DEFINED LSM3D_BLOCKSERIAL_IMP)
    set(LSM3D_BLOCKSERIAL_IMP "${OPENST_DEFAULT_LSM3D_BLOCKSERIAL}")
endif()

# set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib")

if(DEFINED OPENST_PATH_LIB_SUFFIX)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/${OPENST_PATH_LIB_SUFFIX}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
    "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${OPENST_PATH_LIB_SUFFIX}")
endif()

# ======================
# ==== CONFIG START ====
# ======================

# COMPILER FLAGS
message(STATUS "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER} \
[${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}]")
message(STATUS "CMAKE_AR: ${CMAKE_AR}")
message(STATUS "CMAKE_RANLIB: ${CMAKE_RANLIB}")

# Enable OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${OpenMP_C_FLAGS} ${CMAKE_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
else()
    message(WARNING "CMake reports no OpenMP support!")
endif()

if(CMAKE_COMPILER_IS_GNUCC)
    # enable LTO flags
    if(GCC_LTO)
        set(CMAKE_C_FLAGS "-flto ${CMAKE_C_FLAGS}")
    endif()
    # warnings
    set(CMAKE_C_FLAGS "-Winline -Wall -Wextra -Wconversion \
-Wfatal-errors -Wdouble-promotion -fmessage-length=0 \
-pedantic ${CMAKE_C_FLAGS}")
    if(${CMAKE_C_COMPILER_VERSION} VERSION_GREATER "4.9.0")
        set(CMAKE_C_FLAGS "-fdiagnostics-color=always ${CMAKE_C_FLAGS}")
    endif()
endif()

message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_C_FLAGS_DEBUG: ${CMAKE_C_FLAGS_DEBUG}")
message(STATUS "CMAKE_C_FLAGS_RELEASE: ${CMAKE_C_FLAGS_RELEASE}")

if(${CMAKE_BUILD_TYPE} MATCHES Debug)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS}")
elseif(${CMAKE_BUILD_TYPE} MATCHES Release)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS_RELEASE} ${CMAKE_C_FLAGS}")
endif()

message(STATUS "CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")

# SET BUILDINFO
# WARNING! can be confusing when using generators for IDEs
add_definitions(-DBUILDINFO_C_COMPILER="${CMAKE_C_COMPILER}")
add_definitions(-DBUILDINFO_C_COMPILER_ID="${CMAKE_C_COMPILER_ID}")
add_definitions(-DBUILDINFO_C_COMPILER_VERSION="${CMAKE_C_COMPILER_VERSION}")
add_definitions(-DBUILDINFO_C_FLAGS="${CMAKE_C_FLAGS}")

# GATHER ALL SOURCES
include_directories(BEFORE include)

# 'evil' globbing. remember to rerun cmake when changing project structure.
aux_source_directory(src SRC_LIST)
file(GLOB_RECURSE SRC_LIST RELATIVE "${PROJECT_SOURCE_DIR}" "src/*.c" "include/*.h")
# message(STATUS "SRC_LIST: ${SRC_LIST}")

# ==============================
# === IMPLEMENTATIONS CHOICE ===
# ==============================
# FSM3D
SelectImplementation("${SRC_LIST}" "src/openst/eikonal/fsm/fsm3d_computepartial" "${FSM3D_IMP}")
SelectImplementation("${SRC_LIST}" "src/openst/eikonal/fsm/fsm3d_blockserial" "${FSM3D_BLOCKSERIAL_IMP}")
# LSM3D
SelectImplementation("${SRC_LIST}" "src/openst/eikonal/lsm/lsm3d_computepartial" "${LSM3D_IMP}")
SelectImplementation("${SRC_LIST}" "src/openst/eikonal/lsm/lsm3d_blockserial" "${LSM3D_BLOCKSERIAL_IMP}")

message(STATUS "FINAL SRC_LIST: ${SRC_LIST}")

# SET OUTPUT PATHS AND NAMES
message(STATUS
    "CMAKE_RUNTIME_OUTPUT_DIRECTORY: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
message(STATUS
    "CMAKE_LIBRARY_OUTPUT_DIRECTORY: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
message(STATUS
    "CMAKE_ARCHIVE_OUTPUT_DIRECTORY: ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}")

# STATIC library target
add_library(OPENST_STATIC STATIC ${SRC_LIST})
set_target_properties(OPENST_STATIC PROPERTIES
    VERSION ${OPENST_VERSION} SOVERSION ${OPENST_SOVERSION})
set_property(TARGET OPENST_STATIC PROPERTY C_STANDARD 99)
if(NOT WIN32)
set_target_properties(OPENST_STATIC PROPERTIES OUTPUT_NAME "openst")
else()
set_target_properties(OPENST_STATIC PROPERTIES
    OUTPUT_NAME "openst-${OPENST_SOVERSION}-static")
endif()

# SHARED library target
add_library(OPENST_SHARED SHARED ${SRC_LIST})
target_compile_definitions(OPENST_SHARED PRIVATE -DOPENST_LINK_SHARED)
set_target_properties(OPENST_SHARED PROPERTIES
    VERSION ${OPENST_VERSION} SOVERSION ${OPENST_SOVERSION})
set_property(TARGET OPENST_SHARED PROPERTY C_STANDARD 99)
if(NOT WIN32)
set_target_properties(OPENST_SHARED PROPERTIES OUTPUT_NAME "openst")
else()
set_target_properties(OPENST_SHARED PROPERTIES
    OUTPUT_NAME "openst-${OPENST_SOVERSION}")
endif()

# ==== LIBRARIES ====
# system libs (math etc.)
if(UNIX)
  find_library(lib_math m)
  target_link_libraries(OPENST_SHARED ${lib_math})
  target_link_libraries(OPENST_STATIC ${lib_math})
endif()

# ==== TESTS ====
enable_testing()
add_subdirectory(test build/test)

# ==== MISC ====
# workaround for Qt Creator project tree not showing includes etc.
# doesn't affect build process
file(GLOB_RECURSE DUMMY_HEADERS "${PROJECT_SOURCE_DIR}" *.h *.c)
add_custom_target(dummy ALL SOURCES ${DUMMY_HEADERS})
