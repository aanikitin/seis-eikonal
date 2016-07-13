project(OpenST)
set(OPENST_VERSION "0.1.0")
set(OPENST_SOVERSION "0.1")
message(STATUS "OPENST_SOVERSION " ${OPENST_SOVERSION})
message(STATUS "OPENST_VERSION " ${OPENST_VERSION})
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_SKIP_RPATH TRUE)

if(NOT MSVC)
set(OPENST_DEFAULT_FSM3D "fsm3d_bfsm_openmp_v1.c")
set(OPENST_DEFAULT_LSM3D "lsm3d_blsm_openmp_v1.c")
elseif(MSVC)
set(OPENST_DEFAULT_FSM3D "fsm3d_bfsm_openmp_v2.c")
set(OPENST_DEFAULT_LSM3D "lsm3d_blsm_openmp_v2.c")
endif()

cmake_minimum_required(VERSION 2.8)
# set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_VERBOSE_MAKEFILE OFF)

# set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib")

if(DEFINED OPENST_PATH_LIB_SUFFIX)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/${OPENST_PATH_LIB_SUFFIX}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${OPENST_PATH_LIB_SUFFIX}")
endif()

# ======================
# ==== CONFIG START ====
# ======================

# COMPILER FLAGS
message(STATUS "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER} [${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}]")

if(NOT MSVC)
	if(NOT DEFINED C_FLAGS_OPT)
          set(C_FLAGS_OPT "-O3" CACHE STRING "Compiler optimization flags")
	endif()
	set(CMAKE_C_FLAGS "${C_FLAGS_OPT} -Winline -Wall -Wextra -fmessage-length=0 -std=c99 -pedantic -fopenmp ${CMAKE_C_FLAGS}")
	if(CMAKE_COMPILER_IS_GNUCC)
		if(${CMAKE_C_COMPILER_VERSION} VERSION_GREATER "4.9.0")
			set(CMAKE_C_FLAGS "-fdiagnostics-color=always ${CMAKE_C_FLAGS}")
		endif()
	endif()
elseif(MSVC)
	if(NOT DEFINED C_FLAGS_OPT)
	  set(C_FLAGS_OPT "/O2" CACHE STRING "Compiler optimization flags")
	endif()
	set(CMAKE_C_FLAGS "${C_FLAGS_OPT} /openmp")
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
file(GLOB_RECURSE SRC_LIST RELATIVE "${PROJECT_SOURCE_DIR}" "src/*.c")
# message(STATUS "SRC_LIST: ${SRC_LIST}")

# ==============================
# === IMPLEMENTATIONS CHOICE ===
# ==============================
# TODO: replace with function()

# FSM3D
# remove FSM3D implementation
file(GLOB_RECURSE FSM3D_IMP_LIST RELATIVE "${PROJECT_SOURCE_DIR}" "src/openst/eikonal/fsm/fsm3d_imp/*.c")
message(STATUS "FSM3D_IMP_LIST: ${FSM3D_IMP_LIST}")
list(REMOVE_ITEM SRC_LIST ${FSM3D_IMP_LIST})
# set FSM3D implementation
if(NOT DEFINED FSM3D_IMP)
  set(FSM3D_IMP_PATH ${OPENST_DEFAULT_FSM3D} CACHE FILEPATH "FSM3D implementation .c filepath")
  get_filename_component(FSM3D_IMP "${FSM3D_IMP_PATH}" NAME)
endif()
add_definitions(-DFSM3D_IMP="${FSM3D_IMP}")
message(STATUS "FSM3D_IMP: ${FSM3D_IMP}")
list(APPEND SRC_LIST "src/openst/eikonal/fsm/fsm3d_imp/${FSM3D_IMP}")

# LSM3D
# remove LSM3D implementation
file(GLOB_RECURSE LSM3D_IMP_LIST RELATIVE "${PROJECT_SOURCE_DIR}" "src/openst/eikonal/lsm/lsm3d_imp/*.c")
message(STATUS "LSM3D_IMP_LIST: ${LSM3D_IMP_LIST}")
list(REMOVE_ITEM SRC_LIST ${LSM3D_IMP_LIST})
# set LSM3D implementation
if(NOT DEFINED LSM3D_IMP)
  set(LSM3D_IMP_PATH ${OPENST_DEFAULT_LSM3D} CACHE FILEPATH "LSM3D implementation .c filepath")
  get_filename_component(LSM3D_IMP "${LSM3D_IMP_PATH}" NAME)
endif()
add_definitions(-DLSM3D_IMP="${LSM3D_IMP}")
message(STATUS "LSM3D_IMP: ${LSM3D_IMP}")
list(APPEND SRC_LIST "src/openst/eikonal/lsm/lsm3d_imp/${LSM3D_IMP}")

message(STATUS "FINAL SRC_LIST: ${SRC_LIST}")

# SET OUTPUT PATHS AND NAMES

message(STATUS "CMAKE_RUNTIME_OUTPUT_DIRECTORY: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
message(STATUS "CMAKE_LIBRARY_OUTPUT_DIRECTORY: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
message(STATUS "CMAKE_ARCHIVE_OUTPUT_DIRECTORY: ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}")

# STATIC library target
add_library(OPENST_STATIC STATIC ${SRC_LIST})
set_target_properties(OPENST_STATIC PROPERTIES VERSION ${OPENST_VERSION} SOVERSION ${OPENST_SOVERSION})
if(NOT WIN32)
set_target_properties(OPENST_STATIC PROPERTIES OUTPUT_NAME "openst")
else()
set_target_properties(OPENST_STATIC PROPERTIES OUTPUT_NAME "openst-${OPENST_SOVERSION}-static")
endif()

# SHARED library target
add_library(OPENST_SHARED SHARED ${SRC_LIST})
target_compile_definitions(OPENST_SHARED PRIVATE -DOPENST_LINK_SHARED)
set_target_properties(OPENST_SHARED PROPERTIES VERSION ${OPENST_VERSION} SOVERSION ${OPENST_SOVERSION})
if(NOT WIN32)
set_target_properties(OPENST_SHARED PROPERTIES OUTPUT_NAME "openst")
else()
set_target_properties(OPENST_SHARED PROPERTIES OUTPUT_NAME "openst-${OPENST_SOVERSION}")
endif()

# ==== LIBRARIES ====
# system libs (math etc.)
if(UNIX)
  find_library(lib_math m)
  target_link_libraries(OPENST_SHARED ${lib_math})
  target_link_libraries(OPENST_STATIC ${lib_math})
endif()

# ==== TESTS ====
add_subdirectory(test build/test)

# ==== MISC ====
# workaround for Qt Creator project tree not showing includes etc.
# doesn't affect build process
file(GLOB_RECURSE DUMMY_HEADERS "${PROJECT_SOURCE_DIR}" *.h *.c)
add_custom_target(dummy ALL SOURCES ${DUMMY_HEADERS})
