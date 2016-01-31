project(FSM3DOPENMP)
cmake_minimum_required(VERSION 2.8)
# set(CMAKE_VERBOSE_MAKEFILE OFF)

# set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
message(STATUS "CMAKE_RUNTIME_OUTPUT_DIRECTORY: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

# ======================
# ==== CONFIG START ====
# ======================

# ==== COMPILER FLAGS ====
message(STATUS "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER} [${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}]")

if(UNIX)
	if(NOT DEFINED C_FLAGS_OPT)
	  set(C_FLAGS_OPT "-O3" CACHE STRING "Compile optimization flags")
	endif()
	set(CMAKE_C_FLAGS "${C_FLAGS_OPT} -Wall -Wextra -fmessage-length=0 \
	-std=c99 -pedantic -fopenmp ${CMAKE_C_FLAGS}")
	if(CMAKE_COMPILER_IS_GNUCC)
		if(${CMAKE_C_COMPILER_VERSION} VERSION_GREATER "4.9.0")
			set(CMAKE_C_FLAGS "-fdiagnostics-color=always ${CMAKE_C_FLAGS}")
		endif()
	endif()
elseif(WIN32)
	if(NOT DEFINED C_FLAGS_OPT)
	  set(C_FLAGS_OPT "/O2" CACHE STRING "Compile optimization flags")
	endif()
	set(CMAKE_C_FLAGS "${C_FLAGS_OPT} /openmp")
endif()

message(STATUS "CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")

# ==== SOURCE CODE CONFIG ====
# set source directory
include_directories(src)

# set FSM implementation
if(NOT DEFINED IMP_FSM3D)
  set(IMP_FSM3D_PATH fsm_openmp_v3.c CACHE FILEPATH "FSM3D implementation .c filepath")
  get_filename_component(IMP_FSM3D "${IMP_FSM3D_PATH}" NAME)
endif()
add_definitions(-DIMP_FSM3D="${IMP_FSM3D}")
message(STATUS "IMP_FSM3D: ${IMP_FSM3D}")

# set source directory
aux_source_directory(src SRC_LIST)
aux_source_directory(src/fsm_stable IMP_FSM3DSTABLE)
list(APPEND SRC_LIST src/fsm/${IMP_FSM3D})
list(APPEND SRC_LIST ${IMP_FSM3DSTABLE})

# ==== TARGET CONFIG ====
# resulting binary file name
set(BIN_NAME "${PROJECT_NAME}_${IMP_FSM3D}")
add_definitions(-DBIN_NAME="${BIN_NAME}")
message(STATUS "BIN_NAME: ${BIN_NAME}")
# set output executable target
add_executable(${BIN_NAME} ${SRC_LIST} ${HEADER_LIST})

# ==== LIBRARIES ====
# system libs (math etc.)
if(UNIX)
find_library(lib_math m)
target_link_libraries(${BIN_NAME} ${lib_math})
endif()

# ==== BUILD_INFO definitions ====
# useful for compiler/optimization flags tracking during performance testing
add_definitions(-DBUILD_INFO_CMAKE_C_COMPILER="${CMAKE_C_COMPILER}")
add_definitions(-DBUILD_INFO_CMAKE_C_COMPILER_ID="${CMAKE_C_COMPILER_ID}")
add_definitions(-DBUILD_INFO_CMAKE_C_COMPILER_VERSION="${CMAKE_C_COMPILER_VERSION}")
add_definitions(-DBUILD_INFO_CMAKE_C_FLAGS="${CMAKE_C_FLAGS}")

# ==== MISC ====
# workaround for Qt Creator project tree not showing includes etc.
# doesn't affect build process
file(GLOB_RECURSE DUMMY_HEADERS src *.h src/fsm *.c src/fsm *.h)
add_custom_target(dummy ALL SOURCES ${DUMMY_HEADERS})
