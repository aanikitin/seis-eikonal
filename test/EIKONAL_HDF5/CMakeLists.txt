cmake_minimum_required(VERSION 3.6.1)

set(EIKONAL_HDF5_SRC src/main.c src/fio_hdf5.c)
add_executable(EIKONAL_HDF5 ${EIKONAL_HDF5_SRC})
    
set(_SAVED_HDF5_DIR ${HDF5_DIR})

# try to find HDF5 supplied module
# CMake built-in FindHDF5 is bogus on Windows - tries to link with dll directly
FIND_PACKAGE(HDF5 "1.8.13" NO_MODULE)

IF(NOT HDF5_FOUND)
    message(STATUS "Fallback on FindHDF5 supplied with CMake")
    
    set(HDF5_DIR ${_SAVED_HDF5_DIR} CACHE PATH "HDF5 install dir" FORCE)
    FIND_PACKAGE(HDF5 "1.8.13" REQUIRED)

    message(STATUS "HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}")
    message(STATUS "HDF5_LIBRARIES: ${HDF5_LIBRARIES}")
    message(STATUS "HDF5_HL_LIBRARIES: ${HDF5_HL_LIBRARIES}")

    target_include_directories(EIKONAL_HDF5 SYSTEM PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(EIKONAL_HDF5 ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
    
ELSE()
    message(STATUS "Found HDF5 ${HDF5_VERSION_STRING}")
    message(STATUS "HDF5_INCLUDE_DIR: ${HDF5_INCLUDE_DIR}")

    find_library(HDF5_LIB hdf5)
    find_library(HDF5_HL_LIB hdf5_hl)
    message(STATUS "HDF5_LIB: ${HDF5_LIB}")
    message(STATUS "HDF5_HL_LIB: ${HDF5_HL_LIB}")

    target_include_directories(EIKONAL_HDF5 SYSTEM PUBLIC ${HDF5_INCLUDE_DIR})
    target_link_libraries(EIKONAL_HDF5 ${HDF5_LIB})
    target_link_libraries(EIKONAL_HDF5 ${HDF5_HL_LIB})
ENDIF()

    target_compile_definitions(EIKONAL_HDF5 PRIVATE H5_BUILT_AS_DYNAMIC_LIB)
    
if(WIN32)
    target_link_libraries(EIKONAL_HDF5 OPENST_SHARED)
    target_compile_definitions(EIKONAL_HDF5 PRIVATE -DOPENST_LINK_SHARED)
else()
    target_link_libraries(EIKONAL_HDF5 OPENST_STATIC)
endif()
