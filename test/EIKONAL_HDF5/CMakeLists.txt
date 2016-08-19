cmake_minimum_required(VERSION 3.6.1)

find_package(HDF5)

if(HDF5_FOUND)
    set(EIKONAL_HDF5_SRC src/main.c src/fio_hdf5.c)
    add_executable(EIKONAL_HDF5 ${EIKONAL_HDF5_SRC})

    message(STATUS "HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}")
    message(STATUS "HDF5_LIBRARIES: ${HDF5_LIBRARIES}")
    message(STATUS "HDF5_HL_LIBRARIES: ${HDF5_HL_LIBRARIES}")

    target_include_directories(EIKONAL_HDF5 SYSTEM PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(EIKONAL_HDF5 ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

    if(WIN32)
        target_link_libraries(EIKONAL_HDF5 OPENST_SHARED)
        target_compile_definitions(EIKONAL_HDF5 PRIVATE -DOPENST_LINK_SHARED)
    else()
        target_link_libraries(EIKONAL_HDF5 OPENST_STATIC)
    endif()    
else()
    message(WARNING "HDF5 library not found. EIKONAL_HDF5 is disabled.")
endif()
