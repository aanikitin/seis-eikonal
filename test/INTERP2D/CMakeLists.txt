set(INTERP2D_SRC src/main.c)
add_executable(INTERP2D ${INTERP2D_SRC})
if(WIN32)
    target_link_libraries(INTERP2D OPENST_SHARED)
    target_compile_definitions(INTERP2D PRIVATE -DOPENST_LINK_SHARED)
else()
    target_link_libraries(INTERP2D OPENST_STATIC)
endif()

add_test(NAME INTERP2D COMMAND $<TARGET_FILE:INTERP2D>)
