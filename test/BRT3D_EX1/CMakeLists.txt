set(BRT3D_EX1_SRC src/main.c)
add_executable(BRT3D_EX1 ${BRT3D_EX1_SRC})
if(WIN32)
    target_link_libraries(BRT3D_EX1 OPENST_SHARED)
	target_compile_definitions(BRT3D_EX1 PRIVATE -DOPENST_LINK_SHARED)
else()
    target_link_libraries(BRT3D_EX1 OPENST_STATIC)
endif()

add_test(NAME BRT3D_EX1 COMMAND $<TARGET_FILE:BRT3D_EX1>)
