set(BRT3D_EX1_SRC src/main.c)
add_executable(BRT3D_EX1 ${BRT3D_EX1_SRC})
target_link_libraries(BRT3D_EX1 OPENST_STATIC)

add_test(NAME BRT3D_EX1 COMMAND $<TARGET_FILE:BRT3D_EX1>)
