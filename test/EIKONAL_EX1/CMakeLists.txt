set(EIKONAL_EX1_SRC src/main.c)
add_executable(EIKONAL_EX1 ${EIKONAL_EX1_SRC})
target_link_libraries(EIKONAL_EX1 OPENST_STATIC)
