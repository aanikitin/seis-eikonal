# this function removes all files from "${PROJECT_SOURCE_DIR}"/IMP_DIR in SRC_LIST and leaves only IMP_NAME
function(SelectImplementation SRC_LIST IMP_DIR IMP_NAME)
    file(GLOB_RECURSE IMP_LIST RELATIVE "${PROJECT_SOURCE_DIR}" "${IMP_DIR}/*")
    list(REMOVE_ITEM SRC_LIST ${IMP_LIST})
    list(APPEND SRC_LIST "${IMP_DIR}/${IMP_NAME}")
    set(SRC_LIST "${SRC_LIST}" PARENT_SCOPE)
endfunction(SelectImplementation)
