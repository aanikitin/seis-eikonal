#ifndef OPENST_COMMON_ERROR_H
#define OPENST_COMMON_ERROR_H

#ifdef __cplusplus
extern "C" {
#endif

typedef enum OPENST_ERR_enum{
    OPENST_ERR_SUCCESS = 0,
    OPENST_ERR_UNDEFINED,
    OPENST_ERR_MEMORY,
    OPENST_ERR_PARAM_INVALID
} OPENST_ERR;

#ifdef __cplusplus
}
#endif

#endif /* OPENST_COMMON_ERROR_H */
