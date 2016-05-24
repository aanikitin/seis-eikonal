/*!
 * \file dynarr.h Simple type-agnostic dynamic array implementation using standard C library.
 * \author Alexandr Nikitin
 */

#ifndef OPENST_COMMON_DYNARR_H
#define OPENST_COMMON_DYNARR_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <string.h>
#include "openst/common/macros.h"

/*!
 * Capacity multiplication factor used when growing array.
 */
#define DYNARR_GROWTH_FACTOR 2u

/*!
 * Dynamic array structure.
 */
struct OpenST_DYNARR{
    void *data; /*!< beginning of memory buffer for storing elements. */
    size_t num; /*!< current number of elements in array. */
    size_t capacity; /*!< current size of memory buffer in elements. */
    size_t element_size; /*!< element size in bytes. */
};

/*!
 * \brief Initialize dynamic array.
 * \param[in,out] arr pointer to dynamic array structure.
 * \param[in] capacity initial capacity.
 * \param[in] element_size array element size.
 * \return pointer to allocated memory buffer or NULL on error.
 */
OPENST_API void *OpenST_DYNARR_Init(struct OpenST_DYNARR *arr, size_t capacity,
                            size_t element_size);

/*!
 * \brief Free data memory buffer of dynamic array.
 * \param[in,out] arr pointer to dynamic array structure.
 */
OPENST_API void OpenST_DYNARR_Free(struct OpenST_DYNARR *arr);

/*!
 * \brief Grow dynamic array.
 * \note Previous memory buffer is left untoched on error (realloc behavior).
 * \param[in,out] arr pointer to dynamic array structure.
 * \return pointer to reallocated memory buffer or NULL on error.
 */
OPENST_API void *OpenST_DYNARR_Grow(struct OpenST_DYNARR *arr);

/*!
 * \brief Shrink dynamic array.
 * \note Previous memory buffer is left untoched on error (realloc behavior).
 * \param[in,out] arr pointer to dynamic array structure.
 * \return pointer to reallocated memory buffer or NULL on error.
 */
OPENST_API void *OpenST_DYNARR_Shrink(struct OpenST_DYNARR *arr);

/*!
 * \brief Get dynamic array element.
 * \note Warning: no memory boundary checks are made.
 * \param[in,out] arr pointer to dynamic array structure.
 * \param[in] i index of element.
 * \return pointer to array element.
 */
OPENST_API void *OpenST_DYNARR_At(struct OpenST_DYNARR *arr, size_t i);

/*!
 * \brief Push element to the end of dynamic array.
 * \note Automatically grows array if capacity is reached.
 * \param[in,out] arr pointer to dynamic array structure.
 * \param[in] element pointer to element to be copied.
 * \return pointer to added array element or NULL on error.
 */
OPENST_API void *OpenST_DYNARR_Pushback(struct OpenST_DYNARR *arr, void *element);

#ifdef __cplusplus
}
#endif

#endif /* OPENST_COMMON_DYNARR_H */
