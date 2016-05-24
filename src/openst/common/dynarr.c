#include "openst/common/dynarr.h"


void *OpenST_DYNARR_Init(struct OpenST_DYNARR *arr, size_t capacity, size_t element_size){
  arr->data = malloc(element_size * capacity);
  if(arr->data != NULL){
    arr->capacity = capacity;
    arr->element_size = element_size;
  } else {
    arr->capacity = 0;
    arr->element_size = 0;
  }
  arr->num = 0;
  return arr->data;
}


void OpenST_DYNARR_Free(struct OpenST_DYNARR *arr){
  free(arr->data);
  arr->capacity = 0;
  arr->element_size = 0;
  arr->num = 0;
}


void *OpenST_DYNARR_Grow(struct OpenST_DYNARR *arr){
  size_t newcapacity = arr->capacity * DYNARR_GROWTH_FACTOR;
  void *newptr = realloc(arr->data, arr->element_size * newcapacity);
  if(newptr != NULL){
    arr->capacity = newcapacity;
    arr->data = newptr;
  }
  return newptr;
}


void *OpenST_DYNARR_Shrink(struct OpenST_DYNARR *arr){
  size_t newcapacity = arr->num;
  void *newptr = realloc(arr->data, arr->element_size * newcapacity);
  if(newptr != NULL){
    arr->capacity = newcapacity;
    arr->data = newptr;
  }
  return newptr;
}


void *OpenST_DYNARR_At(struct OpenST_DYNARR *arr, size_t i){
  return (void *)((char *)arr->data + i * arr->element_size);
}


void *OpenST_DYNARR_Pushback(struct OpenST_DYNARR *arr, void *element){
  if(arr->num == arr->capacity){
    if(OpenST_DYNARR_Grow(arr) == NULL) return NULL;
  }
  void *dest = OpenST_DYNARR_At(arr, arr->num);
  memcpy(dest, element, arr->element_size);
  ++(arr->num);
  return dest;
}
