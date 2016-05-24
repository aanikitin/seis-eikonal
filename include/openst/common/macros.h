#ifndef OPENST_COMMON_MACROS_H
#define OPENST_COMMON_MACROS_H

#define OPENST_M_XSTR(s) OPENST_M_STR(s)
#define OPENST_M_STR(s) #s

#if defined(_WIN32) && defined(OPENST_LINK_SHARED)
#define OPENST_DLLEXPORT __declspec(dllexport)
#define OPENST_DLLIMPORT __declspec(dllimport)
#else
# define OPENST_DLLEXPORT
# define OPENST_DLLIMPORT
#endif

#if defined(OPENST_SHARED_EXPORTS)
#define OPENST_DECLSPEC OPENST_DLLEXPORT
#else
#define OPENST_DECLSPEC OPENST_DLLIMPORT
#endif

#define OPENST_API OPENST_DECLSPEC

#endif /* OPENST_COMMON_MACROS_H */
