#ifndef EXPORT_H
#define EXPORT_H 1

#if defined(_MSC_VER)
#  define BURKHARDT_EXPORT __declspec(dllexport)
#  define BURKHARDT_LOCAL
#elif defined(__GNUC__) && ((__GNUC__ * 100 + __GNUC_MINOR__) >= 480) && (__cplusplus >= 201103)
#  define BURKHARDT_EXPORT [[gnu::visibility("default")]]
#  define BURKHARDT_LOCAL [[gnu::visibility("hidden")]]
#elif defined(__GNUC__) && ((__GNUC__ * 100 + __GNUC_MINOR__) >= 303)
#  define BURKHARDT_EXPORT __attribute__ ((visibility("default")))
#  define BURKHARDT_LOCAL __attribute__ ((visibility("hidden")))
#else
#  define BURKHARDT_EXPORT
#  define BURKHARDT_LOCAL
#endif

#endif
