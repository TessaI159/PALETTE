#ifndef ALIGN_H
#define ALIGN_H

#if defined(_MSC_VER)
#define ALIGN32 __declspec(align(32))
#else
#define ALIGN32 __attribute__((aligned(32)))
#endif

#if defined(_MSC_VER)
#define ALIGN16 __declspec(align(32))
#else
#define ALIGN16 __attribute__((aligned(16)))
#endif

#include <stddef.h>

#ifdef _WIN32
#include <malloc.h>
static inline void *aligned_malloc(size_t alignment, size_t size) {
	return _aligned_malloc(size, alignment);
}

static inline void aligned_free(void *ptr) {
	_aligned_free(ptr);
}
#else
#include <stdlib.h>
static inline void *aligned_malloc(size_t alignment, size_t size) {
	void *ptr = NULL;
	if (posix_memalign(&ptr, alignment, size) != 0)
		return NULL;
	return ptr;
}

static inline void aligned_free(void *ptr) {
	free(ptr);
}
#endif

#endif
