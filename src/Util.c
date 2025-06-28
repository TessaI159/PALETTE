#include "Util.h"
#include "Color.h"

uint_fast64_t measure_cycles(double (*func)(Color*, Color*), Color* col1,
			     Color* col2) {
	uint64_t start = __rdtsc();
	func(col1, col2);
	uint64_t end = __rdtsc();
	return end - start;
}

uint_fast64_t test_diff_speed(double (*func)(Color*, Color*),
			      const uint64_t RUNS) {
	uint_fast64_t tot = 0;
	for (uint i = 0; i < RUNS; ++i) {
		Color* col1 =
		    Color_create(rand() % 255, rand() % 255, rand() % 255);
		Color* col2 =
		    Color_create(rand() % 255, rand() % 255, rand() % 255);
		tot += measure_cycles(func, col1, col2);
		Color_destroy(col1);
		Color_destroy(col2);
	}
	return tot;
}

void* aligned_malloc(uint64_t alignment, uint64_t size) {
#if defined(_MSC_VER)
	return _aligned_malloc(size, alignment);
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
	return aligned_alloc(alignment, size); // C11
#else
	void* ptr = NULL;
	if (posix_memalign(&ptr, alignment, size) != 0)
		return NULL;
	return ptr;
#endif
}

void aligned_free(void* ptr) {
#if defined(_MSC_VER)
	_aligned_free(ptr);
#else
	free(ptr);
#endif
}



