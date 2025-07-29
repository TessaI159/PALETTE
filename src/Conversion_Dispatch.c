#include <stdbool.h>
#include <stdlib.h>

#include "Color.h"
#include "Conversion.h"
#include "Conversion_Internal.h"
#include "FeatureDetection.h"
#include "Parameters.h"

#define INST_SETS 3
#define SPACES 3

static void (*srgb_to_cielab_intern)(struct Color *restrict);
static void (*srgb_to_oklab_intern)(struct Color *restrict);
static void (*oklab_to_srgb_intern)(struct Color *restrict);
static void (*cielab_to_srgb_intern)(struct Color *restrict);

static void (*conversion_table[INST_SETS][SPACES][SPACES])(
    struct Color *restrict) = {
    [AVX] = {[OK] =
		 {[OK] = NULL, [CIE] = NULL, [SRGB] = _mm256_oklab_to_srgb_ps},
	     [CIE] =
		 {[OK] = NULL, [CIE] = NULL, [SRGB] = _mm256_cielab_to_srgb_ps},
	     [SRGB] = {[OK]   = _mm256_srgb_to_oklab_ps,
		       [CIE]  = _mm256_srgb_to_cielab_ps,
		       [SRGB] = NULL}},
    [SSE]  = {[OK]  = {[OK] = NULL, [CIE] = NULL, [SRGB] = _mm_oklab_to_srgb_ps},
	     [CIE] = {[OK] = NULL, [CIE] = NULL, [SRGB] = _mm_cielab_to_srgb_ps},
	     [SRGB] = {[OK]   = _mm_srgb_to_oklab_ps,
		       [CIE]  = _mm_srgb_to_cielab_ps,
		       [SRGB] = NULL}},
    [FB]  = {[OK]   = {[OK] = NULL, [CIE] = NULL, [SRGB] = oklab_to_srgb_fb},
	     [CIE]  = {[OK] = NULL, [CIE] = NULL, [SRGB] = cielab_to_srgb_fb},
	     [SRGB] = {[OK]   = srgb_to_oklab_fb,
		       [CIE]  = srgb_to_cielab_fb,
		       [SRGB] = NULL}}};

static inline void conversion_init() {
	enum InstSet inst = FB;
	if (!features.initialized) {
		query_features(&features);
	}
	if (features.avx) {
		inst = AVX;
	} else if (features.sse) {
		inst = SSE;
	}
	srgb_to_cielab_intern = conversion_table[inst][SRGB][CIE];
	srgb_to_oklab_intern  = conversion_table[inst][SRGB][OK];
	oklab_to_srgb_intern  = conversion_table[inst][OK][SRGB];
	cielab_to_srgb_intern = conversion_table[inst][CIE][SRGB];
}

static inline void check_initialized() {
	static bool conversion_initialized = false;
	if (!conversion_initialized) {
		conversion_init();
		conversion_initialized = true;
	}
}

void srgb_to_oklab(struct Color *restrict colors) {
	check_initialized();
	return srgb_to_oklab_intern(colors);
}
void srgb_to_cielab(struct Color *restrict colors) {
	check_initialized();
	return srgb_to_cielab_intern(colors);
}
void oklab_to_srgb(struct Color *restrict colors) {
	check_initialized();
	return oklab_to_srgb_intern(colors);
}
void cielab_to_srgb(struct Color *restrict colors) {
	check_initialized();
	return cielab_to_srgb_intern(colors);
}
