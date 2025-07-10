#ifndef TEST_SHARED_H
#define TEST_SHARED_H

#include "Color.h"

extern struct Color	colors[];
extern struct sRGB	linears[];
extern struct cieLAB	cielabs[];
extern struct okLAB	oklabs[];
extern struct Grayscale grayscales[];
#define NUM_REF_COL 32

enum ref_indices {
	INDEX,
	SRGBR,
	SRGBG,
	SRGBB,
	LINR,
	LING,
	LINB,
	LABL,
	LABA,
	LABB,
	OKL,
	OKA,
	OKB,
	GRAY
};

#endif
