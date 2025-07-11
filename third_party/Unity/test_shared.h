#ifndef TEST_SHARED_H
#define TEST_SHARED_H

#include "Color.h"

typedef struct E2000_diff {
	float l[2];
	float a[2];
	float b[2];
	float ap[2];
	float cp[2];
	float hp[2];

	float avghp;
	float g;
	float t;
	float sl;
	float sc;
	float sh;
	float rt;
	float diff;
} E2000_diff;

extern struct Color	colors[];
extern struct sRGB	linears[];
extern struct cieLAB	cielabs[];
extern struct okLAB	oklabs[];

extern float	  oklab_diffs[];
extern float	  cie76_diffs[];
extern float	  cie94_diffs[];
extern E2000_diff e2000_diffs_scanned[];
extern E2000_diff e2000_diffs_calc[];

#define NUM_REF_COL 32
#define NUM_DIF ((NUM_REF_COL + 1) * (NUM_REF_COL)) / 2
#define NUM_E2000_PAIR 34
#define CREATION_DELTA 9e-5 /* 0.00009 */
#define DIFF_DELTA 5e-5	    /* 0.00005 */
#define E2000_DELTA 1e-5    /* 0.00001 */

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

enum diff_indices {
	COLOR1_INDEX,
	COLOR2_INDEX,
	OKLAB_DIFF,
	CIE76_DIFF,
	CIE94_DIFF
};

enum e2000_indices {
	PAIR,
	I,
	L,
	A,
	B,
	AP,
	CP,
	HP,
	AVGHP,
	G,
	T,
	SL,
	SC,
	SH,
	RT,
	DIFF
};

#endif
