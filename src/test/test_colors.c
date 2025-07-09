#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "unity_config.h"

#include "Color.h"
#include "unity.h"

#define DELTA 9e-5

#define NUM_REF_COL 22
#define NUM_FIELDS_DIFF 6
#define NUM_FIELDS_REF 14

#define MAX_LINE 2048
#define NUM_DIFFS ((NUM_REF_COL - 1) * (NUM_REF_COL) / 2)

#define NUM_FIELDS_E2000 16
#define NUM_PAIRS 34

typedef struct Diff {
	int    color1_index;
	int    color2_index;
	double oklab;
	double cie76;
	double cie94;
	double ciede2000;
} Diff;

typedef struct Diff_e2000 {
	int    pair;
	int    i;
	double l;
	double a;
	double b;
	double ap;
	double cp;
	double hp;
	double avghp;
	double g;
	double t;
	double sl;
	double sc;
	double sh;
	double rt;
	double e2000;
} Diff_e2000;

enum diff_csv_indices {
	COLOR1_INDEX = 0,
	COLOR2_INDEX,
	OKLAB_DIFF,
	CIE76_DIFF,
	CIE94_DIFF
};

enum ref_csv_indices {
	INDEX = 0,
	SRGB_R,
	SRGB_G,
	SRGB_B,
	LINEAR_R,
	LINEAR_G,
	LINEAR_B,
	LAB_L,
	LAB_A,
	LAB_B,
	OKLAB_L,
	OKLAB_A,
	OKLAB_B,
	GRAYSCALE
};

enum e2000_csv_indices {
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
	E2000
};

Diff		 diffs[NUM_DIFFS];
Diff_e2000	 ediffs[NUM_PAIRS];
struct Color	 colors[NUM_REF_COL];
struct sRGB	 linears[NUM_REF_COL];
struct cieLAB	 cielabs[NUM_REF_COL];
struct okLAB	 oklabs[NUM_REF_COL];
struct Grayscale grayscales[NUM_REF_COL];

static inline void print_color(const struct Color color) {
	printf("(%f,%f,%f)\n", color.srgb.r, color.srgb.g, color.srgb.b);
}

static bool parse_data_refs(struct Color *colors, struct sRGB *linears,
			    struct cieLAB *cielabs, struct okLAB *oklabs,
			    struct Grayscale *grayscales,
			    const char	     *filename) {
	FILE *fp = fopen(filename, "r");
	if (!fp) {
		fprintf(stderr, "Failed to open %s\n", filename);
		return false;
	}

	char line[MAX_LINE];
	fgets(line, sizeof(line), fp);
	while (fgets(line, sizeof(line), fp)) {
		char *fields[NUM_FIELDS_REF] = {0};
		char *token		     = strtok(line, ",\n");
		int   i			     = 0;
		while (token && i < NUM_FIELDS_REF) {
			fields[i++] = token;
			token	    = strtok(NULL, ",\n");
		}
		colors[atoi(fields[INDEX])] = Color_create_norm(
		    strtod(fields[SRGB_R], NULL), strtod(fields[SRGB_G], NULL),
		    strtod(fields[SRGB_B], NULL));
		linears[atoi(fields[INDEX])].r = strtod(fields[LINEAR_R], NULL);
		linears[atoi(fields[INDEX])].g = strtod(fields[LINEAR_G], NULL);
		linears[atoi(fields[INDEX])].b = strtod(fields[LINEAR_B], NULL);
		cielabs[atoi(fields[INDEX])].l = strtod(fields[LAB_L], NULL);
		cielabs[atoi(fields[INDEX])].a = strtod(fields[LAB_A], NULL);
		cielabs[atoi(fields[INDEX])].b = strtod(fields[LAB_B], NULL);
		oklabs[atoi(fields[INDEX])].l  = strtod(fields[OKLAB_L], NULL);
		oklabs[atoi(fields[INDEX])].a  = strtod(fields[OKLAB_A], NULL);
		oklabs[atoi(fields[INDEX])].b  = strtod(fields[OKLAB_B], NULL);
		grayscales[atoi(fields[INDEX])].l =
		    strtod(fields[GRAYSCALE], NULL);
	}
	fclose(fp);
	return true;
}

static bool parse_data_diffs(Diff *diffs, const char *filename) {
	FILE *fp = fopen(filename, "r");
	if (!fp) {
		fprintf(stderr, "Failed to open %s\n", filename);
		return false;
	}
	char line[MAX_LINE];
	fgets(line, sizeof(line), fp);
	int ind = 0;
	while (fgets(line, sizeof(line), fp)) {
		char *fields[NUM_FIELDS_DIFF] = {0};
		char *token		      = strtok(line, ",\n");
		int   i			      = 0;
		while (token && i < NUM_FIELDS_DIFF) {
			fields[i++] = token;
			token	    = strtok(NULL, ",\n");
		}
		diffs[ind].color1_index = atoi(fields[COLOR1_INDEX]);
		diffs[ind].color2_index = atoi(fields[COLOR2_INDEX]);
		diffs[ind].oklab	= strtod(fields[OKLAB_DIFF], NULL);
		diffs[ind].cie76	= strtod(fields[CIE76_DIFF], NULL);
		diffs[ind].cie94	= strtod(fields[CIE94_DIFF], NULL);
		++ind;
	}
	fclose(fp);
	return true;
}

static bool parse_data_e2000(Diff_e2000 *diffs, const char *filename) {
	FILE *fp = fopen(filename, "r");
	if (!fp) {
		fprintf(stderr, "Failed to open %s\n", filename);
		return false;
	}
	char line[MAX_LINE];
	fgets(line, sizeof(line), fp);
	int ind = 0;
	while (fgets(line, sizeof(line), fp)) {
		char *fields[NUM_FIELDS_E2000];
		char *token = strtok(line, ",\n");
		int   i	    = 0;
		while (token && i < NUM_FIELDS_E2000) {
			fields[i++] = token;
			token	    = strtok(NULL, ",\n");
		}
		ediffs[ind].pair = atoi(fields[PAIR]);
		ediffs[ind].i	 = atoi(fields[I]);
		ediffs[ind].l	 = strtod(fields[L], NULL);
		ediffs[ind].a	 = strtod(fields[A], NULL);
		ediffs[ind].b	 = strtod(fields[B], NULL);
		ediffs[ind].ap	 = strtod(fields[AP], NULL);
		ediffs[ind].cp	 = strtod(fields[CP], NULL);
		ediffs[ind].hp	 = strtod(fields[HP], NULL);
		if (ediffs[ind].i >> 1 == 1) {
			ediffs[ind].g	  = strtod(fields[G], NULL);
			ediffs[ind].t	  = strtod(fields[T], NULL);
			ediffs[ind].sl	  = strtod(fields[SL], NULL);
			ediffs[ind].sc	  = strtod(fields[SC], NULL);
			ediffs[ind].sh	  = strtod(fields[SH], NULL);
			ediffs[ind].rt	  = strtod(fields[RT], NULL);
			ediffs[ind].e2000 = strtod(fields[E2000], NULL);
		}
	}
	fclose(fp);
	return true;
}

void setUp(void) {
	if (!parse_data_diffs(diffs, "sharma_pairwise_differences.csv")) {
		TEST_FAIL_MESSAGE("Could not load csv data.");
	}
	if (!parse_data_refs(colors, linears, cielabs, oklabs, grayscales,
			     "sharma_reference_colors.csv")) {
		TEST_FAIL_MESSAGE("Could not load csv data.");
	}
	if (!parse_data_e2000(ediffs, "sharma_e2000.csv")) {
		TEST_FAIL_MESSAGE("Could not load csv data.");
	}
}

void tearDown(void) {
}

void test_color_create(void) {
	for (int i = 0; i < NUM_REF_COL; ++i) {
		Color_calc_spaces(&colors[i]);
		TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
		    DELTA, linears[i].r, colors[i].srgb.r,
		    "Linearized srgb r is too far out of bounds");
		TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
		    DELTA, linears[i].g, colors[i].srgb.g,
		    "Linearized srgb g is too far out of bounds");
		TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
		    DELTA, linears[i].b, colors[i].srgb.b,
		    "Linearized srgb b is too far out of bounds");

		TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
		    DELTA, cielabs[i].l, colors[i].cielab.l,
		    "cielab l is too far out of bounds");
		TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
		    DELTA, cielabs[i].a, colors[i].cielab.a,
		    "cielab a is too far out of bounds");
		TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
		    DELTA, cielabs[i].b, colors[i].cielab.b,
		    "cielab b is too far out of bounds");

		TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
		    DELTA, oklabs[i].l, colors[i].oklab.l,
		    "oklab l is too far out of bounds");
		TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
		    DELTA, oklabs[i].a, colors[i].oklab.a,
		    "oklab a is too far out of bounds");
		TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
		    DELTA, oklabs[i].b, colors[i].oklab.b,
		    "oklab b is too far out of bounds");

		TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
		    DELTA, grayscales[i].l, colors[i].grayscale.l,
		    "grayscale is too far out of bounds");
	}
}

void test_color_check_flags(void) {
	struct Color ref = Color_create(1, 1, 1);
	for (int i = 0; i < NUM_REF_COL; ++i) {
		TEST_ASSERT_EQUAL_UINT8_MESSAGE(1, colors[i].valid_spaces,
						"Wrong marked colors");
		delta_ok_diff(&colors[i], &ref);
	}
	
	for (int i = 0; i < NUM_REF_COL; ++i) {
	  TEST_ASSERT_EQUAL_UINT8_MESSAGE(3, colors[i].valid_spaces,
						"Wrong marked colors");
	  delta_cie76_diff(&colors[i], &ref);
	}
	
	for (int i = 0; i < NUM_REF_COL; ++i) {
	  TEST_ASSERT_EQUAL_UINT8_MESSAGE(7, colors[i].valid_spaces,
						"Wrong marked colors");
	  convert_srgb_to_grayscale(&colors[i]);
	}

	for (int i = 0; i < NUM_REF_COL; ++i) {
	    TEST_ASSERT_EQUAL_UINT8_MESSAGE(15, colors[i].valid_spaces,
						"Wrong marked colors");
	}
}

void test_cie94_diff(void) {
	int	  diff_num  = 0;
	const int ref_const = NUM_REF_COL + NUM_REF_COL - 1;
	for (int c1 = 0; c1 < NUM_REF_COL; ++c1) {
		for (int c2 = c1 + 1; c2 < NUM_REF_COL; ++c2) {
			diff_num = (c1 * (ref_const - c1)) / 2 + (c2 - c1 - 1);
			Diff cur_diff = diffs[diff_num];
			TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
			    DELTA, cur_diff.cie94,
			    delta_cie94_diff(&colors[cur_diff.color1_index],
					     &colors[cur_diff.color2_index]),
			    "cie94 diff is too far off.");
		}
	}
}

void test_oklab_diff(void) {
	int	  diff_num  = 0;
	const int ref_const = NUM_REF_COL + NUM_REF_COL - 1;
	for (int c1 = 0; c1 < NUM_REF_COL; ++c1) {
		for (int c2 = c1 + 1; c2 < NUM_REF_COL; ++c2) {
			diff_num = (c1 * (ref_const - c1)) / 2 + (c2 - c1 - 1);

			Diff cur_diff = diffs[diff_num];
			TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
			    DELTA, cur_diff.oklab,
			    delta_ok_diff(&colors[cur_diff.color1_index],
					  &colors[cur_diff.color2_index]),
			    "oklab diff is too far off.");
		}
	}
}

void test_cie76_diff(void) {
	int	  diff_num  = 0;
	const int ref_const = NUM_REF_COL + NUM_REF_COL - 1;
	for (int c1 = 0; c1 < NUM_REF_COL; ++c1) {
		for (int c2 = c1 + 1; c2 < NUM_REF_COL; ++c2) {
			diff_num = (c1 * (ref_const - c1)) / 2 + (c2 - c1 - 1);
			Diff cur_diff = diffs[diff_num];
			TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
			    DELTA, cur_diff.cie76,
			    delta_cie76_diff(&colors[cur_diff.color1_index],
					     &colors[cur_diff.color2_index]),
			    "cie76 diff is too far off.");
		}
	}
}

void test_ciede2000_diff(void) {
}
