#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "unity_config.h"

#include "Color.h"
#include "unity.h"

#define LIN_CRE_DLT 1e-10
#define LAB_CRE_DLT 1e-10
#define OKL_CRE_DLT 1e-4
#define GRY_CRE_DLT 1e-10
#define OKL_DIF_DLT 1e-4
#define C76_DIF_DLT 1e-9
#define C94_DIF_DLT 1e-10
#define C20_DIF_DLT 1e-2

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
Diff_e2000	 ediffs[NUM_PAIRS * 2];
struct Color	 colors[NUM_REF_COL];
struct sRGB	 linears[NUM_REF_COL];
struct cieLAB	 cielabs[NUM_REF_COL];
struct okLAB	 oklabs[NUM_REF_COL];
struct Grayscale grayscales[NUM_REF_COL];

static inline void print_color(const struct Color color) {
	printf("(%f,%f,%f)\n", color.srgb.r, color.srgb.g, color.srgb.b);
}

// Returns number of fields if all are valid, or -1 if any field is empty
static int split_csv(char *line, char **fields, int max_fields) {
	int   count = 0;
	char *start = line;

	while (count < max_fields) {
		char *comma = strchr(start, ',');
		if (comma) {
			*comma		= '\0';
			fields[count++] = start;
			start		= comma + 1;
		} else {
			fields[count++] = start;
			break;
		}
	}

	while (count < max_fields) {
		fields[count++] = NULL;
	}

	// Validation: check for NULL or empty strings
	for (int i = 0; i < max_fields; ++i) {
		if (!fields[i] || fields[i][0] == '\0') {
			return -1;
		}
	}
	return max_fields;
}

static bool parse_data_refs(const char *filename) {
	FILE *fp = fopen(filename, "r");
	if (!fp) {
		fprintf(stderr, "Failed to open %s\n", filename);
		return false;
	}
	char line[MAX_LINE];
	fgets(line, sizeof(line), fp); // header
	while (fgets(line, sizeof(line), fp)) {
		char *fields[NUM_FIELDS_REF] = {0};
		line[strcspn(line, "\r\n")]  = '\0';
		if (split_csv(line, fields, NUM_FIELDS_REF) == -1) {
			fprintf(stderr, "Too few fields in line: %s file: %s\n",
				line, filename);
			continue;
		}
		int idx	       = atoi(fields[INDEX]);
		colors[idx]    = Color_create_norm(strtod(fields[SRGB_R], NULL),
						   strtod(fields[SRGB_G], NULL),
						   strtod(fields[SRGB_B], NULL));
		linears[idx].r = strtod(fields[LINEAR_R], NULL);
		linears[idx].g = strtod(fields[LINEAR_G], NULL);
		linears[idx].b = strtod(fields[LINEAR_B], NULL);
		cielabs[idx].l = strtod(fields[LAB_L], NULL);
		cielabs[idx].a = strtod(fields[LAB_A], NULL);
		cielabs[idx].b = strtod(fields[LAB_B], NULL);
		oklabs[idx].l  = strtod(fields[OKLAB_L], NULL);
		oklabs[idx].a  = strtod(fields[OKLAB_A], NULL);
		oklabs[idx].b  = strtod(fields[OKLAB_B], NULL);
		grayscales[idx].l = strtod(fields[GRAYSCALE], NULL);
	}
	fclose(fp);
	return true;
}

static bool parse_data_diffs(const char *filename) {
	FILE *fp = fopen(filename, "r");
	if (!fp) {
		fprintf(stderr, "Failed to open %s\n", filename);
		return false;
	}
	char line[MAX_LINE];
	fgets(line, sizeof(line), fp); // header
	int ind = 0;
	while (fgets(line, sizeof(line), fp)) {
		char *fields[NUM_FIELDS_DIFF] = {0};
		line[strcspn(line, "\r\n")]   = '\0';
		if (split_csv(line, fields, NUM_FIELDS_DIFF) == -1) {
			fprintf(stderr, "Too few fields in line: %s file: %s\n",
				line, filename);
			continue;
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

static bool parse_data_e2000(const char *filename) {
	FILE *fp = fopen(filename, "r");
	if (!fp) {
		fprintf(stderr, "Failed to open %s\n", filename);
		return false;
	}
	char line[MAX_LINE];
	fgets(line, sizeof(line), fp); // header
	int ind = 0;
	while (fgets(line, sizeof(line), fp)) {
		char *fields[NUM_FIELDS_E2000] = {0};
		line[strcspn(line, "\r\n")]    = '\0';
		if (split_csv(line, fields, NUM_FIELDS_E2000) == -1) {
			fprintf(stderr, "Too few fields in line: %s file: %s\n",
				line, filename);
			continue;
		}
		ediffs[ind].pair = atoi(fields[PAIR]);
		ediffs[ind].i	 = atoi(fields[I]);
		ediffs[ind].l	 = strtod(fields[L], NULL);
		ediffs[ind].a	 = strtod(fields[A], NULL);
		ediffs[ind].b	 = strtod(fields[B], NULL);
		ediffs[ind].ap	 = strtod(fields[AP], NULL);
		ediffs[ind].cp	 = strtod(fields[CP], NULL);
		ediffs[ind].hp	 = strtod(fields[HP], NULL);
		if (ediffs[ind].i == 1) {
			ediffs[ind].avghp = strtod(fields[AVGHP], NULL);
			ediffs[ind].g	  = strtod(fields[G], NULL);
			ediffs[ind].t	  = strtod(fields[T], NULL);
			ediffs[ind].sl	  = strtod(fields[SL], NULL);
			ediffs[ind].sc	  = strtod(fields[SC], NULL);
			ediffs[ind].sh	  = strtod(fields[SH], NULL);
			ediffs[ind].rt	  = strtod(fields[RT], NULL);
			ediffs[ind].e2000 = strtod(fields[E2000], NULL);
		}
		++ind;
	}
	fclose(fp);
	return true;
}
void setUp(void) {
	if (!parse_data_diffs("sharma_pairwise_differences.csv")) {
		TEST_FAIL_MESSAGE("Could not load csv data.");
	}
	if (!parse_data_refs("sharma_reference_colors.csv")) {
		TEST_FAIL_MESSAGE("Could not load csv data.");
	}
	if (!parse_data_e2000("sharma_e2000.csv")) {
		TEST_FAIL_MESSAGE("Could not load csv data.");
	}
}

void tearDown(void) {
}

void test_color_create(void) {
	for (int i = 0; i < NUM_REF_COL; ++i) {
		Color_calc_spaces(&colors[i]);
		TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
		    LIN_CRE_DLT, linears[i].r, colors[i].srgb.r,
		    "Linearized srgb r is too far out of bounds");
		TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
		    LIN_CRE_DLT, linears[i].g, colors[i].srgb.g,
		    "Linearized srgb g is too far out of bounds");
		TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
		    LIN_CRE_DLT, linears[i].b, colors[i].srgb.b,
		    "Linearized srgb b is too far out of bounds");

		TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
		    LAB_CRE_DLT, cielabs[i].l, colors[i].cielab.l,
		    "cielab l is too far out of bounds");
		TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
		    LAB_CRE_DLT, cielabs[i].a, colors[i].cielab.a,
		    "cielab a is too far out of bounds");
		TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
		    LAB_CRE_DLT, cielabs[i].b, colors[i].cielab.b,
		    "cielab b is too far out of bounds");

		TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
		    OKL_CRE_DLT, oklabs[i].l, colors[i].oklab.l,
		    "oklab l is too far out of bounds");
		TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
		    OKL_CRE_DLT, oklabs[i].a, colors[i].oklab.a,
		    "oklab a is too far out of bounds");
		TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
		    OKL_CRE_DLT, oklabs[i].b, colors[i].oklab.b,
		    "oklab b is too far out of bounds");

		TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
		    GRY_CRE_DLT, grayscales[i].l, colors[i].grayscale.l,
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
			    C94_DIF_DLT, cur_diff.cie94,
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
			    OKL_DIF_DLT, cur_diff.oklab,
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
			    C76_DIF_DLT, cur_diff.cie76,
			    delta_cie76_diff(&colors[cur_diff.color1_index],
					     &colors[cur_diff.color2_index]),
			    "cie76 diff is too far off.");
		}
	}
}

void test_ciede2000_diff(void) {
	for (int i = 0; i < NUM_PAIRS * 2; i += 2) {
		struct cieLAB_test cols[2] = {0};
		for (int j = 0; j < 2; ++j) {
			cols[j].l = ediffs[i + j].l;
			cols[j].a = ediffs[i + j].a;
			cols[j].b = ediffs[i + j].b;
		}
		delta_ciede2000_diff_fast(&cols[0], &cols[1]);
		for (int j = 0; j < 2; ++j) {
			TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
			    C20_DIF_DLT, ediffs[i + j].ap, cols[j].ap,
			    "message");
			TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
			    C20_DIF_DLT, ediffs[i + j].cp, cols[j].cp,
			    "message");
			TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
			    C20_DIF_DLT, ediffs[i + j].hp, cols[j].hp,
			    "message");
			if (j == 0) {
				TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
				    C20_DIF_DLT, ediffs[i + j].avghp,
				    cols[j].avghp, "message");
				TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
				    C20_DIF_DLT, ediffs[i + j].g, cols[j].g,
				    "message");
				TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
				    C20_DIF_DLT, ediffs[i + j].t, cols[j].t,
				    "message");
				TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
				    C20_DIF_DLT, ediffs[i + j].sl, cols[j].sl,
				    "message");
				TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
				    C20_DIF_DLT, ediffs[i + j].sc, cols[j].sc,
				    "message");
				TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
				    C20_DIF_DLT, ediffs[i + j].sh, cols[j].sh,
				    "message");
				TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
				    C20_DIF_DLT, ediffs[i + j].rt, cols[j].rt,
				    "message");
				TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(
				    C20_DIF_DLT, ediffs[i + j].e2000,
				    cols[j].e2000, "message");
			}
		}
	}
}
