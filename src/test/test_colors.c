#include <stdio.h>

#include "Color.h"
#include "test_shared.h"
#include "unity.h"

static inline void print_refs() {
	for (size_t i = 0; i < NUM_REF_COL; ++i) {
		printf("Calculated: \n");
		Color_print(&colors[i]);
		printf("Scanned: \n");
		printf("Linear: (%f, %f, %f)\n", linears[i].r, linears[i].g,
		       linears[i].b);
		printf("okLAB: (%f, %f, %f)\n", oklabs[i].l, oklabs[i].a,
		       oklabs[i].b);
		printf("cieLAB: (%f, %f, %f)\n", cielabs[i].l, cielabs[i].a,
		       cielabs[i].b);
		printf("Grayscale: %f\n", grayscales[i].l);

		printf("\n\n");
	}
}

static inline void print_diffs() {
	size_t i = 0;
	for (size_t i1 = 0; i1 < NUM_REF_COL; ++i1) {
		for (size_t i2 = i1 + 1; i2 < NUM_REF_COL; ++i2) {
			printf("diff (%zu, %zu): \nCalculated: \noklab: "
			       "%f\ncie76: %f\ncie94 "
			       "%f\n",
			       i1, i2, delta_ok_diff(&colors[i1], &colors[i2]),
			       delta_cie76_diff(&colors[i1], &colors[i2]),
			       delta_cie94_diff(&colors[i1], &colors[i2]));

			printf("Scanned: \n");
			printf("oklab: %f\ncie76 %f\ncie94 "
			       "%f\n\n",
			       oklab_diffs[i], cie76_diffs[i], cie94_diffs[i]);
			++i;
		}
	}
}

/* TODO (Tess): As with the others, print the calculated and the scanned values
 */
static inline void print_e2000s() {

	for (size_t i = 0; i < NUM_E2000_PAIR; ++i) {
		E2000_diff calc_diff;
		calc_diff.l[0] = colors[i].cielab.l;
		calc_diff.a[0] = colors[i].cielab.a;
		calc_diff.b[0] = colors[i].cielab.b;

		calc_diff.l[1] = colors[i + 1].cielab.l;
		calc_diff.a[1] = colors[i + 1].cielab.a;
		calc_diff.b[1] = colors[i + 1].cielab.b;
		delta_ciede2000_diff_fast(&calc_diff);
		printf("Computed: \n");
		printf("Lab1: (%f, %f, %f)\nLab2: (%f, %f, %f)\n",
		       calc_diff.l[0], calc_diff.a[0], calc_diff.b[0],
		       calc_diff.l[1], calc_diff.a[1], calc_diff.b[1]);
	}
}

void test_color_check_flags(void) {
	struct Color test1 = Color_create(124, 191, 171);
	struct Color test2 = Color_create(255, 193, 204);

	TEST_ASSERT_TRUE(Color_has_space(&test1, COLOR_SRGB));
	TEST_ASSERT_TRUE(Color_has_space(&test2, COLOR_SRGB));

	delta_ok_diff(&test1, &test2);

	TEST_ASSERT_TRUE(Color_has_space(&test1, COLOR_OKLAB));
	TEST_ASSERT_TRUE(Color_has_space(&test2, COLOR_OKLAB));

	delta_ciede2000_diff(&test1, &test2);

	TEST_ASSERT_TRUE(Color_has_space(&test1, COLOR_CIELAB));
	TEST_ASSERT_TRUE(Color_has_space(&test2, COLOR_CIELAB));

	convert_srgb_to_grayscale(&test1);
	convert_srgb_to_grayscale(&test2);

	TEST_ASSERT_TRUE(Color_has_space(&test1, COLOR_GRAY));
	TEST_ASSERT_TRUE(Color_has_space(&test2, COLOR_GRAY));
}

void test_color_create(void) {
	for (size_t i = 0; i < NUM_REF_COL; ++i) {
		Color_calc_spaces(&colors[i]);
		const char *tmp_msg = "Color %d %s %c off by %f";
		char	    act_msg[256];
		/* Linear sRGBs */
		sprintf(act_msg, tmp_msg, i, "linear srgb", 'r',
			fabs(linears[i].r - colors[i].srgb.r));
		TEST_ASSERT_FLOAT_WITHIN_MESSAGE(CREATION_DELTA, linears[i].r,
						 colors[i].srgb.r, act_msg);
		sprintf(act_msg, tmp_msg, i, "linear srgb", 'g',
			fabs(linears[i].g - colors[i].srgb.g));
		TEST_ASSERT_FLOAT_WITHIN_MESSAGE(CREATION_DELTA, linears[i].g,
						 colors[i].srgb.g, act_msg);
		sprintf(act_msg, tmp_msg, i, "linear srgb", 'b',
			fabs(linears[i].b - colors[i].srgb.b));
		TEST_ASSERT_FLOAT_WITHIN_MESSAGE(CREATION_DELTA, linears[i].b,
						 colors[i].srgb.b, act_msg);

		/* cieLABs */
		sprintf(act_msg, tmp_msg, i, "cielab", 'l',
			fabs(cielabs[i].l - colors[i].cielab.l));
		TEST_ASSERT_FLOAT_WITHIN_MESSAGE(CREATION_DELTA, cielabs[i].l,
						 colors[i].cielab.l, act_msg);
		sprintf(act_msg, tmp_msg, i, "cielab", 'a',
			fabs(cielabs[i].a - colors[i].cielab.a));
		TEST_ASSERT_FLOAT_WITHIN_MESSAGE(CREATION_DELTA, cielabs[i].a,
						 colors[i].cielab.a, act_msg);
		sprintf(act_msg, tmp_msg, i, "cielab", 'b',
			fabs(cielabs[i].b - colors[i].cielab.b));
		TEST_ASSERT_FLOAT_WITHIN_MESSAGE(CREATION_DELTA, cielabs[i].b,
						 colors[i].cielab.b, act_msg);

		/* okLABs */
		sprintf(act_msg, tmp_msg, i, "oklab", 'l',
			fabs(oklabs[i].l - colors[i].oklab.l));
		TEST_ASSERT_FLOAT_WITHIN_MESSAGE(CREATION_DELTA, oklabs[i].l,
						 colors[i].oklab.l, act_msg);
		sprintf(act_msg, tmp_msg, i, "oklab", 'a',
			fabs(oklabs[i].a - colors[i].oklab.a));
		TEST_ASSERT_FLOAT_WITHIN_MESSAGE(CREATION_DELTA, oklabs[i].a,
						 colors[i].oklab.a, act_msg);
		sprintf(act_msg, tmp_msg, i, "oklab", 'b',
			fabs(oklabs[i].b - colors[i].oklab.b));
		TEST_ASSERT_FLOAT_WITHIN_MESSAGE(CREATION_DELTA, oklabs[i].b,
						 colors[i].oklab.b, act_msg);

		/* Grayscale */
		sprintf(act_msg, tmp_msg, i, "grayscale", 'l',
			fabs(grayscales[i].l - colors[i].grayscale.l));
		TEST_ASSERT_FLOAT_WITHIN_MESSAGE(
		    CREATION_DELTA, grayscales[i].l, colors[i].grayscale.l,
		    act_msg);
	}
}

void test_cie76_diff(void) {
	size_t i = 0;
	for (size_t i1 = 0; i1 < NUM_REF_COL; ++i1) {
		for (size_t i2 = i1 + 1; i2 < NUM_REF_COL; ++i2) {
			float cur_cie76 =
			    delta_cie76_diff(&colors[i1], &colors[i2]);

			const char *tmp_msg =
			    "delta_cie76_diff(%zu, %zu) off by  %f.";
			char act_msg[128];
			sprintf(act_msg, tmp_msg, i1, i2,
				fabs(cur_cie76 - cie76_diffs[i]));
			TEST_ASSERT_FLOAT_WITHIN_MESSAGE(
			    DIFF_DELTA, cie76_diffs[i], cur_cie76, act_msg);
			++i;
		}
	}
}

void test_cie94_diff(void) {
	size_t i = 0;
	for (size_t i1 = 0; i1 < NUM_REF_COL; ++i1) {
		for (size_t i2 = i1 + 1; i2 < NUM_REF_COL; ++i2) {
			float cur_cie94 =
			    delta_cie94_diff(&colors[i1], &colors[i2]);

			const char *tmp_msg =
			    "delta_cie94_diff(%zu, %zu) off by  %f.";
			char act_msg[128];
			sprintf(act_msg, tmp_msg, i1, i2,
				fabs(cur_cie94 - cie94_diffs[i]));
			TEST_ASSERT_FLOAT_WITHIN_MESSAGE(
			    DIFF_DELTA, cie94_diffs[i], cur_cie94, act_msg);
			++i;
		}
	}
}

void test_ok_diff(void) {
	size_t i = 0;
	for (size_t i1 = 0; i1 < NUM_REF_COL; ++i1) {
		for (size_t i2 = i1 + 1; i2 < NUM_REF_COL; ++i2) {
			float cur_ok = delta_ok_diff(&colors[i1], &colors[i2]);

			const char *tmp_msg =
			    "delta_ok_diff(%zu, %zu) off by  %f.";
			char act_msg[128];
			sprintf(act_msg, tmp_msg, i1, i2,
				fabs(cur_ok - oklab_diffs[i]));
			TEST_ASSERT_FLOAT_WITHIN_MESSAGE(
			    DIFF_DELTA, oklab_diffs[i], cur_ok, act_msg);
			++i;
		}
	}
}

void test_ciede2000_diff(void) {
	print_e2000s();
	for (size_t i = 0; i < NUM_E2000_PAIR; ++i) {
	}
}
