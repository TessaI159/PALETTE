#include <stdio.h>
#include <stdlib.h>

#include "Color.h"
#include "Util.h"
#include "test_shared.h"
#include "unity.h"

#define RUNS 4194304

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

static inline void print_e2000s() {
	for (int i = 0; i < NUM_E2000_PAIR; ++i) {
		printf("Pair %d\n", i);
		printf("Calculated: \n");
		printf("\t(l,a,b)\t\t\t\t(ap,cp,hp)\t\t\tavghp\t\t(g,t)\t\t\t("
		       "sl,sc,sh)"
		       "\t\t\trt\t\tdiff\n");
		printf("lab1: "
		       "\t(%f,%f,%f)\t(%f,%f,%f)\t%f\t(%f,%f)\t(%f,%f,%f)\t%"
		       "f\t%f\n",
		       e2000_diffs_calc[i].l[0], e2000_diffs_calc[i].a[0],
		       e2000_diffs_calc[i].b[0], e2000_diffs_calc[i].ap[0],
		       e2000_diffs_calc[i].cp[0], e2000_diffs_calc[i].hp[0],
		       e2000_diffs_calc[i].avghp, e2000_diffs_calc[i].g,
		       e2000_diffs_calc[i].t, e2000_diffs_calc[i].sl,
		       e2000_diffs_calc[i].sc, e2000_diffs_calc[i].sh,
		       e2000_diffs_calc[i].rt, e2000_diffs_calc[i].diff);
		printf("lab1: "
		       "\t(%f,%f,%f)\t(%f,%f,%f)\t%f\t(%f,%f)\t(%f,%f,%f)\t%"
		       "f\t%f\n",
		       e2000_diffs_calc[i].l[1], e2000_diffs_calc[i].a[1],
		       e2000_diffs_calc[i].b[1], e2000_diffs_calc[i].ap[1],
		       e2000_diffs_calc[i].cp[1], e2000_diffs_calc[i].hp[1],
		       e2000_diffs_calc[i].avghp, e2000_diffs_calc[i].g,
		       e2000_diffs_calc[i].t, e2000_diffs_calc[i].sl,
		       e2000_diffs_calc[i].sc, e2000_diffs_calc[i].sh,
		       e2000_diffs_calc[i].rt, e2000_diffs_calc[i].diff);

		printf("Scanned: \n");
		printf("\t(l,a,b)\t\t\t\t(ap,cp,hp)\t\t\tavghp\t\t(g,t)\t\t\t("
		       "sl,sc,sh)"
		       "\t\t\trt\t\tdiff\n");
		printf(
		    "lab1: "
		    "\t(%f,%f,%f)\t(%f,%f,%f)\t%f\t(%f,%f)\t(%f,%f,%f)\t%f\t%"
		    "f\n",
		    e2000_diffs_scanned[i].l[0], e2000_diffs_scanned[i].a[0],
		    e2000_diffs_scanned[i].b[0], e2000_diffs_scanned[i].ap[0],
		    e2000_diffs_scanned[i].cp[0], e2000_diffs_scanned[i].hp[0],
		    e2000_diffs_scanned[i].avghp, e2000_diffs_scanned[i].g,
		    e2000_diffs_scanned[i].t, e2000_diffs_scanned[i].sl,
		    e2000_diffs_scanned[i].sc, e2000_diffs_scanned[i].sh,
		    e2000_diffs_scanned[i].rt, e2000_diffs_scanned[i].diff);
		printf(
		    "lab1: "
		    "\t(%f,%f,%f)\t(%f,%f,%f)\t%f\t(%f,%f)\t(%f,%f,%f)\t%f\t%"
		    "f\n",
		    e2000_diffs_scanned[i].l[1], e2000_diffs_scanned[i].a[1],
		    e2000_diffs_scanned[i].b[1], e2000_diffs_scanned[i].ap[1],
		    e2000_diffs_scanned[i].cp[1], e2000_diffs_scanned[i].hp[1],
		    e2000_diffs_scanned[i].avghp, e2000_diffs_scanned[i].g,
		    e2000_diffs_scanned[i].t, e2000_diffs_scanned[i].sl,
		    e2000_diffs_scanned[i].sc, e2000_diffs_scanned[i].sh,
		    e2000_diffs_scanned[i].rt, e2000_diffs_scanned[i].diff);

		printf("\n");
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
	for (int i = 0; i < NUM_E2000_PAIR; ++i) {

		for (int j = 0; j < 2; ++j) {
			const char *tmp_msg = "Pair %d color %d %s off by %f";
			char	    act_msg[256];
			/* lab values */
			sprintf(act_msg, tmp_msg, i, j, "l",
				fabs(e2000_diffs_scanned[i].l[j] -
				     e2000_diffs_calc[i].l[j]));

			TEST_ASSERT_FLOAT_WITHIN_MESSAGE(
			    E2000_DELTA, e2000_diffs_scanned[i].l[j],
			    e2000_diffs_calc[i].l[j], act_msg);

			sprintf(act_msg, tmp_msg, i, j, "a",
				fabs(e2000_diffs_scanned[i].a[j] -
				     e2000_diffs_calc[i].a[j]));

			TEST_ASSERT_FLOAT_WITHIN_MESSAGE(
			    E2000_DELTA, e2000_diffs_scanned[i].a[j],
			    e2000_diffs_calc[i].a[j], act_msg);

			sprintf(act_msg, tmp_msg, i, j, "b",
				fabs(e2000_diffs_scanned[i].b[j] -
				     e2000_diffs_calc[i].b[j]));

			TEST_ASSERT_FLOAT_WITHIN_MESSAGE(
			    E2000_DELTA, e2000_diffs_scanned[i].b[j],
			    e2000_diffs_calc[i].b[j], act_msg);

			/* (ap, cp, hp ) */
			sprintf(act_msg, tmp_msg, i, j, "ap",
				fabs(e2000_diffs_scanned[i].ap[j] -
				     e2000_diffs_calc[i].ap[j]));

			TEST_ASSERT_FLOAT_WITHIN_MESSAGE(
			    E2000_DELTA, e2000_diffs_scanned[i].ap[j],
			    e2000_diffs_calc[i].ap[j], act_msg);

			sprintf(act_msg, tmp_msg, i, j, "cp",
				fabs(e2000_diffs_scanned[i].cp[j] -
				     e2000_diffs_calc[i].cp[j]));

			TEST_ASSERT_FLOAT_WITHIN_MESSAGE(
			    E2000_DELTA, e2000_diffs_scanned[i].cp[j],
			    e2000_diffs_calc[i].cp[j], act_msg);

			sprintf(act_msg, tmp_msg, i, j, "hp",
				fabs(e2000_diffs_scanned[i].hp[j] -
				     e2000_diffs_calc[i].hp[j]));

			TEST_ASSERT_FLOAT_WITHIN_MESSAGE(
			    E2000_DELTA, e2000_diffs_scanned[i].hp[j],
			    e2000_diffs_calc[i].hp[j], act_msg);
		}
		const char *tmp_msg = "Pair %d %s off by %f";
		char	    act_msg[256];

		sprintf(act_msg, tmp_msg, i, "avghp",
			fabs(e2000_diffs_calc[i].avghp -
			     e2000_diffs_scanned[i].avghp));
		TEST_ASSERT_FLOAT_WITHIN_MESSAGE(
		    E2000_DELTA, e2000_diffs_scanned[i].avghp,
		    e2000_diffs_calc[i].avghp, act_msg);

		sprintf(act_msg, tmp_msg, i, "g",
			fabs(e2000_diffs_calc[i].g - e2000_diffs_scanned[i].g));
		TEST_ASSERT_FLOAT_WITHIN_MESSAGE(
		    E2000_DELTA, e2000_diffs_scanned[i].g,
		    e2000_diffs_calc[i].g, act_msg);

		sprintf(act_msg, tmp_msg, i, "t",
			fabs(e2000_diffs_calc[i].t - e2000_diffs_scanned[i].t));
		TEST_ASSERT_FLOAT_WITHIN_MESSAGE(
		    E2000_DELTA, e2000_diffs_scanned[i].t,
		    e2000_diffs_calc[i].t, act_msg);

		sprintf(
		    act_msg, tmp_msg, i, "sl",
		    fabs(e2000_diffs_calc[i].sl - e2000_diffs_scanned[i].sl));
		TEST_ASSERT_FLOAT_WITHIN_MESSAGE(
		    E2000_DELTA, e2000_diffs_scanned[i].sl,
		    e2000_diffs_calc[i].sl, act_msg);

		sprintf(
		    act_msg, tmp_msg, i, "sc",
		    fabs(e2000_diffs_calc[i].sc - e2000_diffs_scanned[i].sc));
		TEST_ASSERT_FLOAT_WITHIN_MESSAGE(
		    E2000_DELTA, e2000_diffs_scanned[i].sc,
		    e2000_diffs_calc[i].sc, act_msg);

		sprintf(
		    act_msg, tmp_msg, i, "sh",
		    fabs(e2000_diffs_calc[i].sh - e2000_diffs_scanned[i].sh));
		TEST_ASSERT_FLOAT_WITHIN_MESSAGE(
		    E2000_DELTA, e2000_diffs_scanned[i].sh,
		    e2000_diffs_calc[i].sh, act_msg);

		sprintf(
		    act_msg, tmp_msg, i, "rt",
		    fabs(e2000_diffs_calc[i].rt - e2000_diffs_scanned[i].rt));
		TEST_ASSERT_FLOAT_WITHIN_MESSAGE(
		    E2000_DELTA, e2000_diffs_scanned[i].rt,
		    e2000_diffs_calc[i].rt, act_msg);

		sprintf(act_msg, tmp_msg, i, "diff",
			fabs(e2000_diffs_calc[i].diff -
			     e2000_diffs_scanned[i].diff));
		TEST_ASSERT_FLOAT_WITHIN_MESSAGE(
		    E2000_DELTA, e2000_diffs_scanned[i].diff,
		    e2000_diffs_calc[i].diff, act_msg);
	}
}

void speed_test(void) {
#ifdef PALETTE_DEBUG
	srand(0);
	Color col1 = Color_create(rand() % 255, rand() % 255, rand() % 255);
	Color_calc_spaces(&col1);
	printf("srgb_to_cielab took %lu cycles on average.\n",
	       time_conv(convert_srgb_to_cielab, col1, RUNS));
	printf("srgb_to_oklab took %lu cycles on average.\n",
	       time_conv(convert_srgb_to_oklab, col1, RUNS));
	printf("srgb_to_grayscale took %lu cycles on average.\n",
	       time_conv(convert_srgb_to_grayscale, col1, RUNS));
	printf("cielab_to_srgb took %lu cycles on average.\n",
	       time_conv(convert_cielab_to_srgb, col1, RUNS));
#endif
}
