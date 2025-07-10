#include <stdio.h>

#include "Color.h"
#include "test_shared.h"
#include "unity.h"

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


/* for (int i = 0; i < NUM_REF_COL; ++i) { */
/* 		printf("Calculated: \n"); */
/* 		Color_print(&colors[i]); */
/* 		printf("Scanned: \n"); */
/* 		printf("Linear: (%f, %f, %f)\n", linears[i].r, linears[i].g, */
/* 		       linears[i].b); */
/* 		printf("okLAB: (%f, %f, %f)\n", oklabs[i].l, oklabs[i].a, */
/* 		       oklabs[i].b); */
/* 		printf("cieLAB: (%f, %f, %f)\n", cielabs[i].l, cielabs[i].a, */
/* 		       cielabs[i].b); */
/* 		printf("Grayscale: %f\n", grayscales[i].l); */

/* 		printf("\n\n"); */
/* 	} */
