#include <stdbool.h>
#include <stdlib.h>

#include "unity_config.h"

#include "Color.h"
#include "csv_parser.h"
#include "test_shared.h"
#include "unity.h"

#define MAX_LINE 512
#define MAX_FIELDS 64

struct Color	 colors[NUM_REF_COL];
struct sRGB	 linears[NUM_REF_COL];
struct cieLAB	 cielabs[NUM_REF_COL];
struct okLAB	 oklabs[NUM_REF_COL];
struct Grayscale grayscales[NUM_REF_COL];

static inline bool parse_ref(const char *filename, struct Color *colors,
			     struct sRGB *linears, struct cieLAB *cielabs,
			     struct okLAB     *oklabs,
			     struct Grayscale *grayscales) {
	FILE *fp = fopen(filename, "r");
	if (!fp) {
		fprintf(stderr, "Unable to open %s\n", filename);
		return false;
	}
	char  line[MAX_LINE];
	char *fields[MAX_FIELDS];
	fgets(line, sizeof(line), fp);
	while (fgets(line, sizeof(line), fp)) {
		int n = csv_parse_line(line, fields, MAX_FIELDS);
		if (n < 0) {
			fprintf(stderr, "Unable to parse %s\n", filename);
			return false;
		}
		int index = atoi(fields[INDEX]);

		colors[index] = Color_create_norm(strtof(fields[SRGBR], NULL),
						  strtof(fields[SRGBG], NULL),
						  strtof(fields[SRGBB], NULL));
		linears[index].r = strtof(fields[LINR], NULL);
		linears[index].g = strtof(fields[LING], NULL);
		linears[index].b = strtof(fields[LINB], NULL);

		cielabs[index].l = strtof(fields[LABL], NULL);
		cielabs[index].a = strtof(fields[LABA], NULL);
		cielabs[index].b = strtof(fields[LABB], NULL);

		oklabs[index].l = strtof(fields[OKL], NULL);
		oklabs[index].a = strtof(fields[OKA], NULL);
		oklabs[index].b = strtof(fields[OKB], NULL);

		grayscales[index].l = strtof(fields[GRAY], NULL);
	}
	return true;
}

void setUp(void) {
	parse_ref("sharma_reference_colors.csv", colors, linears, cielabs,
		  oklabs, grayscales);
}

void tearDown(void) {
}

extern void test_queue_setup(void);
extern void test_queue_push(void);
extern void test_queue_pop(void);
extern void test_queue_counters(void);
extern void test_color_create(void);
extern void test_cie94_diff(void);
extern void test_cie76_diff(void);
extern void test_oklab_diff(void);
extern void test_ciede2000_diff(void);
extern void test_color_check_flags(void);

int main(void) {
	UNITY_BEGIN();
	RUN_TEST(test_queue_setup);
	RUN_TEST(test_queue_push);
	RUN_TEST(test_queue_pop);
	RUN_TEST(test_queue_counters);
	RUN_TEST(test_color_check_flags);
	/* RUN_TEST(test_color_create); */
	/* RUN_TEST(test_cie76_diff); */
	/* RUN_TEST(test_cie94_diff); */
	/* RUN_TEST(test_ciede2000_diff); */
	/* RUN_TEST(test_oklab_diff); */
	return UNITY_END();
}
