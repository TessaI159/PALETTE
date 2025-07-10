#include <stdbool.h>
#include <stdlib.h>

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

float	   oklab_diffs[NUM_DIF];
float	   cie76_diffs[NUM_DIF];
float	   cie94_diffs[NUM_DIF];
E2000_diff e2000_diffs_scanned[NUM_E2000_PAIR];
E2000_diff e2000_diffs_calc[NUM_E2000_PAIR];

static inline bool parse_ref(const char *filename) {
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

static inline bool parse_diff(const char *filename) {
	FILE *fp = fopen(filename, "r");
	if (!fp) {
		fprintf(stderr, "Unable to open %s\n", filename);
		return false;
	}
	char  line[MAX_LINE];
	char *fields[MAX_FIELDS];
	fgets(line, sizeof(line), fp);
	int i = 0;
	while (fgets(line, sizeof(line), fp)) {
		int n = csv_parse_line(line, fields, MAX_FIELDS);
		if (n < 0) {
			fprintf(stderr, "Unable to parse %s\n", filename);
			return false;
		}
		oklab_diffs[i] = strtof(fields[OKLAB_DIFF], NULL);
		cie76_diffs[i] = strtof(fields[CIE76_DIFF], NULL);
		cie94_diffs[i] = strtof(fields[CIE94_DIFF], NULL);
		++i;
	}
	return true;
}

static inline bool parse_e2000(const char *filename) {
	FILE *fp = fopen(filename, "r");
	if (!fp) {
		fprintf(stderr, "Unable to open %s\n", filename);
		return false;
	}
	char  line[MAX_LINE];
	char *fields[MAX_FIELDS];
	int   index;
	fgets(line, sizeof(line), fp);
	while (fgets(line, sizeof(line), fp)) {
		int n = csv_parse_line(line, fields, MAX_FIELDS);
		if (n < 0) {
			fprintf(stderr, "Unable to parse %s\n", filename);
			return false;
		}
		index = atoi(fields[PAIR]) - 1;

		e2000_diffs_scanned[index].l[0]	 = strtof(fields[L], NULL);
		e2000_diffs_scanned[index].a[0]	 = strtof(fields[A], NULL);
		e2000_diffs_scanned[index].b[0]	 = strtof(fields[B], NULL);
		e2000_diffs_scanned[index].ap[0] = strtof(fields[AP], NULL);
		e2000_diffs_scanned[index].cp[0] = strtof(fields[CP], NULL);
		e2000_diffs_scanned[index].hp[0] = strtof(fields[HP], NULL);

		e2000_diffs_scanned[index].avghp = strtof(fields[AVGHP], NULL);
		e2000_diffs_scanned[index].g	 = strtof(fields[G], NULL);
		e2000_diffs_scanned[index].t	 = strtof(fields[T], NULL);
		e2000_diffs_scanned[index].sl	 = strtof(fields[SL], NULL);
		e2000_diffs_scanned[index].sc	 = strtof(fields[SC], NULL);
		e2000_diffs_scanned[index].sh	 = strtof(fields[SH], NULL);
		e2000_diffs_scanned[index].rt	 = strtof(fields[RT], NULL);
		e2000_diffs_scanned[index].diff	 = strtof(fields[DIFF], NULL);

		fgets(line, sizeof(line), fp);
		n = csv_parse_line(line, fields, MAX_FIELDS);
		if (n < 0) {
			fprintf(stderr, "Unable to parse %s\n", filename);
			return false;
		}
		e2000_diffs_scanned[index].l[1]	 = strtof(fields[L], NULL);
		e2000_diffs_scanned[index].a[1]	 = strtof(fields[A], NULL);
		e2000_diffs_scanned[index].b[1]	 = strtof(fields[B], NULL);
		e2000_diffs_scanned[index].ap[1] = strtof(fields[AP], NULL);
		e2000_diffs_scanned[index].cp[1] = strtof(fields[CP], NULL);
		e2000_diffs_scanned[index].hp[1] = strtof(fields[HP], NULL);
	}
	for (int i = 0; i < index + 1; ++i) {
		for (int j = 0; j < 2; ++j) {
			e2000_diffs_calc[i].l[j] = e2000_diffs_scanned[i].l[j];
			e2000_diffs_calc[i].a[j] = e2000_diffs_scanned[i].a[j];
			e2000_diffs_calc[i].b[j] = e2000_diffs_scanned[i].b[j];
		}
		delta_ciede2000_diff_fast(&e2000_diffs_calc[i]);
	}
	return true;
}

void setUp(void) {
	if (!parse_ref("sharma_reference_colors.csv")) {
		TEST_FAIL_MESSAGE("Could not parse references");
	}
	if (!parse_diff("sharma_pairwise_differences.csv")) {
		TEST_FAIL_MESSAGE("Could not parse differences");
	}
	if (!parse_e2000("sharma_e2000.csv")) {
		TEST_FAIL_MESSAGE("Could not parse E2000");
	}
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
extern void test_ok_diff(void);
extern void test_ciede2000_diff(void);
extern void test_color_check_flags(void);

int main(void) {
	UNITY_BEGIN();
	RUN_TEST(test_queue_setup);
	RUN_TEST(test_queue_push);
	RUN_TEST(test_queue_pop);
	RUN_TEST(test_queue_counters);
	RUN_TEST(test_color_check_flags);
	RUN_TEST(test_color_create);
	RUN_TEST(test_cie76_diff);
	RUN_TEST(test_cie94_diff);
	RUN_TEST(test_ok_diff);
	RUN_TEST(test_ciede2000_diff);
	return UNITY_END();
}
