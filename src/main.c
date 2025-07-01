#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Color.h"
#include "FeatureDetection.h"
#include "Util.h"
#include "Video.h"

/*
 * High level overview:
 * There are four queues, and four threads.
 * The first contains frames passed in from FFMPEG
 * The first thread will pull from this queue, and perform pre-processing if
 * frame_id % fps == 0;
 *
 * Pre-processing consists of down-sampling, generating a saliency map, and pick
 * 5000 or so pixels through a stochastic process using the salience of each
 * pixels as the weight.
 *
 * The first thread then places a stream of 5000 pixels in the second queue.
 *
 * The second thread picks up the stream of pixels from the queue, and runs a
 * kmeans algorithm on it. This way, the kmeans algorithm does not need to take
 * weights into account. K will be picked via sillhouette scores (?)
 * It would be very hard to deal with k being inconsistent accross frames
 *
 * The second thread outputs k centroids.
 *
 * The third thread takes pairs of centroids and interpolates the colors and
 * percentages over n seconds of video.
 *
 * The third thread will then constanly stream interpolation calculated
 * centroids to the fourth queue.
 *
 * The fourth queue will pick up a centroid, create a rectangular buffer of
 * pixels matching the color and weights of each centroid, and write it over top
 * a corner of the original video.
 *
 * These frames are then passed to FFMPEG to write to a new video file.
 */

#define RUNS 131072

static inline int width_from_pixels(int pixels) {
	return (int)(sqrt(pixels) * 4.0f / 3.0f);
}

static inline int height_from_width(int width) {
	return (int)(width * 9.0f / 16.0f);
}

int main(int argc, char **argv) {
	srand(0);

	struct system_features_t features;
	query_features(&features);
	printf("Color is %zu bytes\n", sizeof(struct Color));

	Color col1 = Color_create(rand() * 255, rand() * 255, rand() * 255);
	Color col2 = Color_create(rand() * 255, rand() * 255, rand() * 255);
	Color_calc_spaces(&col1);
	Color_calc_spaces(&col2);
	printf("Euclidean took %lu cycles on average.\n",
	       time_diff(euclidean_diff, col1, col2, RUNS));
	printf("Redmean took %lu cycles on average.\n",
	       time_diff(redmean_diff, col1, col2, RUNS));
	printf("Delta ok took %lu cycles on average.\n",
	       time_diff(delta_ok_diff, col1, col2, RUNS));
	printf("Delta cie76 took %lu cycles on average.\n",
	       time_diff(delta_cie76_diff, col1, col2, RUNS));
	printf("Delta cie94 took %lu cycles on average.\n",
	       time_diff(delta_cie94_diff, col1, col2, RUNS));
	printf("Delta ciede2000 took %lu cycles on average.\n",
	       time_diff(delta_ciede2000_diff, col1, col2, RUNS));
	printf("%" PRIu64 " bytes of l1 cache\n", features.l[0]);
	printf("%" PRIu64 " bytes of l2 cache\n", features.l[1]);
	printf("%" PRIu64 " bytes of l3 cache\n", features.l[2]);

	uint64_t cache	= features.l[0];
	uint64_t pixels = cache / 3;
	int	 width	= width_from_pixels(pixels);
	int	 height = height_from_width(width);
	printf("%" PRIu64 " bytes of l1 cache can hold ~%" PRIu64
	       " pixels. This means that if we want each frame to be held "
	       "entirely in l1 cache at a 16:9 ratio, the frame must be shrunk "
	       "to ~%dx%d\n",
	       cache, pixels, width, height);

	cache += features.l2_shared ? features.l[1]
				    : features.l[1] / features.physical_cores;
	pixels = cache / 3;
	width  = width_from_pixels(pixels);
	height = height_from_width(width);

	printf("%" PRIu64 " bytes of l1 and l2 cache can hold ~%" PRIu64
	       " pixels. This means that if we want each frame to be held "
	       "entirely in l1 and l2 cache at a 16:9 ratio, the frame must be "
	       "shrunk to ~%dx%d\n",
	       cache, pixels, width, height);

	cache += features.l[2] / features.physical_cores;
	pixels = cache / 3;
	width  = width_from_pixels(pixels);
	height = height_from_width(width);

	printf("%" PRIu64 " bytes of l1, l2, and l3 cache can hold ~%" PRIu64
	       " pixels. This means that if we want each frame to be held "
	       "entirely in l1, l2, and l3 cache at a 16:9 ratio, the frame "
	       "must be shrunk to ~%dx%d\n",
	       cache, pixels, width, height);

	/* struct Video in_video; */
	/* open_video_file(&in_video, "vid.webp", "r"); */
	/* while (get_next_frame(&in_video)) { */
	/* 	continue; */
	/* } */
	/* close_video_file(&in_video); */
	return 0;
}
