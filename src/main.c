#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Color.h"
#include "FeatureDetection.h"
#include "Util.h"
#include "Video.h"

/* TODO (Tess): Consider throwing GPU compute capabilities in here somewhere */

/*
 * High level overview:
 * There are three queues, and four threads.
 * The first queue contains frames passed in from FFMPEG
 * FFMPEG will be picky about what frames it queues if:
 * The frame is a key frame
 * It has been >=1 second (in video time) since the last frame
 * The first thread will pull from this queue, and perform pre-processing
 *
 * Pre-processing consists of down-sampling, generating a saliency map, and
 * picking 5000 or so pixels through a stochastic process using the salience of
 * each pixel as the weight.
 *
 * The first thread then queues those 5000 pixels
 *
 * The second and third threads pick up the stream of pixels from the queue, and
 * run kmeans algorithm on them. This way, the kmeans algorithm does not need
 * to take weights into account.
 *
 * TODO (Tess): K might be picked via sillhouette scores (?) It would be very
 * hard to deal with k being inconsistent accross frames. Although, I could just
 * interpolate sizes from 0 -> size or size -> 0, but how would the color be
 * interpolated?
 *
 * The second and third threads queue centroids
 * TODO (Tess): How will order be taken care of here? Will threads 2 and 3 just
 * throw things into the queue willy-nilly and force the third thread to search
 * the queue until it finds the right partner frame? Or will we hang thread 2/3
 * until the latest centroid has the latest id?
 *
 * The third thread takes pairs of centroids and interpolates the colors and
 * percentages over n seconds of video, where n is frame_n.pts - frame_n-1.pts
 *
 * The third thread will then constanly stream interpolation calculated
 * centroids to disk. For instance, if we have CFR of 24, then for each pair of
 * centroids, the third thread will create 22 and output 24 centroids.
 *
 * Close the video stream.
 *
 * After this process is complete, we will open the video again, but this time
 * for a very quick pass. We will decode frame -> add dominant colors to video
 * -> encode frame -> write to disk -> Delete intermediary file.
 *
 *
 * Goal: process 60 frames per second, minimum. ie if the video is 3600 frames
 * long, the whole process should take about 1 minute, INCLUDING WRITING THE NEW
 * VIDEO TO DISK.
 */

#define RUNS 4194304

static inline int width_from_pixels(int pixels) {
	return round(sqrt(pixels) * 4.0f / 3.0f);
}

static inline int height_from_width(int width) {
	return round(9.0f / 16.0f * width);
}

int main(int argc, char **argv) {
	/* srand(0); */

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

	uint64_t cache	= features.l[1];
	uint64_t pixels = cache / 3;
	int	 width	= width_from_pixels(pixels);
	int	 height = height_from_width(width);
	printf("%" PRIu64 " bytes of l1 cache can hold ~%" PRIu64
	       " pixels. This means that if we want each frame to be held "
	       "entirely in l1 cache at a 16:9 ratio, the frame must be shrunk "
	       "to ~%dx%d\n",
	       cache, pixels, width, height);

	cache += features.l2_shared ? features.l[2]
				    : features.l[2] / features.physical_cores;
	pixels = cache / 3;
	width  = width_from_pixels(pixels);
	height = height_from_width(width);

	printf("%" PRIu64 " bytes of l1 and l2 cache can hold ~%" PRIu64
	       " pixels. This means that if we want each frame to be held "
	       "entirely in l1 and l2 cache at a 16:9 ratio, the frame must be "
	       "shrunk to ~%dx%d\n",
	       cache, pixels, width, height);

	cache += features.l[3] / features.physical_cores;
	pixels = cache / 3;
	width  = width_from_pixels(pixels);
	height = height_from_width(width);

	printf("%" PRIu64 " bytes of l1, l2, and l3 cache can hold ~%" PRIu64
	       " pixels. This means that if we want each frame to be held "
	       "entirely in l1, l2, and l3 cache at a 16:9 ratio, the frame "
	       "must be shrunk to ~%dx%d\n",
	       cache, pixels, width, height);

	/* struct Video in_video; */
	/* open_video_file(&in_video, "vid.webm"); */

	return 0;
}
