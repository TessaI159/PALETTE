#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Video.h"

static void get_video_info(struct Video *video, const char *file) {
	char command[512];

	snprintf(command, sizeof(command),
		 "ffprobe -v error -select_streams v:0 -show_entries "
		 "stream=width,height,r_frame_rate -of "
		 "default=noprint_wrappers=1:nokey=0 \"%s\"",
		 file);

	FILE *fp = popen(command, "r");
	if (!fp) {
		perror("popen");
		return;
	}

	char line[128];
	while (fgets(line, sizeof(line), fp)) {
		if (sscanf(line, "width=%d", &video->w) == 1)
			continue;
		if (sscanf(line, "height=%d", &video->h) == 1)
			continue;
		char numer[16];
		char denom[16];
		if (sscanf(line, "r_frame_rate=%15[^/]/%15s", numer, denom) ==
		    2)
			video->fps = atof(numer) / atof(denom);
	}
}

int get_next_frame(struct Video *video) {

	if (!video || !video->read_file || !video->frame) {
		return 0;
	}

	size_t frame_size = video->w * video->h * 3;
	size_t read_total = 0;
	while (read_total < frame_size) {
		size_t read_now =
		    fread(video->frame + read_total, 1, frame_size - read_total,
			  video->read_file);
		if (read_now == 0)
			return 0;
		read_total += read_now;
	}

	return 1;
}

void close_video_file(struct Video *video) {
	if (video->frame) {
		free(video->frame);
	}
	if (video->read_file) {
		pclose(video->read_file);
	}
	if (video->write_file) {
		pclose(video->write_file);
	}
}

void open_video_file(struct Video *video, const char *file, const char *mode) {
	/* usage: ffmpeg [options] [[infile options] -i infile]... {[outfile
	  options] outfile}... */
	memset(video, 0, sizeof(struct Video));
	get_video_info(video, file);

	char command[512];
	snprintf(command, sizeof(command),
		 "ffmpeg -v 1 -i \"%s\" -f rawvideo -pix_fmt rgb24 -", file);
	if (strstr(mode, "w") && strstr(mode, "r")) {
		printf("Cannot open file for both read and write. Exiting.\n");
		exit(1);
	}

	if (!(strstr(mode, "w") || strstr(mode, "r"))) {
		printf("Must specify either read or write mode. Exiting.\n");
		exit(1);
	}

	if (strstr(mode, "r")) {
		video->read_file = popen(command, mode);
	}

	if (strstr(mode, "w")) {
		video->write_file = popen(command, mode);
	}

	if (!video->read_file && !video->write_file) {
		perror("popen");
		return;
	}

	size_t frame_size = video->w * video->h * 3;
	video->frame	  = malloc(frame_size);
	if (!video->frame) {
		perror("malloc");
		video->frame = NULL;
		if (strstr(mode, "r")) {
			pclose(video->read_file);
		} else {
			pclose(video->write_file);
		}
		return;
	}
}
