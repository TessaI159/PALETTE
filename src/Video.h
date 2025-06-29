#ifndef VIDEO_H
#define VIDEO_H
#include <stdint.h>
#include <stdio.h>

struct Video {
	FILE *read_file;
	FILE *write_file;

	uint8_t *frame;
	int	 w;
	int	 h;
	float	 fps;
};

/* Keep video file names under ~350 characters, please */
void open_video_file(struct Video *video, const char *file,
			const char *mode);
int  get_next_frame(struct Video *video);
void close_video_file(struct Video *video);

#endif
