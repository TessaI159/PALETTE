#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libavutil/imgutils.h>
#include <libavutil/timestamp.h>
#include <libswscale/swscale.h>

#include "Video.h"
#include "libavutil/avutil.h"

void open_video_file(struct Video *video, const char *filename) {
	memset(video, 0, sizeof(struct Video));

	if (avformat_open_input(&video->fmt_ctx, filename, NULL, NULL) < 0) {
		fprintf(stderr, "Could not open %s", filename);
		exit(1);
	}

	if (avformat_find_stream_info(video->fmt_ctx, NULL) < 0) {
		fprintf(stderr, "Could not find stream info.\n");
		exit(1);
	}

	video->video_stream_index = -1;
	video->audio_stream_index = -1;
	for (uint8_t i = 0; i < video->fmt_ctx->nb_streams; ++i) {
		if (video->fmt_ctx->streams[i]->codecpar->codec_type ==
		    AVMEDIA_TYPE_VIDEO) {
			video->video_stream_index = i;
		}
		if (video->fmt_ctx->streams[i]->codecpar->codec_type ==
		    AVMEDIA_TYPE_AUDIO) {
			video->audio_stream_index = i;
		}
	}

	if (video->video_stream_index == -1) {
		fprintf(stderr, "No video stream found.\n");
		exit(1);
	}

	if (video->audio_stream_index == -1) {
		fprintf(stderr,
			"No audio stream found. Continuing without audio.\n");
	}

	/* Format wrapper data */
	video->video_stream =
	    video->fmt_ctx->streams[video->video_stream_index];
	AVCodecParameters *par = video->video_stream->codecpar;

	video->w	 = par->width;
	video->h	 = par->height;
	video->duration	 = video->fmt_ctx->duration;
	video->nb_frames = video->video_stream->nb_frames;

	video->avg_frame_rate_num = video->video_stream->avg_frame_rate.num;
	video->avg_frame_rate_den = video->video_stream->avg_frame_rate.den;

	video->r_frame_rate_num = video->video_stream->r_frame_rate.num;
	video->r_frame_rate_den = video->video_stream->r_frame_rate.den;

	video->codec_name = avcodec_get_name(par->codec_id);

	video->decoder = avcodec_find_decoder(par->codec_id);

	if (!video->decoder) {
		fprintf(stderr, "Decoder not found.\n");
		exit(1);
	}

	video->cod_ctx = avcodec_alloc_context3(video->decoder);
	if (!video->cod_ctx) {
		fprintf(stderr, "Could not allocate codec context.\n");
		exit(1);
	}

	avcodec_parameters_to_context(video->cod_ctx, par);

	if (avcodec_open2(video->cod_ctx, video->decoder, NULL) < 0) {
		fprintf(stderr, "Could not open codec context.\n");
		exit(1);
	}
	video->time_base_num = video->video_stream->time_base.num;
	video->time_base_den = video->video_stream->time_base.den;
}
