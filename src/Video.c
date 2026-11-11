#include "Video.h"
#include "Conversion.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libavutil/imgutils.h>
#include <libavutil/timestamp.h>
#include <libswscale/swscale.h>

bool video_open(struct Video *video, const char *filename) {
  if (!video || !filename)
    return false;

  memset(video, 0, sizeof(struct Video));

  if (avformat_open_input(&video->fmt_ctx, filename, NULL, NULL) < 0) {
    fprintf(stderr, "Error: Could not open video file %s\n", filename);
    return false;
  }

  if (avformat_find_stream_info(video->fmt_ctx, NULL) < 0) {
    fprintf(stderr, "Error: could not find stream info for %s\n", filename);
    return false;
  }

  video->video_stream_index = -1;
  video->audio_stream_index = -1;
  for (unsigned int i = 0; i < video->fmt_ctx->nb_streams; ++i) {
    if (video->fmt_ctx->streams[i]->codecpar->codec_type ==
            AVMEDIA_TYPE_VIDEO &&
        video->video_stream_index == -1) {
      video->video_stream_index = (int)i;
    }
    if (video->fmt_ctx->streams[i]->codecpar->codec_type ==
            AVMEDIA_TYPE_AUDIO &&
        video->audio_stream_index == -1) {
      video->audio_stream_index = (int)i;
    }
  }

  if (video->video_stream_index == -1) {
    fprintf(stderr, "Error: No video stream found in %s.\n", filename);
    return false;
  }

  if (video->audio_stream_index == -1) {
    fprintf(stderr,
            "Warning: No audio stream found. Continuing without audio.\n");
  }

  AVCodecParameters *par = video->video_stream->codecpar;

  video->w = par->width;
  video->h = par->height;
  video->duration = video->fmt_ctx->duration;
  video->nb_frames = video->video_stream->nb_frames;

  video->avg_frame_rate_num = (int16_t)video->video_stream->avg_frame_rate.num;
  video->avg_frame_rate_den = (int16_t)video->video_stream->avg_frame_rate.den;

  video->r_frame_rate_num = (int16_t)video->video_stream->r_frame_rate.num;
  video->r_frame_rate_den = (int16_t)video->video_stream->r_frame_rate.den;

  video->codec_name = avcodec_get_name(par->codec_id);
  video->decoder = avcodec_find_decoder(par->codec_id);

  if (!video->decoder) {
    fprintf(stderr, "Error: Decoder not found for codec %s.\n",
            video->codec_name);
    video_close(video);
    return false;
  }

  video->cod_ctx = avcodec_alloc_context3(video->decoder);
  if (!video->cod_ctx) {
    fprintf(stderr, "Error: Could not allocate codec context.\n");
    video_close(video);
    return false;
  }

  if (avcodec_open2(video->cod_ctx, video->decoder, NULL) < 0) {
    fprintf(stderr, "Error: Could not open codec context.\n");
    return false;
  }
  video->time_base_num = (int16_t)video->video_stream->time_base.num;
  video->time_base_den = (int16_t)video->video_stream->time_base.den;

  video->raw_frame = av_frame_alloc();
  video->rgb_frame = av_frame_alloc();
  video->packet = av_frame_alloc();

  if (!video->raw_frame || !video->rgb_frame || !video->packet) {
    fprintf(stderr, "Error: Failed to allocate frames or packet.\n");
    video_close(video);
    return false;
  }

  int buf_size =
      av_image_alloc(video->rgb_frame->data, video->rgb_frame->linesize,
                     video->w, video->h, AV_PIX_FMT_RGB24, 32);
  if (buf_size < 0) {
    fprintf(stderr, "Error: Could not allocate RGB image buffer.\n");
    video_close(video);
    return false;
  }

  video->rgb_buffer = video->rgb_frame->data[0];

  video->sws_cts = sws_getContext(video->w, video->h, video->cod_cts->fix_fmt,
                                  video->w, video->h, AV_PIX_FMT_RGB24,
                                  SWS_BILINEAR, NULL, NULL, NULL);
  if (!video->sws_ctx) {
    fprintf(stderr, "Error: Could not initialize sws converter context.\n");
    video_close(video);
    return false;
  }

  video->current_frame = 0;
  video->flushing = false;
  return true;
}

void video_close(struct Video *video) {
  if (!video)
    return;
  if (video->sws_ctx) {
    sws_freeContext(video->sws_ctx);
    video->sws_context = NULL;
  }
  if (video->rgb_buffer) {
    av_freep(&video->rgb_frame->data[0]);
    video->rgb_buffer = NULL;
  }
  if (video->rgb_frame) {
    av_frame_free(&video->rgb_frame);
  }
  if (video->raw_frame) {
    av_frame_free(&video->raw_frame);
  }
  if (video->packet) {
    av_packet_free(&video->packet);
  }
  if (video->cod_ctx) {
    avcodec_free_context(&video->cod_ctx);
  }
  if (video->fmt_ctx) {
    avformat_close_input(&video->fmt_ctx);
  }
  memset(video, 0, sizeof(struct Video));
}
