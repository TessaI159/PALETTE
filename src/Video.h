#ifndef VIDEO_H
#define VIDEO_H
#include <stdbool.h>
#include <stdint.h>

/* TODO (Tess): Document */

typedef struct AVFormatContext AVFormatContext;
typedef struct AVCodecContext  AVCodecContext;
typedef struct AVStream        AVStream;
typedef struct AVCodec         AVCodec;
typedef struct AVFrame         AVFrame;
typedef struct AVPacket        AVPacket;
typedef struct SwsContext      SwsContext;

struct Frame {
  uint8_t *rgb_data;
  uint32_t width;
  uint32_t height;
  int64_t  pts;
  bool     key_frame;
  uint32_t frame_index;
};

/* I call it a Video struct, but really it's just basic Video metadata as well
   as a single frame */
struct Video {
  AVFormatContext *fmt_ctx;
  AVCodecContext  *cod_ctx;
  AVStream        *video_stream;
  AVStream        *audio_stream;
  const AVCodec   *decoder;
  struct Frame     frame;
  int              video_stream_index;
  int              audio_stream_index;
  int              w;
  int              h;
  int64_t          duration;
  int64_t          nb_frames;
  const char      *codec_name;
  int16_t          avg_frame_rate_num;
  int16_t          avg_frame_rate_den;
  int16_t          r_frame_rate_num;
  int16_t          r_frame_rate_den;
  int16_t          time_base_num;
  int16_t          time_base_den;
};

void open_video_file(struct Video *video, const char *file);

#endif
