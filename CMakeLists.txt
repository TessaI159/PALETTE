cmake_minimum_required(VERSION 3.16)
project(PALETTE C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

add_compile_options(
  -g
  -Wall
  -Wextra
  -Werror
  -Wno-error=unused-variable
  -Wno-error=unused-function
  -Wno-error=unused-parameter
)

# Automatically check for new source files
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS src/*.c)


find_package(PkgConfig REQUIRED)
pkg_check_modules(AVFORMAT REQUIRED libavformat)
pkg_check_modules(AVCODEC REQUIRED libavcodec)
pkg_check_modules(AVUTIL REQUIRED libavutil)
pkg_check_modules(SWSCALE REQUIRED libswscale)

# Include local headers
include_directories(
  ${PROJECT_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}
  ${AVFORMAT_INCLUDE_DIRS}
  ${AVCODEC_INCLUDE_DIRS}
  ${AVUTIL_INCLUDE_DIRS}
  ${SWSCALE_INCLUDE_DIRS}
)

link_directories(
  ${AVFORMAT_LIBRARY_DIRS}
  ${AVCODEC_LIBRARY_DIRS}
  ${AVUTIL_LIBRARY_DIRS}
  ${SWSCALE_LIBRARY_DIRS}
)


add_executable(PALETTE ${SOURCES})

target_link_libraries(
  PALETTE
  m
  ${AVFORMAT_LIBRARIES}
  ${AVCODEC_LIBRARIES}
  ${AVUTIL_LIBRARIES}
  ${SWSCALE_LIBRARIES}
)
