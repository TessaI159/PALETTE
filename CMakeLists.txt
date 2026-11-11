cmake_minimum_required(VERSION 3.16)
project(PALETTE C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Detect if cross-compiling to Windows
if(CMAKE_TOOLCHAIN_FILE)
  message(STATUS "Using toolchain: ${CMAKE_TOOLCHAIN_FILE}")
endif()

# ---- Compile Options ----
add_compile_options(
  -g
  -Wall
  -Wextra
  -Werror
  -Wno-error=unused-variable
  -Wno-error=unused-function
  -Wno-error=unused-parameter
)

# ---- Option for Tests ----
option(PALETTE_BUILD_TESTS "Build unit tests" ON)

# ---- Optional FFmpeg Dependency ----
find_package(PkgConfig QUIET)
set(FFMPEG_FOUND FALSE)

if(PkgConfig_FOUND AND NOT WIN32)
  pkg_check_modules(AVFORMAT libavformat)
  pkg_check_modules(AVCODEC libavcodec)
  pkg_check_modules(AVUTIL libavutil)
  pkg_check_modules(SWSCALE libswscale)
  if(AVFORMAT_FOUND AND AVCODEC_FOUND AND AVUTIL_FOUND AND SWSCALE_FOUND)
    set(FFMPEG_FOUND TRUE)
  endif()
elseif(WIN32)
  set(FFMPEG_WIN_ROOT "$ENV{HOME}/.local/share/ffmpeg-windows-install")

  if(EXISTS "${FFMPEG_WIN_ROOT}/include/libavformat/avformat.h")
    set(FFMPEG_FOUND TRUE)
    set(AVFORMAT_INCLUDE_DIRS "${FFMPEG_WIN_ROOT}/include")
    set(AVFORMAT_LIBRARY_DIRS "${FFMPEG_WIN_ROOT}/lib")
    set(AVFORMAT_LIBRARIES avformat)
    set(AVCODEC_INCLUDE_DIRS "${FFMPEG_WIN_ROOT}/include")
    set(AVCODEC_LIBRARY_DIRS "${FFMPEG_WIN_ROOT}/lib")
    set(AVCODEC_LIBRARIES avcodec)
    set(AVUTIL_INCLUDE_DIRS "${FFMPEG_WIN_ROOT}/include")
    set(AVUTIL_LIBRARY_DIRS "${FFMPEG_WIN_ROOT}/lib")
    set(AVUTIL_LIBRARIES avutil)
    set(SWSCALE_INCLUDE_DIRS "${FFMPEG_WIN_ROOT}/include")
    set(SWSCALE_LIBRARY_DIRS "${FFMPEG_WIN_ROOT}/lib")
    set(SWSCALE_LIBRARIES swscale)
  endif()
endif()

if(FFMPEG_FOUND)
  message(STATUS "FFmpeg libraries found. Video decoding support enabled.")
else()
  message(STATUS "FFmpeg libraries not found. Building PALETTE core engine without video decoding.")
endif()

# ---- Architecture Flags (set once; applies to every target using these files) ----
# File names match the on-disk case exactly (Linux is case-sensitive).
set(AVX_SOURCES
  src/Average_AVX.c
  src/Difference_AVX.c
  src/conversion_AVX.c
)
set(SSE_SOURCES
  src/Average_SSE.c
  src/Difference_SSE.c
  src/Conversion_SSE.c
  src/Math_SSE.c
)

if(MSVC)
  set_source_files_properties(${AVX_SOURCES} PROPERTIES COMPILE_FLAGS "/arch:AVX2")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|i[3-6]86")
  set_source_files_properties(${AVX_SOURCES} PROPERTIES COMPILE_FLAGS "-mavx2 -mfma")
  set_source_files_properties(${SSE_SOURCES} PROPERTIES COMPILE_FLAGS "-msse4.1")
endif()

# ---- Main App Sources ----
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS src/*.c)
list(FILTER SOURCES EXCLUDE REGEX "/\\.#[^/]*$")
list(FILTER IMPL_SOURCES EXCLUDE REGEX "/\\.#[^/]*$")
list(FILTER SOURCES EXCLUDE REGEX "/src/test/.*\\.c$")
list(FILTER SOURCES EXCLUDE REGEX "/src/trig_test\\.c$")
list(FILTER SOURCES EXCLUDE REGEX "/src/pow_speed_test\\.c$")



add_executable(PALETTE ${SOURCES})

target_include_directories(PALETTE PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(PALETTE m)

if(FFMPEG_FOUND)
  target_compile_definitions(PALETTE PRIVATE PALETTE_HAS_FFMPEG=1)
  target_include_directories(PALETTE PRIVATE
    ${AVFORMAT_INCLUDE_DIRS}
    ${AVCODEC_INCLUDE_DIRS}
    ${AVUTIL_INCLUDE_DIRS}
    ${SWSCALE_INCLUDE_DIRS}
  )
  target_link_directories(PALETTE PRIVATE
    ${AVFORMAT_LIBRARY_DIRS}
    ${AVCODEC_LIBRARY_DIRS}
    ${AVUTIL_LIBRARY_DIRS}
    ${SWSCALE_LIBRARY_DIRS}
  )
  target_link_libraries(PALETTE
    ${AVFORMAT_LIBRARIES}
    ${AVCODEC_LIBRARIES}
    ${AVUTIL_LIBRARIES}
    ${SWSCALE_LIBRARIES}
  )
endif()

# ---- Unit Tests ----
if(PALETTE_BUILD_TESTS)
  enable_testing()

  set(UNITY_DIR ${PROJECT_SOURCE_DIR}/third_party/Unity)
  add_library(unity STATIC ${UNITY_DIR}/unity.c)
  target_include_directories(unity PUBLIC ${UNITY_DIR})
  target_compile_definitions(unity PRIVATE UNITY_INCLUDE_DOUBLE)

  # Implementation sources: everything except the app's main() and test/bench files.
  file(GLOB_RECURSE IMPL_SOURCES CONFIGURE_DEPENDS src/*.c)
  list(FILTER SOURCES EXCLUDE REGEX "/\\.#[^/]*$")
  list(FILTER IMPL_SOURCES EXCLUDE REGEX "/\\.#[^/]*$")
  list(FILTER IMPL_SOURCES EXCLUDE REGEX "/src/trig_test\\.c$")
  list(FILTER IMPL_SOURCES EXCLUDE REGEX "/src/pow_speed_test\\.c$")
  list(FILTER IMPL_SOURCES EXCLUDE REGEX "/src/main\\.c$")
  list(FILTER IMPL_SOURCES EXCLUDE REGEX "/src/test/.*\\.c$")


  file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS src/test/*.c)
  list(FILTER SOURCES EXCLUDE REGEX "/\\.#[^/]*$")
  list(FILTER IMPL_SOURCES EXCLUDE REGEX "/\\.#[^/]*$")
  message(STATUS "Test sources: ${TEST_SOURCES}")

  add_executable(unit_tests ${TEST_SOURCES} ${IMPL_SOURCES})

  target_compile_definitions(unit_tests PRIVATE
    $<$<CONFIG:Debug>:PALETTE_DEBUG>)

  target_include_directories(unit_tests PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${UNITY_DIR}
  )

  target_link_libraries(unit_tests m unity)

  if(FFMPEG_FOUND)
    target_compile_definitions(unit_tests PRIVATE PALETTE_HAS_FFMPEG=1)
    target_include_directories(unit_tests PRIVATE
      ${AVFORMAT_INCLUDE_DIRS}
      ${AVCODEC_INCLUDE_DIRS}
      ${AVUTIL_INCLUDE_DIRS}
      ${SWSCALE_INCLUDE_DIRS}
    )
    target_link_directories(unit_tests PRIVATE
      ${AVFORMAT_LIBRARY_DIRS}
      ${AVCODEC_LIBRARY_DIRS}
      ${AVUTIL_LIBRARY_DIRS}
      ${SWSCALE_LIBRARY_DIRS}
    )
    target_link_libraries(unit_tests
      ${AVFORMAT_LIBRARIES}
      ${AVCODEC_LIBRARIES}
      ${AVUTIL_LIBRARIES}
      ${SWSCALE_LIBRARIES}
    )
  endif()

  add_test(NAME all_tests COMMAND unit_tests)
endif()
