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
option(PALETTE_BUILD_TESTS "Build unit tests" OFF)

# ---- Dependencies ----
find_package(PkgConfig REQUIRED)
if(NOT WIN32)
  pkg_check_modules(AVFORMAT REQUIRED libavformat)
  pkg_check_modules(AVCODEC REQUIRED libavcodec)
  pkg_check_modules(AVUTIL REQUIRED libavutil)
  pkg_check_modules(SWSCALE REQUIRED libswscale)
else()
  # Manual paths for cross-compiled FFmpeg (adjust as needed)
  set(FFMPEG_WIN_ROOT "$ENV{HOME}/.local/share/ffmpeg-windows-install")
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

# ---- Main App Sources ----
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS src/*.c)
list(FILTER SOURCES EXCLUDE REGEX "src/test/.*\.c$")

add_executable(PALETTE ${SOURCES})

if(MSVC)
  set_source_files_properties(src/Average_AVX.c PROPERTIES COMPILE_FLAGS
    "/arch:AVX2")
  
else()
  set_source_files_properties(src/FeatureDetection.c PROPERTIES COMPILE_FLAGS
  "-mxsave")

  
  set_source_files_properties(src/Average_AVX.c PROPERTIES COMPILE_FLAGS
    "-mavx2")
  
  set_source_files_properties(src/Average_SSE3.c PROPERTIES COMPILE_FLAGS
    "-msse2")
  
endif()

# ---- Includes and Linking ----
target_include_directories(PALETTE PRIVATE
  ${PROJECT_SOURCE_DIR}/src
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
  m
  ${AVFORMAT_LIBRARIES}
  ${AVCODEC_LIBRARIES}
  ${AVUTIL_LIBRARIES}
  ${SWSCALE_LIBRARIES}
)

# ---- Unit Tests ----
if (PALETTE_BUILD_TESTS)
  enable_testing()

  set(UNITY_DIR ${PROJECT_SOURCE_DIR}/third_party/Unity)
  add_library(unity STATIC ${UNITY_DIR}/unity.c)
  target_include_directories(unity PUBLIC ${UNITY_DIR})

  file(GLOB_RECURSE IMPL_SOURCES CONFIGURE_DEPENDS src/*.c)
  list(FILTER IMPL_SOURCES EXCLUDE REGEX "src/main.c")
  list(FILTER IMPL_SOURCES EXCLUDE REGEX "src/test/.*\.c$")
  file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS src/test/*.c)

  add_executable(unit_tests ${TEST_SOURCES} ${IMPL_SOURCES})
  target_compile_options(unit_tests PRIVATE -fsanitize=address,undefined
    -fno-omit-frame-pointer)
  
  target_link_options(unit_tests PRIVATE -fsanitize=address,undefined)

  target_compile_definitions(unity PRIVATE UNITY_INCLUDE_DOUBLE)
  target_compile_definitions(unit_tests PRIVATE
    $<$<CONFIG:Debug>:PALETTE_DEBUG>)
  

  if (MSVC)
    set_source_files_properties(src/Average_AVX.c PROPERTIES COMPILE_FLAGS
      "/arch:AVX")
    
  else()
    set_source_files_properties(src/Average_AVX.c PROPERTIES COMPILE_FLAGS
      "-mavx")
    
    set_source_files_properties(src/Average_SSE3.c PROPERTIES COMPILE_FLAGS
      "-msse3")
    
  endif()

  target_include_directories(unit_tests PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${UNITY_DIR}
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
    m
    unity
    ${AVFORMAT_LIBRARIES}
    ${AVCODEC_LIBRARIES}
    ${AVUTIL_LIBRARIES}
    ${SWSCALE_LIBRARIES}
  )

  add_test(NAME all_tests COMMAND unit_tests)
endif()
