cmake_minimum_required(VERSION 3.16)
project(PALETTE C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

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
# find_package(OpenMP REQUIRED)
pkg_check_modules(AVFORMAT REQUIRED libavformat)
pkg_check_modules(AVCODEC REQUIRED libavcodec)
pkg_check_modules(AVUTIL REQUIRED libavutil)
pkg_check_modules(SWSCALE REQUIRED libswscale)

# ---- Main App Sources ----
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS src/*.c)
list(FILTER SOURCES EXCLUDE REGEX "src/test/.*\\.c$")

add_executable(PALETTE ${SOURCES})

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
  # OpenMP::OpenMP_C
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
  list(FILTER IMPL_SOURCES EXCLUDE REGEX "src/test/.*\\.c$")
  file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS src/test/*.c)

  add_executable(unit_tests ${TEST_SOURCES} ${IMPL_SOURCES})
  target_link_libraries(unit_tests unity)

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
    # OpenMP::OpenMP_C
    ${AVFORMAT_LIBRARIES}
    ${AVCODEC_LIBRARIES}
    ${AVUTIL_LIBRARIES}
    ${SWSCALE_LIBRARIES}
  )

  add_test(NAME all_tests COMMAND unit_tests)
endif()
