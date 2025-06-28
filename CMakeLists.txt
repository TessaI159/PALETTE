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

# Include local headers
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR})
if(${APPLE})
include_directories("/usr/local/include")
endif()
find_package(SDL3 REQUIRED)

add_executable(PALETTE ${SOURCES})

target_link_libraries(PALETTE m SDL3)
