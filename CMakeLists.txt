cmake_minimum_required(VERSION 3.0)
project(TUTools)

enable_language(CXX)
set(CMAKE_CXX_STANDARD		14)
set(CMAKE_CXX_STANDARD_REQUIRED	ON)
set(CMAKE_CXX_EXTENSIONS	OFF)
set(CMAKE_CXX_FLAGS_RELEASE	"-Wall -O3")
set(CMAKE_CXX_FLAGS_DEBUG	"-g")
set(CMAKE_BUILD_TYPE		Release)
add_definitions(-DNDEBUG)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake_modules")
include(FindSIMD)
FindSIMD()

if(AVX2_FOUND)
  add_definitions(-DAVX2)
  set(CMAKE_CXX_FLAGS -mavx2)
elseif(AVX_FOUND)
  add_definitions(-DAVX)
  set(CMAKE_CXX_FLAGS -mavx)
elseif(SSE4_2_FOUND)
  add_definitions(-DSSE4)
  set(CMAKE_CXX_FLAGS "-msse4.2 -mpopcnt")
elseif(SSSE3_FOUND)
  add_definitions(-DSSSE3)
  set(CMAKE_CXX_FLAGS "-msse3 -mssse3")
elseif(SSE3_FOUND)
  add_definitions(-DSSE3)
  set(CMAKE_CXX_FLAGS -msse3)
elseif(SSE2_FOUND)
  add_definitions(-DSSE2)
  set(CMAKE_CXX_FLAGS -msse2)
elseif(NEON_FOUND)
  add_definitions(-DNEON)
  set(CMAKE_CXX_FLAGS "-mcpu=cortex-a7 -mtune=cortex-a7 -mfpu=neon -funsafe-math-optimizations")
endif()

if(APPLE)
  set(MACOSX_RPATH TRUE)
endif()

find_package(Boost REQUIRED)
include_directories(${PROJECT_SOURCE_DIR} ${Boost_INCLUDE_DIRS})

file(GLOB sources *.cc)
add_library(${PROJECT_NAME} SHARED ${sources})

install(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION lib)

add_subdirectory(TU)
add_subdirectory(cmake_modules)

option(BUILD_TUTOOLSPP_TEST_PROGRAMS "Building test programs for TUTools" OFF)

if(BUILD_TUTOOLSPP_TEST_PROGRAMS)
  add_subdirectory(Testprog)
endif()

