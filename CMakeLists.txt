cmake_minimum_required(VERSION 2.6)
project(vslab CXX)

set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

# check if the user provided CXXFLAGS on the command line
if (CMAKE_CXX_FLAGS)
    set(CXXFLAGS_PROVIDED true)
endif ()
if (CXXFLAGS_PROVIDED)
  set(CMAKE_CXX_FLAGS_DEBUG          "")
  set(CMAKE_CXX_FLAGS_MINSIZEREL     "")
  set(CMAKE_CXX_FLAGS_RELEASE        "")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "")
else ()
  set(IGNORED_WARNINGS "-Wno-gnu-folding-constant")
  set(CMAKE_CXX_FLAGS                "-std=c++11 -Wextra -Wall -fpermissive ${IGNORED_WARNINGS}")
  set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g")
  set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os")
  set(CMAKE_CXX_FLAGS_RELEASE        "-O4")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
endif ()

if(NOT CAF_ROOT_DIR)
  set(CAF_ROOT_DIR "/home/otto/vs/actor-framework/build/")
endif()
find_package(CAF COMPONENTS core io REQUIRED)
if (CAF_FOUND)
  include_directories(${CAF_INCLUDE_DIRS})
endif ()

# we need this to be able to include headers produced by uic in our code
# (CMAKE_BINARY_DIR holds a path to the build directory, while INCLUDE_DIRECTORIES() works just like INCLUDEPATH
# from qmake)
include_directories(. include  ${CMAKE_CURRENT_SOURCE_DIR}/third_party ${LIBCAF_INCLUDE_DIRS})

set(SOURCES
    src/main.cpp)

file(GLOB HEADERS "include/*.hpp")

add_executable(vslab ${SOURCES} ${HEADERS})
target_link_libraries(vslab ${CMAKE_DL_LIBS} ${CAF_LIBRARIES})
add_dependencies(vslab all_programs)

