# Mandatory information that CMAKE needs
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0)

# Project name
PROJECT (CMEA_2d-linFEM)

find_package (Threads REQUIRED)

# Prepare for Hunter (currently unused)
include("cmake/CommonHeader.cmake")

# Set some compiler flags
include ("cmake/Compiler.cmake")

# get Eigen
include("cmake/GetEigen.cmake")

# get IGL
ExternalProject_Add(
    IGL
    URL https://github.com/libigl/libigl/archive/refs/tags/v2.3.0.zip
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/IGL
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/IGL/include)
INCLUDE_DIRECTORIES(".")

FILE(GLOB SRC *.cpp *.hpp)
ADD_EXECUTABLE(fem2d ${SRC})
ADD_DEPENDENCIES(fem2d Eigen IGL)
ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
ADD_DEFINITIONS(-DCMEA_DATA_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data/")

TARGET_LINK_LIBRARIES(fem2d ${THREADS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

add_subdirectory(unittest)