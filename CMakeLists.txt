cmake_minimum_required(VERSION 3.10)

project(v_umlaut VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(Eigen3 REQUIRED)

add_definitions(-D_USE_MATH_DEFINES)

include_directories(${PROJECT_SOURCE_DIR}/external)
include_directories(${PROJECT_SOURCE_DIR})

add_executable(v_umlaut
    main.cpp
    Utils.cpp
    SyntheticData.cpp
    Solver.cpp
    external/poselib/univariate.cc
    external/poselib/relpose_5pt.cc
    external/poselib/relpose_7pt.cc
)

target_link_libraries(v_umlaut PUBLIC Eigen3::Eigen)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")