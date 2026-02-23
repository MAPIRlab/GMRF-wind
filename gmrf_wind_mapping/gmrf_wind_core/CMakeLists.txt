project(gmrf_wind_core)

if(NOT ${CMAKE_BUILD_TYPE})
    set(CMAKE_BUILD_TYPE "Release") 
endif()

set(CMAKE_CXX_FLAGS "-fPIC")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 ${CMAKE_CXX_FLAGS}")

find_package(Eigen3 REQUIRED) 

add_library(gmrf_wind_core SHARED src/gmrf_map.cpp)

target_include_directories(gmrf_wind_core PUBLIC 
    include
    ${EIGEN3_INCLUDE_DIRS}
)