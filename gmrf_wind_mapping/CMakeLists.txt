cmake_minimum_required(VERSION 3.10)
project(gmrf_wind_mapping)


set(CMAKE_BUILD_TYPE "Release") 
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 ${CMAKE_CXX_FLAGS}")

# DEPENDENCIES
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED) 
find_package(angles REQUIRED) 
find_package(Eigen3 REQUIRED) 
find_package(tf2_geometry_msgs REQUIRED) 
find_package(tf2_ros REQUIRED) 
find_package(std_msgs REQUIRED) 
find_package(nav_msgs REQUIRED) 
find_package(sensor_msgs REQUIRED) 
find_package(visualization_msgs REQUIRED) 
find_package(olfaction_msgs REQUIRED) 
find_package(geometry_msgs REQUIRED) 
find_package(gmrf_msgs REQUIRED) 
find_package(OpenCV REQUIRED)


FILE(GLOB_RECURSE MYFILES_CPP "src/*.cpp")

# INCLUDES
include_directories(
  include
  ${EIGEN3_INCLUDE_DIRS}
  ${OpenCV_INCLUDE_DIRS}
)

# Install launch files and related config files.
install(DIRECTORY
  launch
  DESTINATION share/${PROJECT_NAME}
)

# BUILDING
add_executable(gmrf_wind_mapping_node src/gmrf_node.cpp src/gmrf_map.cpp)
target_include_directories(gmrf_wind_mapping_node PUBLIC ${OpenCV_INCLUDE_DIRS})
target_link_libraries(gmrf_wind_mapping_node ${OpenCV_LIBS})


ament_target_dependencies(gmrf_wind_mapping_node
    rclcpp
    angles
    Eigen3
    tf2_ros
    tf2_geometry_msgs
    std_msgs
    nav_msgs
    sensor_msgs
    visualization_msgs
    olfaction_msgs
    geometry_msgs
    gmrf_msgs
)

install(
  TARGETS gmrf_wind_mapping_node
  DESTINATION lib/${PROJECT_NAME}
)

ament_package()