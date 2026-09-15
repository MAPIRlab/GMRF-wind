cmake_minimum_required(VERSION 3.10)
project(gmrf_wind_mapping)


set(CMAKE_BUILD_TYPE "Release") 
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 ${CMAKE_CXX_FLAGS}")

# DEPENDENCIES
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED) 
find_package(angles REQUIRED) 
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
find_package(yaml-cpp REQUIRED)
find_package(Ceres REQUIRED)

# INCLUDES
include_directories(
  include
  ${OpenCV_INCLUDE_DIRS}
  ${CERES_INCLUDE_DIRS}
)

# The GMRF-WIND-CORE library
add_subdirectory(gmrf_wind_core)

# Install launch files and related config files.
install(DIRECTORY
  launch
  IEA_Annex_20_dataset
  DESTINATION share/${PROJECT_NAME}
)

# GMRF-NODE
add_executable(gmrf_wind_mapping_node src/gmrf_node.cpp)
target_include_directories(gmrf_wind_mapping_node PUBLIC ${OpenCV_INCLUDE_DIRS})
target_link_libraries(gmrf_wind_mapping_node gmrf_wind_core ${OpenCV_LIBS} yaml-cpp)

# GMRF-VALIDATION
add_executable(gmrf_validation src/gmrf_validation.cpp)
target_include_directories(gmrf_validation PUBLIC ${OpenCV_INCLUDE_DIRS} ${CERES_INCLUDE_DIRS})
target_link_libraries(gmrf_validation gmrf_wind_core ${OpenCV_LIBS} ${CERES_LIBRARIES} yaml-cpp)

# DEPENDENCIES
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

ament_target_dependencies(gmrf_validation
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

# INSTALLATION
install(
  TARGETS 
  gmrf_wind_mapping_node
  gmrf_validation
  DESTINATION lib/${PROJECT_NAME}
)

install(
  TARGETS 
  gmrf_wind_core
  DESTINATION lib
)

install(
  DIRECTORY
  gmrf_wind_core/include/gmrf_wind_core
  DESTINATION include
)

ament_export_include_directories(include)
ament_export_libraries(gmrf_wind_core)
ament_package()