cmake_minimum_required(VERSION 3.10)
project(gmrf_wind_mapping)
set(CMAKE_BUILD_TYPE "None") 



set(DEBUG OFF)

if(DEBUG)
  add_definitions(-DDEBUG)
  set(CMAKE_CXX_FLAGS "-g -O0 ${CMAKE_CXX_FLAGS}")
else()
  set(CMAKE_CXX_FLAGS "-O3 ${CMAKE_CXX_FLAGS}")
endif(DEBUG)


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
find_package(rosidl_default_generators REQUIRED) 

## System dependencies are found with CMake's conventions


rosidl_generate_interfaces(${PROJECT_NAME}
    srv/WindEstimation.srv
    DEPENDENCIES std_msgs
)
ament_export_dependencies(rosidl_default_runtime)



FILE(GLOB_RECURSE MYFILES_CPP "src/*.cpp")

include_directories(
  include
)

add_executable(gmrf_wind_mapping_node src/gmrf_node.cpp src/gmrf_map.cpp)

rosidl_get_typesupport_target(cpp_typesupport_target ${PROJECT_NAME} "rosidl_typesupport_cpp")
target_link_libraries(gmrf_wind_mapping_node "${cpp_typesupport_target}") 


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
)

install(
  TARGETS gmrf_wind_mapping_node
  DESTINATION lib/${PROJECT_NAME}
)

ament_package()