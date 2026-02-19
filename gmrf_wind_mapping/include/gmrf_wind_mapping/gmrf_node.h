/*********************************************************************
* gmrf_node.h - GMRF Wind-Distribution Mapping ROS2 Node
* Description: Implements the Gaussian Markov random field mapping algorithm for the
*              estimation of the windflow from a set of sparse 2D wind measurements.
* Authors: Javier Monroy, Pepe Ojeda 
*         University of Malaga (Spain)
*         http://mapir.isa.uma.es
*********************************************************************/

#include "nav_msgs/msg/occupancy_grid.hpp"
#include <sensor_msgs/msg/point_cloud2.hpp>
#include <visualization_msgs/msg/marker_array.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/math/constants/constants.hpp>
#include <tf2_ros/transform_listener.h>
#include <tf2_ros/buffer.h>
#include <angles/angles.h>
#include "olfaction_msgs/msg/anemometer.hpp"

// Core GMRF map
#include "gmrf_wind_core/gmrf_map.h"

// Services
#include "gmrf_msgs/srv/wind_estimation.hpp"
#include "gmrf_msgs/srv/add_wind_observation.hpp"
using WindEstimation = gmrf_msgs::srv::WindEstimation;
using AddWindObservation = gmrf_msgs::srv::AddWindObservation;

// Cgmrf class (ROS2 Node)
class Cgmrf : public rclcpp::Node
{
public:
    Cgmrf();
    ~Cgmrf();
    void update();
    void publishMaps();
    bool get_wind_value_srv(WindEstimation::Request::SharedPtr req, WindEstimation::Response::SharedPtr res);
    bool add_wind_observation_srv(AddWindObservation::Request::SharedPtr req, AddWindObservation::Response::SharedPtr res);
    double exec_freq;
    bool module_init;
    bool verbose;
    bool visualize_gmrf;
protected:
    void sensorCallback(const olfaction_msgs::msg::Anemometer::SharedPtr msg);
    void mapCallback(const nav_msgs::msg::OccupancyGrid::SharedPtr msg);
    inline void ReadMap();
    void initialize();

    // Subscriptions & Publishers
    rclcpp::Subscription<olfaction_msgs::msg::Anemometer>::SharedPtr sub_sensor;
    rclcpp::Subscription<nav_msgs::msg::OccupancyGrid>::SharedPtr occupancyMap_sub;
    rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr wind_array_pub;
    rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr wind_std_array_pub;

    // GMRF Maps
    std::unique_ptr<CGMRF_map> my_map;         // The GMRF Map being estimated
    nav_msgs::msg::OccupancyGrid occupancyMap; // Occupancy GridMap of the environment

    // Node Parameters
    std::string frame_id;       // frame where to plot the map, usually (map)
    std::string sensor_topic;   // Topic where the anemometer measurements are published
    std::string mapFilePath;    // Path to the map image file (grayscale) to use as occupancy map (if empty, use the one from map_topic)
    std::string map_topic;      // Topic where the occupancy gridmap is published    
    double cell_size;
    
    // GMRF parameters
    double GMRF_lambdaPrior_mass_conservation; // Weight for mass conservation law prior
    double GMRF_lambdaPrior_vorticity;         // Weight for vorticity prior
    double GMRF_lambdaPrior_obstacles;         // Weight for wind close to obstacles prior -->cells close to obstacles has only tangencial wind
    double GMRF_lambdaPrior_reg;               // Weight for regularization prior (should be very small)
    double observation_var_wind_speed;         // Variance of the wind speed measurement (m/s)^2
    double observation_var_wind_direction;     // Variance of the wind direction measurement (rad)^2
    
    // Variables
    tf2_ros::Buffer::SharedPtr tf_buffer;
    std::shared_ptr<tf2_ros::TransformListener> tf_listener;
    boost::mutex mutex_anemometer;
    boost::mutex mutex_position;
    double reading_speed;     // m/s
    double reading_direction; // rad
};