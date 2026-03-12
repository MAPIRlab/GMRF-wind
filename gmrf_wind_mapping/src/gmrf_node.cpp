//========================================================================================
//	GMRF_node - Wind Estimation
//	Description: Implements the Gaussian Markov random field mapping algorithm for the
//               estimation of the windflow from a set of sparse 2D wind measurements.
//
//	topics subscribed: Anemometer measurements, Occupancy GridMap
//  topics published: - 
//	services: WindEstimation - service to get wind estimation at a given 2D location
//
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//	Revision log:
//	version: 2.0	20/11/2025
//========================================================================================

#include "gmrf_wind_mapping/gmrf_node.h"
#include "gmrf_wind_mapping/utils.h"

#include <geometry_msgs/msg/pose_stamped.hpp>
#include <geometry_msgs/msg/transform_stamped.hpp>

#include <chrono>
#include <tf2/time.h>
#include <yaml-cpp/yaml.h>
#include <cmath>     // Required for mathematical operations

using namespace std::placeholders;
using namespace gmrfw; 

Cgmrf::Cgmrf()
    : Node("GMRF_wind")
{
    printf("\n=================================================================");
    printf("\n=	        GMRF Wind-Distribution Mapping ROS2 Node              =");
    printf("\n=================================================================\n");
    
    // Load Parameters
    //------------------
    frame_id = declare_parameter<std::string>("frame_id", "map");
    sensor_topic = declare_parameter<std::string>("sensor_topic", "/anemometer");
    mapFilePath = declare_parameter<std::string>("map_yaml_file", "");
    map_topic = declare_parameter<std::string>("map_topic", "map");
    exec_freq = declare_parameter<double>("exec_freq", 2.0);
    cell_size = declare_parameter<double>("cell_size", 0.5);
    verbose = declare_parameter<bool>("verbose", false);
    visualize_gmrf = declare_parameter<bool>("visualize_gmrf", true);

    // Lambda/weights for the different priors
    GMRF_lambdaPrior_advection = declare_parameter<double>("GMRF_lambdaPrior_advection", 100.0);
    GMRF_lambdaPrior_mass_conservation = declare_parameter<double>("GMRF_lambdaPrior_mass_conservation", 100.0);
    GMRF_lambdaPrior_diffusion = declare_parameter<double>("GMRF_lambdaPrior_diffusion", 100.0);
    GMRF_lambdaPrior_obstacles = declare_parameter<double>("GMRF_lambdaPrior_obstacles", 10.0);

    // Observation variances (Gaussian likelihood)
    observation_var_wind_speed = declare_parameter<double>("observation_var_wind_speed", 0.0001);           // Variance of the wind speed measurement (m/s)^2
    observation_var_wind_direction = declare_parameter<double>("observation_var_wind_direction", 0.0001);       // Variance of the wind direction measurement (rad)^2

    
    // Subscriptions
    //----------------------------------
    sub_sensor = create_subscription<olfaction_msgs::msg::Anemometer>(sensor_topic, 10, std::bind(&Cgmrf::sensorCallback, this, _1));
    if (mapFilePath == "")
    {
        // Subscribe to MapServer topic only if we don't read a map from file
        occupancyMap_sub = create_subscription<nav_msgs::msg::OccupancyGrid>(map_topic, rclcpp::QoS(1).transient_local().reliable(), std::bind(&Cgmrf::mapCallback, this, _1));
        
        // Init state (we need to wait the Occupancy map to initialize the GMRF)
        module_init = false;
    }
    else
    {
        // Read map from file and Initialize GMRF
        ReadMap();
        initialize();
    }

    // Publishers
    //----------------------------------
    wind_array_pub = create_publisher<visualization_msgs::msg::MarkerArray>("wind_array_pub", 1);
    wind_std_array_pub = create_publisher<visualization_msgs::msg::Marker>("wind_std_array_pub", 1);
        
    // TF2 Listener
    tf_buffer = std::make_unique<tf2_ros::Buffer>(get_clock(), tf2::Duration(std::chrono::seconds(30)));
    tf_listener = std::make_shared<tf2_ros::TransformListener>(*tf_buffer);
}


Cgmrf::~Cgmrf()
{
}


//--------------------------
// CALLBACK - OCCUPANCY MAP
//--------------------------
void Cgmrf::mapCallback(const nav_msgs::msg::OccupancyGrid::SharedPtr msg)
{
    if (module_init)
        return;

    occupancyMap = *msg;
    RCLCPP_INFO(get_logger(), "[GMRF-node] Using OccupancyMap from topic '%s'", occupancyMap_sub->get_topic_name());

    initialize();
}

// Read Occupancy GridMap from a YAML file
inline void Cgmrf::ReadMap()
{
    // we can choose to read a map file directly from disk, if we don't want to use the one published by map_server
    // this is often useful when we need to alter the occupancy map to include inlets/outlets, which should be non-occupied
    // for GMRF but which may not be navigable (think of windows). 
    try
    {
        // Load YAML file
        YAML::Node yaml = YAML::LoadFile(mapFilePath);

        // if the image path is relative interpret it as relative to the YAML, not the working directory
        std::filesystem::path imagePath(yaml["image"].as<std::string>());
        if (imagePath.is_relative())
            imagePath = std::filesystem::path(mapFilePath).parent_path() / imagePath;

        // Load Data
        //-----------------------
        cv::Mat mapImage = cv::imread(imagePath, cv::IMREAD_GRAYSCALE);
        size_t width = mapImage.size().width;
        size_t height = mapImage.size().height;
        occupancyMap.data.resize(width * height);

        double free_thresh = 1 - yaml["free_thresh"].as<double>(); // Unused, this is for trinary maps. We only consider cells free or occupied, not "unknown"
        double occupied_thresh = 1 - yaml["occupied_thresh"].as<double>();

        for (int y = 0; y < height; y++)
            for (int x = 0; x < width; x++)
            {
                double value = mapImage.at<uint8_t>(y, x) / 255.;
                int8_t& occupancy = occupancyMap.data.at(width * (height - y - 1) + x);
                if (value > free_thresh)
                    occupancy = 0;
                else if (value < occupied_thresh)
                    occupancy = 100;
                else
                    occupancy = -1;
            }

        // Metadata
        //-------------------
        occupancyMap.header.frame_id = "map";
        occupancyMap.info.map_load_time = now();
        occupancyMap.info.resolution = yaml["resolution"].as<double>();
        auto origin_array = yaml["origin"].as<std::array<double, 3>>();
        occupancyMap.info.origin.position.x = origin_array[0];
        occupancyMap.info.origin.position.y = origin_array[1];
        occupancyMap.info.origin.orientation = Utils::createQuaternionMsgFromYaw(origin_array[2]);

        occupancyMap.info.height = height;
        occupancyMap.info.width = width;

        // Publish the actual map we are using, in case it is different from the one published by map_server
        static rclcpp::Publisher<nav_msgs::msg::OccupancyGrid>::SharedPtr map_republisher =
        create_publisher<nav_msgs::msg::OccupancyGrid>("gmrf_occupancy", rclcpp::QoS(1).transient_local());
        map_republisher->publish(occupancyMap);
    }   
    catch (const std::exception& e)
    {
        RCLCPP_ERROR(get_logger(), "Exception caught: '%s'", e.what());
    }
}



void Cgmrf::initialize()
{
    // Parse OccupancyGridmap to internal GMRF format
    TOccupancyMap occMap;
    occMap.data = occupancyMap.data;
    occMap.resolution = occupancyMap.info.resolution;      // cell size (m)
    occMap.width = occupancyMap.info.width;        // number of cells in X direction
    occMap.height = occupancyMap.info.height;       // number of cells in Y direction
    occMap.origin_x = occupancyMap.info.origin.position.x;        // world coordinates of the origin of the map
    occMap.origin_y = occupancyMap.info.origin.position.y;        // world coordinates of the origin of the map
    
    // Create the GMRF-Map and initialize its Prior Factors
    // Set only factors with a Lambda_prior > 0, (avoids adding unnecessary factors to the graph)
    CGMRF_map::Parameters params{
        cell_size,
        GMRF_lambdaPrior_advection,
        GMRF_lambdaPrior_mass_conservation,
        GMRF_lambdaPrior_diffusion,
        GMRF_lambdaPrior_obstacles,
    };
    my_map = std::make_unique<CGMRF_map>(occMap, 
                                        params,
                                        verbose,
                                        true // estimateTiming
                                        );
    my_map->update_lambdas(GMRF_lambdaPrior_advection, GMRF_lambdaPrior_mass_conservation, GMRF_lambdaPrior_diffusion, GMRF_lambdaPrior_obstacles);
    RCLCPP_INFO(get_logger(), "[GMRF-node] GMRF Initialized");
    module_init = true;


    // DEBUG - ADD SOME RANDOM OBSERVATIONS
    //      insertObservation_GMRF(double wind_speed, double wind_direction, double x_pos, double y_pos, double lambdaObs)
    //my_map->insertObservation_GMRF(12.0, 3.14, 1.15, 1.8, 100.0);
    //my_map->insertObservation_GMRF(12.0, 3.14, 1.15, 3.0, 100.0);
}


//----------------------------------
// CALLBACK - NEW SENSOR OBSERVATION
//----------------------------------
void Cgmrf::sensorCallback(const olfaction_msgs::msg::Anemometer::SharedPtr msg)
{
    // 1. Get wind measurement
    double downwind_direction_map;
    mutex_anemometer.lock();
    try
    {
        reading_speed = msg->wind_speed;         // (m/s)
        reading_direction = msg->wind_direction; // (rad) This is the Upwind direction with respect the Anemometer ref system (standard measurement)

        // RCLCPP_INFO(get_logger(), "Speed:%f Direction:%f", reading_speed, reading_direction);
        //  We need to transform this Upwind direction in the Anemometer ref system---- to ---- DownWind direction in the MAP ref system
        if (reading_speed != 0.0)
        {
            // Transform from anemometer ref_system to the map ref_system using TF
            // This requires the robot localization (tf tree) to be up to date
            geometry_msgs::msg::PoseStamped anemometer_upWind_pose, map_upWind_pose;
            try
            {
                anemometer_upWind_pose.header = msg->header;
                anemometer_upWind_pose.pose.position.x = 0.0;
                anemometer_upWind_pose.pose.position.y = 0.0;
                anemometer_upWind_pose.pose.position.z = 0.0;
                anemometer_upWind_pose.pose.orientation = Utils::createQuaternionMsgFromYaw(msg->wind_direction);

                // lookuptransform (target_frame, target_time, pose_in, fixed_frame, pose_out)
                tf_buffer->transform(anemometer_upWind_pose, map_upWind_pose, frame_id.c_str());

                downwind_direction_map = angles::normalize_angle(Utils::getYaw(map_upWind_pose.pose.orientation) + 3.14159);
            }
            catch (tf2::TransformException& ex)
            {
                RCLCPP_ERROR(get_logger(), "[GMRF-node] - %s - Error: %s", __FUNCTION__, ex.what());
            }
        }
        else
        {
            downwind_direction_map = 0.0;
        }
    }
    catch (std::exception e)
    {
        RCLCPP_ERROR(get_logger(), "[GMRF-node] Exception at new Obs: %s ", e.what());
    }
    mutex_anemometer.unlock();
    // RCLCPP_INFO(get_logger(), "[GMRF-node] New wind observation! %.2f m/s  %.2f rad (DownWind in the map ref system)",msg->wind_speed,
    // downwind_direction_map);


    // 2. Get pose of the sensor in the map reference system
    geometry_msgs::msg::TransformStamped transform;
    bool know_sensor_pose = true;
    try
    {
        // lookuptransform (target_frame, source_frame, result_tf)
        transform = tf_buffer->lookupTransform(frame_id.c_str(), msg->header.frame_id.c_str(), msg->header.stamp);
    }
    catch (tf2::TransformException ex)
    {
        RCLCPP_ERROR(get_logger(), "[GMRF] exception when reading observation: %s", ex.what());
        know_sensor_pose = false;
    }

    // 3. Add observation to the GMRF
    if (module_init)
    {
        // Current sensor pose in the map
        float x_pos = transform.transform.translation.x;
        float y_pos = transform.transform.translation.y;

        mutex_anemometer.lock();
        // TODO: Read from parameter or add to the message the variance of the measurement
        double var_wind_speed = observation_var_wind_speed; // (m/s)^2
        double var_wind_direction = observation_var_wind_direction; // (rad)^2
        my_map->insertObservation_GMRF(reading_speed, downwind_direction_map, var_wind_speed, var_wind_direction, x_pos, y_pos);
        mutex_anemometer.unlock();
    }
}


void Cgmrf::publishMaps()
{
    if (!visualize_gmrf)
        return;
    
    visualization_msgs::msg::MarkerArray wind_array;
    visualization_msgs::msg::Marker wind_std_array;
    Utils::createWindMarkerArrayFromGMRF(*my_map, frame_id, wind_array, wind_std_array);
    wind_array_pub->publish(wind_array);
    wind_std_array_pub->publish(wind_std_array);
}


void Cgmrf::update()
{
    // Update Lambda parameters (from parameter server)
    GMRF_lambdaPrior_advection = get_parameter("GMRF_lambdaPrior_advection").as_double();
    GMRF_lambdaPrior_mass_conservation = get_parameter("GMRF_lambdaPrior_mass_conservation").as_double();
    GMRF_lambdaPrior_diffusion = get_parameter("GMRF_lambdaPrior_diffusion").as_double();
    GMRF_lambdaPrior_obstacles = get_parameter("GMRF_lambdaPrior_obstacles").as_double();
    my_map->update_lambdas(GMRF_lambdaPrior_advection, GMRF_lambdaPrior_mass_conservation, GMRF_lambdaPrior_diffusion, GMRF_lambdaPrior_obstacles);

    // Update GMRF estimation
     my_map->MAP_estimation_GMRF();
}


bool Cgmrf::get_wind_value_srv(WindEstimation::Request::SharedPtr req, WindEstimation::Response::SharedPtr res)
{
    if (!module_init)
    {
        RCLCPP_ERROR(get_logger(), "Trying to query GMRF wind, but it is not initialized yet (probably has not received the occupancy map)");
        return false;
    }

    Eigen::Vector2i dimensions = my_map->map_size();
    res->map_width = dimensions.x();
    // an empty request means get all the points
    if (req->x.empty())
    {
        size_t size = dimensions.x() * dimensions.y();

        res->u.reserve(size);
        res->v.reserve(size);
        res->var_u.reserve(size);
        res->var_v.reserve(size);
        res->cov_uv.reserve(size);

        for (int i = 0; i < size; i++)
        {
            WindVector r = my_map->getEstimation(i);
            res->u.push_back(r.x);
            res->v.push_back(r.y);
            res->var_u.push_back(r.varX);
            res->var_v.push_back(r.varY);
            res->cov_uv.push_back(r.covXY);
        }

        return true;
    }

    // Since the wind fields are identical among different instances, return just the information from instance[0]
    for (int i = 0; i < req->x.size(); i++)
    {
        WindVector r = my_map->getEstimation(req->x[i], req->y[i]);
        res->u.push_back(r.x);
        res->v.push_back(r.y);
        res->var_u.push_back(r.varX);
        res->var_v.push_back(r.varY);
        res->cov_uv.push_back(r.covXY);
    }
    return true;
}

bool Cgmrf::add_wind_observation_srv(AddWindObservation::Request::SharedPtr req, AddWindObservation::Response::SharedPtr res)
{
    if (!module_init)
    {
        RCLCPP_ERROR(get_logger(), "Trying to add GMRF wind observation, but it is not initialized yet (probably has not received the occupancy map)");
        return false;
    }

    for (int i = 0; i < req->wind_speed.size(); i++)
    {
        // TODO: Read from parameter or add to the message the variance of the measurement
        double var_wind_speed = 0.001; // (m/s)^2
        double var_wind_direction = 0.0001; // (rad)^2
        my_map->insertObservation_GMRF(req->wind_speed[i], req->wind_direction[i], var_wind_speed, var_wind_direction, req->x_pos[i], req->y_pos[i]);
        if (verbose)
            RCLCPP_INFO(get_logger(), "[GMRF-node] New wind observation added via service at (%.2f, %.2f): %.2f m/s, %.2f rad", req->x_pos[i], req->y_pos[i], req->wind_speed[i], req->wind_direction[i]);
    }
    return true;
}


//-----------------------------------------------------------------------------
//                                    MAIN
//----------------------------------------------------------------------------
int main(int argc, char** argv)
{
    // Initialize ROS2
    rclcpp::init(argc, argv);

    // Create the GMRF-wind node
    auto my_gmrf_map = std::make_shared<Cgmrf>();

    // Create and Offer the AddWindObservation service
    auto add_wind_observation_service = my_gmrf_map->create_service<AddWindObservation>("AddWindObservation", std::bind(&Cgmrf::add_wind_observation_srv, my_gmrf_map.get(), _1, _2));

    // Create and Offer the WindEstimation service
    auto service = my_gmrf_map->create_service<WindEstimation>("WindEstimation", std::bind(&Cgmrf::get_wind_value_srv, my_gmrf_map.get(), _1, _2));

    // Main loop
    RCLCPP_INFO(my_gmrf_map->get_logger(), "[gmrf-wind] MAIN LOOP....");
    rclcpp::Time last_publication_time = my_gmrf_map->now();
    rclcpp::Rate loop_rate(my_gmrf_map->exec_freq);

    while (rclcpp::ok())
    {
        rclcpp::spin_some(my_gmrf_map);     // Callbacks & Services

        if (my_gmrf_map->module_init)
        {
            // Update Estimation
            my_gmrf_map->update();
           
            // Publish Map as markers (RVIZ2)
            my_gmrf_map->publishMaps();

            // Info about actual rate
            if (my_gmrf_map->verbose)
                RCLCPP_INFO(my_gmrf_map->get_logger(), "[gmrf] Updating every %f seconds. Intended preiod was %f",
                            (my_gmrf_map->now() - last_publication_time).seconds(), 1.0 / my_gmrf_map->exec_freq);
            last_publication_time = my_gmrf_map->now();
        }
        else
        {
            if (my_gmrf_map->verbose)
                RCLCPP_INFO(my_gmrf_map->get_logger(), "[gmrf] Waiting for initialization (Occupancy map of the environment).");
        }

        // Keep the loop rate
        loop_rate.sleep();
    }
}
