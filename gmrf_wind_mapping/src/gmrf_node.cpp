//========================================================================================
//	GMRF_node - Wind Estimation
//	Description: Implements the Gaussian Markov random field mapping algorithm for the
//               estimation of the windflow from a set of sparse 2D wind measurements.
//
//	topics subscribed:
//  topics published:
//	services:
//
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//	Revision log:
//	version: 1.0	23/02/2017
//========================================================================================

#include "gmrf_node.h"
#include <geometry_msgs/msg/pose_stamped.hpp>
#include <geometry_msgs/msg/transform_stamped.hpp>
#include "Utils.h"

using namespace std::placeholders;

Cgmrf::Cgmrf() : Node("GMRF_wind")
{
	printf("\n=================================================================");
	printf("\n=	             GMRF Wind-Distribution Mapping Node              =");
	printf("\n=================================================================\n");

	//------------------
	// Load Parameters
	//------------------
	frame_id = declare_parameter<std::string>("frame_id", "map");
	sensor_topic = declare_parameter<std::string>("sensor_topic", "/anemometer");
	exec_freq = declare_parameter<double>("exec_freq", 2.0);
	cell_size = declare_parameter<double>("cell_size", 0.5);

	GMRF_lambdaPrior_reg = declare_parameter<double>("GMRF_lambdaPrior_reg", 1);                                 // Weight for regularization prior -> neighbour cells have similar wind vectors
	GMRF_lambdaPrior_mass_conservation = declare_parameter<double>("GMRF_lambdaPrior_mass_conservation", 10000); // Weight for mass conservation law prior
	GMRF_lambdaPrior_obstacles = declare_parameter<double>("GMRF_lambdaPrior_obstacles", 10);                    // Weight for wind close to obstacles prior -->cells close to obstacles has only tangencial wind
	GMRF_lambdaObs = declare_parameter<double>("GMRF_lambdaObs", 10.0);                                          // [GMRF model] The initial weight (Lambda) of each observation
	GMRF_lambdaObsLoss = declare_parameter<double>("GMRF_lambdaObsLoss", 0.0);                                   // [GMRF model] The loss of information (Lambda) of the observations with each iteration (see AppTick)

	colormap = declare_parameter<std::string>("colormap", "jet");
	max_pclpoints_cell = declare_parameter<int>("max_pclpoints_cell", 20);
	min_sensor_val = declare_parameter<double>("min_sensor_val", 0.0);
	max_sensor_val = declare_parameter<double>("max_sensor_val", 0.0);

	suggest_next_location_sensor_th = declare_parameter<double>("suggest_next_location_sensor_th", 0.1);

	//----------------------------------
	// Subscriptions
	//----------------------------------
	sub_sensor = create_subscription<olfaction_msgs::msg::Anemometer>(sensor_topic, 1, std::bind(&Cgmrf::sensorCallback, this, _1));
	ocupancyMap_sub = create_subscription<nav_msgs::msg::OccupancyGrid>("map", 1, std::bind(&Cgmrf::mapCallback, this, _1));
	//----------------------------------
	// Publishers
	//----------------------------------
	wind_array_pub = create_publisher<visualization_msgs::msg::MarkerArray>("wind_array_pub", 1);
	//----------------------------------
	// Services
	//----------------------------------
	// rclcpp::ServiceServer service = param_n.advertiseService("suggestNextObservationLocation", suggestNextObservationLocation);

	verbose = declare_parameter<bool>("verbose", false);

	module_init = false;
}

Cgmrf::~Cgmrf() {}

//--------------------------
// CALLBACK - OCCUPANCY MAP
//--------------------------
void Cgmrf::mapCallback(const nav_msgs::msg::OccupancyGrid::SharedPtr msg)
{
	if (module_init)
		return;

	RCLCPP_DEBUG(get_logger(), "[GMRF-node] %s - Map of the environment!", __FUNCTION__);
	occupancyMap = *msg;

	// Set GasMap dimensions as the OccupancyMap
	double map_min_x = msg->info.origin.position.x;
	double map_max_x = msg->info.origin.position.x + msg->info.width * msg->info.resolution;
	double map_min_y = msg->info.origin.position.y;
	double map_max_y = msg->info.origin.position.y + msg->info.height * msg->info.resolution;


	// Create GMRF-Map and init
	my_map = std::make_unique<CGMRF_map>(this, occupancyMap, cell_size, GMRF_lambdaPrior_reg, GMRF_lambdaPrior_mass_conservation, GMRF_lambdaPrior_obstacles, colormap, max_pclpoints_cell, verbose);
	RCLCPP_INFO(get_logger(), "[GMRF-node] GMRF GridMap initialized");

	module_init = true;
}

//----------------------------------
// CALLBACK - NEW WIND OBSERVATION
//----------------------------------
void Cgmrf::sensorCallback(const olfaction_msgs::msg::Anemometer::SharedPtr msg)
{
	static auto tf_buffer = std::make_unique<tf2_ros::Buffer>(get_clock());
	static auto tf_listener = std::make_shared<tf2_ros::TransformListener>(*tf_buffer);

	// 1. Get wind measurement
	double downwind_direction_map;
	mutex_anemometer.lock();
	try
	{
		reading_speed = msg->wind_speed;         // (m/s)
		reading_direction = msg->wind_direction; // (rad) This is the Upwind direction with respect the Anemometer ref system (standard measurement)

		// We need to transform this Upwind direction in the Anemometer ref system---- to ---- DownWind direction in the MAP ref system
		if (reading_speed != 0.0)
		{
			// Transform from anemometer ref_system to the map ref_system using TF
			geometry_msgs::msg::PoseStamped anemometer_upWind_pose, map_upWind_pose;
			try
			{
				anemometer_upWind_pose.header.frame_id = msg->header.frame_id.c_str();
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
				RCLCPP_ERROR(get_logger(), "[GMRF] - %s - Error: %s", __FUNCTION__, ex.what());
			}
		}
		else
		{
			downwind_direction_map = 0.0;
		}
	}
	catch (std::exception e)
	{
		RCLCPP_ERROR(get_logger(), "[GMRF] Exception at new Obs: %s ", e.what());
	}
	mutex_anemometer.unlock();
	// RCLCPP_INFO(get_logger(), "[GMRF-node] New wind observation! %.2f m/s  %.2f rad (DownWind in the map ref system)",msg->wind_speed, downwind_direction_map);

	// 2. Get pose of the sensor in the map reference system
	geometry_msgs::msg::TransformStamped transform;
	bool know_sensor_pose = true;
	try
	{
		// lookuptransform (target_frame, source_frame, result_tf)
		transform = tf_buffer->lookupTransform(frame_id.c_str(), msg->header.frame_id.c_str(), rclcpp::Time(0));
	}
	catch (tf2::TransformException ex)
	{
		RCLCPP_ERROR(get_logger(), "[GMRF] exception when reading observation: %s", ex.what());
		know_sensor_pose = false;
	}

	// 3. Add observation to the GMRF map
	if (module_init)
	{
		// Current sensor pose in the map
		float x_pos = transform.transform.translation.x;
		float y_pos = transform.transform.translation.y;

		mutex_anemometer.lock();
		// RCLCPP_INFO(get_logger(), "[GMRF] New obs: %.2f m/s, %.2f rad at (%.2f,%.2f)", reading_speed,reading_direction,x_pos,y_pos);
		my_map->insertObservation_GMRF(reading_speed, downwind_direction_map, x_pos, y_pos, GMRF_lambdaObs);
		mutex_anemometer.unlock();
	}
}

void Cgmrf::publishMaps()
{
	visualization_msgs::msg::MarkerArray wind_array;
	my_map->get_as_markerArray(wind_array, frame_id);
	wind_array_pub->publish(wind_array);
}

bool Cgmrf::get_wind_value_srv(WindEstimation::Request::SharedPtr req, WindEstimation::Response::SharedPtr res)
{
	// Since the wind fields are identical among different instances, return just the information from instance[0]
	for (int i = 0; i < req->x.size(); i++)
	{
		Eigen::Vector3d r = my_map->getEstimation(req->x[i], req->y[i]);
		res->u.push_back(r.x());
		res->v.push_back(r.y());
		res->stdev_angle.push_back(r.z());
	}
	return true;
}

//-----------------------------------------------------------------------------
//                                    MAIN
//----------------------------------------------------------------------------
int main(int argc, char** argv)
{
	rclcpp::init(argc, argv);
	auto my_gmrf_map = std::make_shared<Cgmrf>();

	auto service = my_gmrf_map->create_service<WindEstimation>("WindEstimation", std::bind(&Cgmrf::get_wind_value_srv, my_gmrf_map.get(), _1, _2));
	RCLCPP_INFO(my_gmrf_map->get_logger(), "[gmrf] LOOP....");
	rclcpp::Time last_publication_time = my_gmrf_map->now();
	rclcpp::Rate loop_rate(my_gmrf_map->exec_freq);

	while (rclcpp::ok())
	{
		rclcpp::spin_some(my_gmrf_map); // Callbacks & Services

		if (my_gmrf_map->module_init)
		{
			// Update and Publish maps
			my_gmrf_map->my_map->updateMapEstimation_GMRF(my_gmrf_map->GMRF_lambdaObsLoss);
			my_gmrf_map->publishMaps();

			// Info about actual rate
			if (my_gmrf_map->verbose)
				RCLCPP_INFO(my_gmrf_map->get_logger(), "[gmrf] Updating every %f seconds. Intended preiod was %f", (my_gmrf_map->now() - last_publication_time).seconds(), 1.0 / my_gmrf_map->exec_freq);
			last_publication_time = my_gmrf_map->now();
		}
		else
		{
			if (my_gmrf_map->verbose)
				RCLCPP_INFO(my_gmrf_map->get_logger(), "[gmrf] Waiting for initialization (Map of environment).");
		}
		loop_rate.sleep();
	}
}
