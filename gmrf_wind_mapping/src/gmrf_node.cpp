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

// -------------------
// Cgmrf Constructor
//--------------------
Cgmrf::Cgmrf()
{
    printf("\n=================================================================");
    printf("\n=	             GMRF Wind-Distribution Mapping Node              =");
    printf("\n=================================================================\n");

    //------------------
    // Load Parameters
    //------------------
    ros::NodeHandle param_n("~");
    param_n.param<std::string>("frame_id", frame_id, "map");
    param_n.param<std::string>("sensor_topic", sensor_topic, "/anemometer");
    param_n.param<double>("exec_freq", exec_freq, 2.0);
    param_n.param<double>("cell_size", cell_size, 0.5);

    param_n.param<double>("GMRF_lambdaPrior_reg", GMRF_lambdaPrior_reg, 1);   // Weight for regularization prior -> neighbour cells have similar wind vectors
    param_n.param<double>("GMRF_lambdaPrior_mass_conservation", GMRF_lambdaPrior_mass_conservation, 10000);   // Weight for mass conservation law prior
    param_n.param<double>("GMRF_lambdaPrior_obstacles", GMRF_lambdaPrior_obstacles, 10);   //Weight for wind close to obstacles prior -->cells close to obstacles has only tangencial wind
    param_n.param<double>("GMRF_lambdaObs", GMRF_lambdaObs, 10.0);       // [GMRF model] The initial weight (Lambda) of each observation
    param_n.param<double>("GMRF_lambdaObsLoss", GMRF_lambdaObsLoss, 0.0);// [GMRF model] The loss of information (Lambda) of the observations with each iteration (see AppTick)


    param_n.param<std::string>("colormap", colormap, "jet");
    param_n.param<int>("max_pclpoints_cell", max_pclpoints_cell, 20);
    param_n.param<double>("min_sensor_val", min_sensor_val, 0.0);
    param_n.param<double>("max_sensor_val", max_sensor_val, 0.0);

    param_n.param<double>("suggest_next_location_sensor_th", suggest_next_location_sensor_th, 0.1);



    //----------------------------------
    // Subscriptions
    //----------------------------------
    sub_sensor = param_n.subscribe(sensor_topic, 1, &Cgmrf::sensorCallback, this);
    ocupancyMap_sub = param_n.subscribe("map", 1, &Cgmrf::mapCallback, this);
    //----------------------------------
    // Publishers
    //----------------------------------
    wind_array_pub = param_n.advertise<visualization_msgs::MarkerArray>("wind_array_pub", 0);
    //----------------------------------
    // Services
    //----------------------------------
    //ros::ServiceServer service = param_n.advertiseService("suggestNextObservationLocation", suggestNextObservationLocation);
    

    module_init = false;
}


Cgmrf::~Cgmrf(){}


//--------------------------
// CALLBACK - OCCUPANCY MAP
//--------------------------
void Cgmrf::mapCallback(const nav_msgs::OccupancyGrid::ConstPtr& msg)
{
    if(module_init)
        return;

    ROS_DEBUG("[GMRF-node] %s - Map of the environment!", __FUNCTION__);
    occupancyMap = *msg;

    //Set GasMap dimensions as the OccupancyMap
    double map_min_x = msg->info.origin.position.x;
    double map_max_x =  msg->info.origin.position.x + msg->info.width*msg->info.resolution;
    double map_min_y =  msg->info.origin.position.y;
    double map_max_y =  msg->info.origin.position.y + msg->info.height*msg->info.resolution;

    //Create GMRF-Map and init
    my_map = new CGMRF_map(occupancyMap, cell_size, GMRF_lambdaPrior_reg, GMRF_lambdaPrior_mass_conservation, GMRF_lambdaPrior_obstacles, colormap, max_pclpoints_cell);
    ROS_INFO("[GMRF-node] GMRF GridMap initialized");

    module_init = true;
}


//----------------------------------
// CALLBACK - NEW WIND OBSERVATION
//----------------------------------
void Cgmrf::sensorCallback(const olfaction_msgs::anemometerConstPtr msg)
{
    //ROS_INFO("[GMRF-node] New wind observation! %.2f m/s  %.2f rad (Upwind in the Anemometer ref system)",msg->wind_speed, msg->wind_direction);

    //1. Get wind measurement
    double downwind_direction_map;
    mutex_anemometer.lock();
    try
    {
        reading_speed = msg->wind_speed;            // (m/s)
        reading_direction = msg->wind_direction;    // (rad) This is the Upwind direction with respect the Anemometer ref system (standard measurement)


        //We need to transform this Upwind direction in the Anemometer ref system---- to ---- DownWind direction in the MAP ref system
        if (reading_speed != 0.0)
        {
            //Transform from anemometer ref_system to the map ref_system using TF
            geometry_msgs::PoseStamped anemometer_upWind_pose, map_upWind_pose;
            try
            {
                anemometer_upWind_pose.header.frame_id = msg->header.frame_id.c_str();
                anemometer_upWind_pose.pose.position.x = 0.0;
                anemometer_upWind_pose.pose.position.y = 0.0;
                anemometer_upWind_pose.pose.position.z = 0.0;
                anemometer_upWind_pose.pose.orientation = tf::createQuaternionMsgFromYaw(msg->wind_direction);

                //lookuptransform (target_frame, target_time, pose_in, fixed_frame, pose_out)
                tf_listener.transformPose(frame_id.c_str(), anemometer_upWind_pose, map_upWind_pose);

                downwind_direction_map = angles::normalize_angle( tf::getYaw(map_upWind_pose.pose.orientation) + 3.14159 );
            }
            catch(tf::TransformException &ex)
            {
                ROS_ERROR("[GMRF] - %s - Error: %s", __FUNCTION__, ex.what());
            }
        }
        else
        {
            downwind_direction_map = 0.0;
        }
    }
    catch(std::exception e){
        ROS_ERROR("[GMRF] Exception at new Obs: %s ", e.what() );
    }
    mutex_anemometer.unlock();
    //ROS_INFO("[GMRF-node] New wind observation! %.2f m/s  %.2f rad (DownWind in the map ref system)",msg->wind_speed, downwind_direction_map);


    //2. Get pose of the sensor in the map reference system
    tf::StampedTransform transform;
    bool know_sensor_pose = true;
    try
    {
        //lookuptransform (target_frame, source_frame, result_tf)
        tf_listener.lookupTransform(frame_id.c_str(),msg->header.frame_id.c_str(), ros::Time(0), transform);
    }
    catch (tf::TransformException ex)
    {
        ROS_ERROR("[GMRF] exception when reading observation: %s",ex.what());
        know_sensor_pose = false;
        ros::Duration(1.0).sleep();
    }


    //3. Add observation to the GMRF map
    if (module_init)
    {
        //Current sensor pose in the map
        float x_pos = transform.getOrigin().x();
        float y_pos = transform.getOrigin().y();

        mutex_anemometer.lock();
        //ROS_INFO("[GMRF] New obs: %.2f m/s, %.2f rad at (%.2f,%.2f)", reading_speed,reading_direction,x_pos,y_pos);
        my_map->insertObservation_GMRF(reading_speed, downwind_direction_map, x_pos, y_pos, GMRF_lambdaObs);
        mutex_anemometer.unlock();
    }
}




void Cgmrf::publishMaps()
{
    visualization_msgs::MarkerArray wind_array;
    my_map->get_as_markerArray(wind_array, frame_id);
    wind_array_pub.publish(wind_array);
}


bool Cgmrf::get_wind_value_srv(gmrf_wind_mapping::WindEstimation::Request  &req, gmrf_wind_mapping::WindEstimation::Response &res)
{
    //Since the wind fields are identical among different instances, return just the information from instance[0]
    for(int i=0;i<req.x.size();i++){
        Eigen::Vector3d r = my_map->getEstimation(req.x[i], req.y[i]);
        res.u.push_back(r.x());
        res.v.push_back(r.y());
        res.stdevAngle.push_back(r.z());
    }
    return true;
}

//-----------------------------------------------------------------------------
//                                    MAIN
//----------------------------------------------------------------------------
int main(int argc, char **argv)
{
	ros::init(argc, argv, "gmrf_node");
    ros::NodeHandle n;
	Cgmrf my_gmrf_map;
    ros::ServiceServer service = n.advertiseService("WindEstimation", &Cgmrf::get_wind_value_srv, &my_gmrf_map);
	ROS_INFO("[gmrf] LOOP....");
	ros::Time last_publication_time = ros::Time::now();
	ros::Rate loop_rate(my_gmrf_map.exec_freq);
	
	while (ros::ok())
	{
		ros::spinOnce();                    //Callbacks & Services

		if (my_gmrf_map.module_init)
		{
			// Update and Publish maps
			my_gmrf_map.my_map->updateMapEstimation_GMRF(my_gmrf_map.GMRF_lambdaObsLoss);
            my_gmrf_map.publishMaps();
			
			// Info about actual rate
			ROS_INFO("[gmrf] Updating every %f seconds. Intended preiod was %f", (ros::Time::now()-last_publication_time).toSec(), 1.0/my_gmrf_map.exec_freq);
			last_publication_time = ros::Time::now();
		}
        else
        {
			ROS_INFO("[gmrf] Waiting for initialization (Map of environment).");
		}
        	loop_rate.sleep();
	}
}







