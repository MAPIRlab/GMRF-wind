/*********************************************************************
*
* Software License Agreement (BSD License)
*
*  Copyright (c)  2015, Örebro University, Sweden
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.

* Author: Javier G. Monroy
* KernelDM+V Implementation:  Victor Hernandez
*********************************************************************/

//-----------------------------------------------------
// 2D Wind estimation with GMRF
//-----------------------------------------------------
#include "ros/ros.h"
#include "std_msgs/Float32.h"
#include "nav_msgs/Odometry.h"
#include "nav_msgs/OccupancyGrid.h"
#include <sensor_msgs/PointCloud2.h>
#include "geometry_msgs/PoseWithCovarianceStamped.h"
#include <visualization_msgs/MarkerArray.h>
#include <boost/thread/mutex.hpp>
#include <boost/math/constants/constants.hpp>
#include <tf/transform_listener.h>
#include <angles/angles.h>
#include "olfaction_msgs/anemometer.h"

#include "gmrf_map.h"


//Services
#include "olfaction_msgs/suggestNextObservationLocation.h"
#include "gmrf_wind_mapping/WindEstimation.h"



class Cgmrf
{
public:
    Cgmrf();
    ~Cgmrf();
    void publishMaps();
    bool get_wind_value_srv(gmrf_wind_mapping::WindEstimation::Request  &req, gmrf_wind_mapping::WindEstimation::Response &res);

    //GMRF variables    
    CGMRF_map                               *my_map;			// The Online Gas Distribution Map being generated
    nav_msgs::OccupancyGrid                 occupancyMap;       //Occupancy GridMap of the environment

    //Node Params
    std::string                             sensor_topic;
    std::string                             frame_id;           //frame where to plot the map, usually (/map)
    double                                  cell_size;
    double                                  exec_freq;
    std::string                             colormap;
    int                                     max_pclpoints_cell;
    double                                  max_sensor_val;
    double                                  min_sensor_val;
    double                                  suggest_next_location_sensor_th;


    double                                  GMRF_lambdaPrior_reg;                   // Weight for regularization prior -> neighbour cells have similar wind vectors
    double                                  GMRF_lambdaPrior_mass_conservation;     // Weight for mass conservation law prior
    double                                  GMRF_lambdaPrior_obstacles;             // Weight for wind close to obstacles prior -->cells close to obstacles has only tangencial wind
    double                                  GMRF_lambdaObs;         // [GMRF model] The initial information (Lambda) of each observation (this information will decrease with time)
    double                                  GMRF_lambdaObsLoss;     // [GMRF model] The loss of information (Lambda) of the observations with each iteration (see AppTick)



    //Variables
    bool module_init;
    boost::mutex                            mutex_anemometer;
    boost::mutex                            mutex_position;    
    double                                  reading_speed;      // m/s
    double                                  reading_direction;  // rad
    bool                                    new_data_position;
    float                                   curr_x;
    float                                   curr_y;



protected:
    ros::NodeHandle                         n;
    tf::TransformListener                   tf_listener; //Do not put inside the callback    

    //Subscriptions & Publishers
    ros::Subscriber sub_sensor, ocupancyMap_sub;
    ros::Publisher wind_array_pub;

    //Callbacks
    void sensorCallback(const olfaction_msgs::anemometerConstPtr msg);
    void mapCallback(const nav_msgs::OccupancyGrid::ConstPtr& msg);
};






//-------------------------------------------------------
//	Variables
//-------------------------------------------------------


/*
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
*/
