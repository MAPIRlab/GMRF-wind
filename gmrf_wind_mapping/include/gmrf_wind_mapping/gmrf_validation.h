/*********************************************************************
* gmrf_validation.h - GMRF Wind-Distribution Mapping ROS2 validation Node
* Description: Implements the Gaussian Markov random field mapping algorithm for the
*              estimation of the windflow from a set of sparse 2D wind measurements.
*              This node is used to validate the GMRF_wind_mapping module, comparing its 
*              results with ground-truth data from a CFD dataset (GADEN, VGR dataset).
* Authors: Javier Monroy, Pepe Ojeda 
*         University of Malaga (Spain)
*         http://mapir.isa.uma.es
*********************************************************************/

#include "rclcpp/rclcpp.hpp"
#include "std_msgs/msg/float32.hpp"
#include "nav_msgs/msg/odometry.hpp"
#include "nav_msgs/msg/occupancy_grid.hpp"
#include <sensor_msgs/msg/point_cloud2.hpp>
#include "geometry_msgs/msg/pose_with_covariance_stamped.hpp"
#include <visualization_msgs/msg/marker_array.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/math/constants/constants.hpp>
#include <angles/angles.h>


// Core GMRF map
#include "gmrf_wind_mapping/gmrf_map.h"


struct WindVectorXY
{
    double x;
    double y;
};

// Cvalgt class (ROS2 Node)
class Cvalgt : public rclcpp::Node
{
public:
    Cvalgt();
    ~Cvalgt();
    void update_parameters();
    void update();
    void publishMaps();
    void compute_performance_metrics();
    bool module_init;
    bool verbose;
    bool visualize_gmrf;
protected:
    inline void ReadMap();
    inline void ReadGroundTruthWindMap(const std::string& filename);
    void SimulateWindObservations();
    void initialize();

    // Publishers (visualization)
    rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr wind_array_pub;
    rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr wind_std_array_pub;
    rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr wind_gt_array_pub;
    rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr metric_pub;

    // GMRF Map
    std::unique_ptr<CGMRF_map> gmrf_map;        // The GMRF Map being estimated
    std::vector<WindVectorXY> gt_map;           // GT wind map from CFD data
    nav_msgs::msg::OccupancyGrid occupancy_map; // Occupancy GridMap of the environment

    // Node Parameters
    std::string mapFilePath;    // Path to the map image file (grayscale) to use as occupancy map
    std::string cfdFilePath;    // Path to the CSV file containing the CFD ground-truth wind data
    double cell_size;
    
    // GMRF parameters
    double GMRF_lambdaPrior_reg;               // Weight for regularization prior -> neighbour cells have similar wind vectors
    double GMRF_lambdaPrior_mass_conservation; // Weight for mass conservation law prior
    double GMRF_lambdaPrior_obstacles;         // Weight for wind close to obstacles prior -->cells close to obstacles has only tangencial wind
    double GMRF_lambdaObs;     // [GMRF model] The initial information (Lambda) of each observation (this information will decrease with time)
    double GMRF_lambdaObsLoss; // [GMRF model] The loss of information (Lambda) of the observations with each iteration (see AppTick)

    // Variables
    boost::mutex mutex_anemometer;
    boost::mutex mutex_position;
    std::string metrics_filename;
};