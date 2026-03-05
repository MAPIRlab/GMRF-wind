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
#include "ceres/ceres.h"


// Core GMRF map
#include "gmrf_wind_core/gmrf_map.h"


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

    void update();
    void publishMaps();
    std::array<double, 4> compute_performance_metrics() const;
    void SimulateWindObservations(size_t N_obs);
    void SimulateFixedWindObservations();
    void update_lambdas(double lambda_adv, double lambda_mass, double lambda_diff, double lambda_vort, double lambda_obst);
    void read_lambdas(double &lambda_adv, double &lambda_mass, double &lambda_diff, double &lambda_vort, double &lambda_obst);
    void saveGMRFEstimationToCSV(const std::string& file_name);
    bool module_init;
    bool verbose;
    bool visualize_gmrf;

    // CERES: The cost function for NLPD optimization
    template <typename T>
    bool evaluate_cost(const T* const params, T* residual) const
    {
        // 1. Update GMRF lambdas (Advection, Mass conservation, Diffusion, Vorticity, Obstacles) with the current optimization parameters
        gmrf_map->update_lambdas( static_cast<double>(params[0]),
                                  static_cast<double>(params[1]),
                                  static_cast<double>(params[2]),
                                  static_cast<double>(params[3]), 
                                  static_cast<double>(params[4])
                                );
        
        // 2. Run MAP estimation and uncertainty computation
        gmrf_map->MAP_estimation_GMRF(num_iterations_MAP);
        gmrf_map->computeUncertainty_GMRF();

        // 3. Compute performance metrics (AAE, RMSE, ANSP, NLPD)
        std::array<double, 4> metrics = this->compute_performance_metrics();
        double AAE = metrics[0];  // AAE is the first element
        double RMSE = metrics[1]; // RMSE is the second element
        double ANSP = (-metrics[2]+1); // ANSP is the third element (normalized dot product, we want to maximize it. adding offset to allow negative values)
        double ANLPD = metrics[3] + 20; // NLPD is the fourth element (adding offset to allow negative values)

        // Metric to minimize
        residual[0] = static_cast<T>(ANLPD);
        return true;
    }


    inline void ReadMap();
    inline void ReadGroundTruthWindMap(const std::string& filename);
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
    
    // GMRF precision parameters
    double GMRF_lambdaPrior_advection;         // Weight for advection prior (wind flow follows the flow lines)
    double GMRF_lambdaPrior_mass_conservation; // Weight for mass conservation (divergence free)
    double GMRF_lambdaPrior_diffusion;         // Weight for diffusion prior (smoothness)
    double GMRF_lambdaPrior_vorticity;         // Weight for vorticity prior (curl free) This one is dynamic, so this is the max value.
    double GMRF_lambdaPrior_obstacles;         // Weight for wind close to obstacles prior -->cells close to obstacles has only tangencial wind    
    double observation_var_wind_speed;         // Variance of the wind speed measurement (m/s)^2
    double observation_var_wind_direction;     // Variance of the wind direction measurement (rad)^2
    int num_iterations_MAP;                 // Maximum number of iterations for the MAP estimation optimization
    int experiment_number;                     // Number of the experiment to run
    
    // Variables
    boost::mutex mutex_anemometer;
    boost::mutex mutex_position;
    std::string metrics_filename;
};

// CERES Cost Functor for NLPD optimization
struct CostFunctor
{
    Cvalgt* gmrf_validation_node;
    CostFunctor(Cvalgt* instance) : gmrf_validation_node(instance) {}
    
    template <typename T>
    bool operator()(const T* const params, T* residual) const {
        // Delegate the work to the class member function
        return gmrf_validation_node->evaluate_cost(params, residual);
    }
};