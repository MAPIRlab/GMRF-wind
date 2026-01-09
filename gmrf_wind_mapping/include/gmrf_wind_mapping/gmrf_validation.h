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

    void update();
    void publishMaps();
    std::array<double, 3> compute_performance_metrics() const;
    void SimulateWindObservations(size_t N_obs);
    void update_lambdas(double lambda_reg, double lambda_flux, double lambda_obstacles, double lambda_obs);
    void read_lambdas(double &lambda_reg, double &lambda_flux, double &lambda_obstacles, double &lambda_obs);
    bool module_init;
    bool verbose;
    bool visualize_gmrf;

    // CERES: The cost function for NLPD optimization
    template <typename T>
    bool evaluate_cost(const T* const params, T* residual) const
    {
        // 1. Update GMRF lambdas (Reg, Flux, Obstacles, Observations)
        gmrf_map->update_lambdas( static_cast<double>(params[0]),
                                  static_cast<double>(params[1]),
                                  static_cast<double>(params[2]),
                                  static_cast<double>(params[3]) );
        
        // 2. Run MAP estimation and uncertainty computation
        gmrf_map->MAP_estimation_GMRF();
        gmrf_map->computeUncertainty_GMRF();

        // 3. Compute performance metrics (AAE, RMSE, NLPD)
        std::array<double, 3> metrics = this->compute_performance_metrics();
        double ANLPD = metrics[2]; // NLPD is the third element

        // 4. Return Average NLPD as residual (we optimize for NLPD minimization)
        residual[0] = static_cast<T>(ANLPD);
        return true;
    }

protected:
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
    
    // GMRF parameters
    double GMRF_lambdaPrior_reg;               // Weight for regularization prior -> neighbour cells have similar wind vectors
    double GMRF_lambdaPrior_flux_conservation; // Weight for flux conservation law prior
    double GMRF_lambdaPrior_obstacles;         // Weight for wind close to obstacles prior -->cells close to obstacles has only tangencial wind
    double GMRF_lambdaObs;     // [GMRF model] The initial information (Lambda) of each observation (this information will decrease with time)
    
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