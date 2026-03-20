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
    void clearEstimation();
    void update_lambdas(double lambda_adv, double lambda_mass, double lambda_diff, double lambda_obst);
    void read_lambdas(double &lambda_adv, double &lambda_mass, double &lambda_diff, double &lambda_obst);
    void saveGMRFEstimationToCSV(const std::string& file_name);
    bool module_init;
    bool verbose;
    bool visualize_gmrf;

    // CERES: The cost function for optimization
    //------------------------------------------
    template <typename T>
    bool evaluate_cost(const std::string& optimization_metric, const T* const alpha, const T* const log_k, T* residual) const
    {
        // 1. Softmax transform
        double exp_alpha[4];
        double sum_exp = 0.0;
        for(int i = 0; i < 4; ++i) {
            exp_alpha[i] = std::exp(alpha[i]);
            sum_exp += exp_alpha[i];
        }

        // 2. Scale factors (k) 
        double k = std::exp(log_k[0]);  // Ensure k is positive
        double lambda[4];               // GMRF precision parameters to be updated in the GMRF map
        for(int i = 0; i < 4; ++i) {
            lambda[i] = k * (exp_alpha[i] / sum_exp);
        }

        // 3. Update GMRF lambdas (Advection, Mass conservation, Diffusion, Obstacles)
        gmrf_map->update_lambdas( lambda[0], lambda[1], lambda[2], lambda[3] );

        // 4. Run MAP estimation and uncertainty computation (when needed)
        gmrf_map->MAP_estimation_GMRF(num_iterations_MAP);
        if (optimization_metric == "anlpd")
            gmrf_map->computeUncertainty_GMRF();

        // 5. Compute performance metrics (AAE, RMSE, ANSP, ANLPD)
        std::array<double, 4> metrics = this->compute_performance_metrics();
        double res_AAE = metrics[0];  // AAE (Average Angular Error) [0,1] --> 0 is the best.
        double res_RMSE = metrics[1]; // RMSE (Root Mean Square Error) [0,inf) --> 0 is the best.
        double res_ANSP = (1-metrics[2]); // ANSP [-1,1] (Average Normalized Scalar Product) --> [0,2] where 0 is the best.
        double res_ANLPD = metrics[3] + 20; // ANLPD (Average Negative Log Predictive Density) --> adding offset to allow negative values

        
        // Metric to minimize (residual for optimization)
        if (optimization_metric == "aae")
            residual[0] = static_cast<T>(res_AAE);
        else if (optimization_metric == "rmse")
            residual[0] = static_cast<T>(res_RMSE);
        else if (optimization_metric == "ansp")
        {
            // check
            if (std::isnan(res_ANSP) || std::isinf(res_ANSP)) {
                std::cerr << "Numerical instability detected in ANSP: " << res_ANSP << std::endl;
                return false; // Indica a Ceres que reduzca el tamaño del paso
            }
            residual[0] = static_cast<T>(res_ANSP);
            //std::cerr << "ANSP: " << metrics[2] << std::endl;
        }
        else if (optimization_metric == "anlpd")
        {
            // check
            if (std::isnan(res_ANLPD) || std::isinf(res_ANLPD)) {
                std::cerr << "Numerical instability detected in ANPD: " << res_ANLPD << std::endl;
                return false; // Indica a Ceres que reduzca el tamaño del paso
            }
            residual[0] = static_cast<T>(res_ANLPD);
            //std::cerr << "ANLPD: " << metrics[3] << std::endl;
        }

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
    std::unique_ptr<gmrfw::CGMRF_map> gmrf_map;        // The GMRF Map being estimated
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
    double GMRF_lambdaPrior_obstacles;         // Weight for wind close to obstacles prior -->cells close to obstacles has only tangencial wind    
    double observation_var_wind_speed;         // Variance of the wind speed measurement (m/s)^2
    double observation_var_wind_direction;     // Variance of the wind direction measurement (rad)^2
    int num_iterations_MAP;                    // Maximum number of iterations for the MAP estimation optimization
    int experiment_number;                     // Number of the experiment to run
    
    // Variables
    boost::mutex mutex_anemometer;
    boost::mutex mutex_position;
    std::string metrics_filename;
};

// CERES Cost Functor for optimization
struct CostFunctor
{
    Cvalgt* gmrf_validation_node;
    std::string optimization_metric;
    CostFunctor(Cvalgt* instance, const std::string metric) : gmrf_validation_node(instance), optimization_metric(metric) {}
    
    template <typename T>
    // bool operator()(const T* const params, T* residual) const
    bool operator()(const T* const alpha, const T* const log_k, T* residual) const
    {
        // Delegate the work to the class member function
        return gmrf_validation_node->evaluate_cost(optimization_metric, alpha, log_k, residual);
    }
};