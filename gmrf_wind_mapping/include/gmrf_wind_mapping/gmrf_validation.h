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
    std::vector<double> compute_performance_metrics(const std::string& metric) const;
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
        if (optimization_metric == "NLPD" || optimization_metric == "ANLPD")
            gmrf_map->computeUncertainty_GMRF();
        
        // 5. Residuals for optimization (per-cell) or average metric (single value)
        std::vector<double> selected_residuals = this->compute_performance_metrics(optimization_metric);
        if (selected_residuals.empty()) {
            std::cerr << "Error: compute_performance_metrics returned empty residuals." << std::endl;
            return false; 
        }

        // 6. Fill the Ceres residual array
        for (size_t i = 0; i < selected_residuals.size(); ++i)
        {
            if (std::isnan(selected_residuals[i]) || std::isinf(selected_residuals[i])) 
            {
                std::cerr << "[CERES ABORT] MAP Solver yielded NaN at residual " << i 
                            << ". Alphas: [" << alpha[0] << ", " << alpha[1] << ", " 
                            << alpha[2] << ", " << alpha[3] << "]\n";
                return false; // Tell Ceres this step caused numerical instability
            }
            residual[i] = static_cast<T>(selected_residuals[i]);
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
    std::unique_ptr<CGMRF_map> gmrf_map;        // The GMRF Map being estimated
    std::vector<WindVectorXY> gt_map;           // GT wind map from CFD data
    nav_msgs::msg::OccupancyGrid occupancy_map; // Occupancy GridMap of the environment
    std::vector<int> free_cell_indices;         // List of free cell indices for faster iteration (used in estimation and metrics computation)

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
    std::string optimization_metric;           // Metric to optimize when tuning the GMRF hyperparameters (lambda weights and observation variances). Options: "AE", "ME", "RSE", "NSP", "NLPD", "ME&AE"

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
    
    // -------------------------------------------------------------------
    // SIGNATURE 1: For DynamicNumericDiffCostFunction (Array of pointers)
    // -------------------------------------------------------------------
    template <typename T>
    bool operator()(T const* const* parameters, T* residual) const
    {
        // Unpack the parameter block array
        const T* const alpha = parameters[0];
        const T* const log_k = parameters[1];

        // Delegate the work using the member variable 'optimization_metric'
        return gmrf_validation_node->evaluate_cost(optimization_metric, alpha, log_k, residual);
    }

    // -------------------------------------------------------------------
    // SIGNATURE 2: For standard NumericDiffCostFunction (Variadic args)
    // -------------------------------------------------------------------
    template <typename T>
    bool operator()(const T* const alpha, const T* const log_k, T* residual) const
    {
        // Delegate the work using the member variable 'optimization_metric'
        return gmrf_validation_node->evaluate_cost(optimization_metric, alpha, log_k, residual);
    }
};