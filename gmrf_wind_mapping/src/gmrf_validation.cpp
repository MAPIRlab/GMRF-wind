//========================================================================================
//	GMRF_validation - GMRF Wind-Distribution Mapping ROS2 validation Node
//  Description: Implements the Gaussian Markov random field mapping algorithm for the
//              estimation of the windflow from a set of sparse 2D wind measurements.
//              This node is used to validate the GMRF_wind_mapping module, comparing its 
//              results with ground-truth data from a CFD dataset (GADEN, VGR dataset).
//
//
//
//----------------------------------------------------------------------------------------
//	Revision log:
//	version: 2.0	20/11/2025
//========================================================================================

#include "gmrf_wind_mapping/gmrf_validation.h"
#include "gmrf_wind_mapping/utils.h"

#include <chrono>
#include <tf2/time.h>
#include <yaml-cpp/yaml.h>
#include <random>
#include <numeric>   // Required for std::inner_product
#include <cmath>     // Required for mathematical operations
#include <filesystem>

using namespace std::placeholders;
using namespace gmrfw;

Cvalgt::Cvalgt()
    : Node("GMRF_validation")
{
    // Load Parameters
    //------------------
    mapFilePath = declare_parameter<std::string>("map_yaml_file", "");
    cfdFilePath = declare_parameter<std::string>("cfd_csv_file", "");
    cell_size = declare_parameter<double>("cell_size", 0.5);
    verbose = declare_parameter<bool>("verbose", false);
    visualize_gmrf = declare_parameter<bool>("visualize_gmrf", true);

    // Observation variances (Gaussian likelihood)
    observation_var_wind_speed = declare_parameter<double>("observation_var_wind_speed", 0.0001);           // Variance of the wind speed measurement (m/s)^2
    observation_var_wind_direction = declare_parameter<double>("observation_var_wind_direction", 0.0001);       // Variance of the wind direction measurement (rad)^2

    // Lambda/weights for the different priors
    GMRF_lambdaPrior_advection = declare_parameter<double>("GMRF_lambdaPrior_advection", 1.0);
    GMRF_lambdaPrior_mass_conservation = declare_parameter<double>("GMRF_lambdaPrior_mass_conservation", 1.0);
    GMRF_lambdaPrior_diffusion = declare_parameter<double>("GMRF_lambdaPrior_diffusion", 1.0);
    GMRF_lambdaPrior_obstacles = declare_parameter<double>("GMRF_lambdaPrior_obstacles", 1.0);
    num_iterations_MAP = declare_parameter<int>("num_iterations_MAP", 100);

    // Experiment
    experiment_number = declare_parameter<int>("experiment_number", 0);
    optimization_metric = declare_parameter<std::string>("optimization_metric", "RSE");

    // Publishers
    //----------------------------------
    wind_array_pub = create_publisher<visualization_msgs::msg::MarkerArray>("wind_array_pub", 1);
    wind_std_array_pub = create_publisher<visualization_msgs::msg::Marker>("wind_std_array_pub", 1);
    metric_pub = create_publisher<visualization_msgs::msg::Marker>("metric_pub", 1);

    rclcpp::QoS qos_profile(1); // Keep the history depth low (1 is often sufficient for static data)
    qos_profile.durability(RMW_QOS_POLICY_DURABILITY_TRANSIENT_LOCAL);
    wind_gt_array_pub = create_publisher<visualization_msgs::msg::MarkerArray>("wind_gt_array_pub", qos_profile);
    
    
    // 1. Read OccupancyGridMap from file
    ReadMap();

    // 2. Initialize GMRF estimation
    initialize();

    // 3. Read Ground-Truth Wind Map
    ReadGroundTruthWindMap(cfdFilePath);
}


Cvalgt::~Cvalgt()
{
}


// Read Occupancy GridMap from a YAML file
inline void Cvalgt::ReadMap()
{
    // we can choose to read a map file directly from disk, if we don't want to use the one published by map_server
    // this is often useful when we need to alter the occupancy map to include inlets/outlets, which should be non-occupied
    // for GMRF but which may not be navigable (think of windows). 
    try
    {
        if (verbose)
            RCLCPP_INFO(get_logger(), "[Cvalgt] Reading Occupancy GridMap from file: %s", mapFilePath.c_str());

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
        occupancy_map.data.resize(width * height);

        double free_thresh = 1 - yaml["free_thresh"].as<double>(); // Unused, this is for trinary maps. We only consider cells free or occupied, not "unknown"
        double occupied_thresh = 1 - yaml["occupied_thresh"].as<double>();
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double value = mapImage.at<uint8_t>(y, x) / 255.0;
                // Calculate the 1D index (Handles the ROS bottom-left origin flip)
                int cell_idx = width * (height - y - 1) + x;
                int8_t& occupancy = occupancy_map.data.at(cell_idx);
                if (value > free_thresh) {
                    occupancy = 0;  // Free cell
                }
                else if (value < occupied_thresh) {
                    occupancy = 100;  // Occupied cell
                }
                else {
                    occupancy = -1; // Unknown
                }
            }
        }

        // Metadata
        //-------------------
        occupancy_map.header.frame_id = "map";
        occupancy_map.info.map_load_time = now();
        occupancy_map.info.resolution = yaml["resolution"].as<double>();
        auto origin_array = yaml["origin"].as<std::array<double, 3>>();
        occupancy_map.info.origin.position.x = origin_array[0];
        occupancy_map.info.origin.position.y = origin_array[1];
        occupancy_map.info.origin.orientation = Utils::createQuaternionMsgFromYaw(origin_array[2]);

        occupancy_map.info.height = height;
        occupancy_map.info.width = width;

        // Publish the actual map we are using, in case it is different from the one published by map_server
        static rclcpp::Publisher<nav_msgs::msg::OccupancyGrid>::SharedPtr map_republisher =
        create_publisher<nav_msgs::msg::OccupancyGrid>("gmrf_occupancy", rclcpp::QoS(1).transient_local());
        map_republisher->publish(occupancy_map);
    }   
    catch (const std::exception& e)
    {
        RCLCPP_ERROR(get_logger(), "[Cvalgt] Exception caught: '%s'", e.what());
    }
}



void Cvalgt::initialize()
{
    // Parse OccupancyGridmap to internal GMRF format
    TOccupancyMap occMap;
    occMap.data = occupancy_map.data;
    occMap.resolution = occupancy_map.info.resolution;              // cell size (m)
    occMap.width = occupancy_map.info.width;                        // number of cells in X direction
    occMap.height = occupancy_map.info.height;                      // number of cells in Y direction
    occMap.origin_x = occupancy_map.info.origin.position.x;         // world coordinates of the origin of the map
    occMap.origin_y = occupancy_map.info.origin.position.y;         // world coordinates of the origin of the map
    
    // Create the GMRF-Map and initialize its Prior Factors with the parameters from the ROS2 params server
    // Set only factors with a Lambda_prior > 0, (avoids adding unnecessary factors to the graph)
    CGMRF_map::Parameters params{
        cell_size,
        GMRF_lambdaPrior_advection,
        GMRF_lambdaPrior_mass_conservation,
        GMRF_lambdaPrior_diffusion,
        GMRF_lambdaPrior_obstacles,
    };
    gmrf_map = std::make_unique<CGMRF_map>(occMap, 
                                        params,
                                        verbose,
                                        true // estimateTiming
                                        );
    gmrf_map->update_lambdas(GMRF_lambdaPrior_advection, GMRF_lambdaPrior_mass_conservation, GMRF_lambdaPrior_diffusion, GMRF_lambdaPrior_obstacles);
    RCLCPP_INFO(get_logger(), "[GMRF-validation] GMRF Initialized");
    module_init = true;

    // =======================================
    // PRE-COMPUTE VALID CELLS FOR OPTIMIZATION (only free cells that are not boundaries)
    // =======================================
    free_cell_indices.clear();
    Eigen::Vector2i dimensions = gmrf_map->map_size();
    size_t N = dimensions.x() * dimensions.y();

    for (size_t i = 0; i < N; ++i) {
        if (gmrf_map->is_cell_free(i) && !gmrf_map->is_cell_boundary(i)) {
            free_cell_indices.push_back(i);
        }
    }
    RCLCPP_INFO(get_logger(), "[gmrf-validation] Pre-calculated %zu valid free cells for optimization.", free_cell_indices.size());
}



inline void Cvalgt::ReadGroundTruthWindMap(const std::string& filename)
{
    // From Gaden_preprocessing::openFoam_to_gaden
    try
    {
        if (verbose)
            RCLCPP_INFO(get_logger(), "[Cvalgt] Reading Ground-Truth Wind Map from CFD file: %s", filename.c_str());

        // let's parse the file
        std::ifstream infile(filename.c_str());
        std::string line;
        struct ParsedLine
        {
            float point[3];
            float windVector[3];
        };
        ParsedLine parsedLine;

        // Depending on the verion of Paraview used to export the file, lines might be (Point, vector) OR (vector, Point)
        // so we need to check the header before we know where to put what
        float* firstPartOfLine;
        float* secondPartOfLine;
        {
            std::getline(infile, line);
            size_t pos = line.find(",");
            std::string firstElement = line.substr(0, pos);

            if (firstElement.find("Points") != std::string::npos)
            {
                firstPartOfLine = parsedLine.point;
                secondPartOfLine = parsedLine.windVector;
            }
            else
            {
                firstPartOfLine = parsedLine.windVector;
                secondPartOfLine = parsedLine.point;
            }
        }

        // Keep a 2D grid of the environment (same size as the GMRF estimation)
        Eigen::Vector2i dimensions = gmrf_map->map_size();                              // (x,y) cells
        std::array<double,4> dimensions_meters = gmrf_map->map_dimensions_meters();     // (x_min,x_max,y_min,y_max) in meters  
        
        gt_map.resize(dimensions.x()*dimensions.y());

        int x_idx = 0;
        int y_idx = 0;
        int z_idx = 0;
        while (std::getline(infile, line))
        {
            if (line.length() != 0)
            {
                for (int i = 0; i < 3; i++)
                {
                    size_t pos = line.find(",");
                    firstPartOfLine[i] = atof(line.substr(0, pos).c_str());
                    line.erase(0, pos + 1);
                }

                for (int i = 0; i < 3; i++)
                {
                    size_t pos = line.find(",");
                    secondPartOfLine[i] = atof(line.substr(0, pos).c_str());
                    line.erase(0, pos + 1);
                }

                // assign each of the points we have information about to the nearest cell
                x_idx = (parsedLine.point[0] - dimensions_meters[0]) / cell_size;
                y_idx = (parsedLine.point[1] - dimensions_meters[2]) / cell_size;
                z_idx = (parsedLine.point[2]) / cell_size;

                // Keep only 2D data (XY) --> z_idx == 0
                if (z_idx == 0)
                {
                    size_t idx = x_idx + y_idx * dimensions.x();
                    // Ensure cell is free
                    if (gmrf_map->is_cell_free(idx))
                    {                    
                        gt_map[idx].x = parsedLine.windVector[0];
                        gt_map[idx].y = parsedLine.windVector[1];
                    }
                }
            }
        }

        infile.close();
        if (verbose)
            RCLCPP_INFO(get_logger(), "[Cvalgt] Ground-Truth Wind Map loaded from CFD data. Total cells with wind data: %zu", gt_map.size());

        // Publish GT wind map as markers (RVIZ2) - Only Once
        if (visualize_gmrf)
        {
            if (verbose)
                RCLCPP_INFO(get_logger(), "[Cvalgt] Publishing Ground-Truth Wind Map as Markers in RVIZ2");

            using visualization_msgs::msg::Marker;
            using visualization_msgs::msg::MarkerArray;
            MarkerArray wind_gt_array;

            // Get GT map dimensions
            size_t N = gt_map.size();
            size_t width = static_cast<size_t>(dimensions.x());
            size_t height = static_cast<size_t>(dimensions.y());
            if (N == 0)
                return;

            // Get max wind vector in the GT map (to normalize the plot)
            double max_module = 0.0;
            for (size_t i = 0; i < N; ++i)
            {
                WindVectorXY w = gt_map[i];
                double module = sqrt(pow(w.x,2) + pow(w.y,2));
                if (module > max_module)
                    max_module = module;
            }

            // One ARROW marker per cell for wind_vector
            wind_gt_array.markers.reserve(N);
            for (size_t i = 0; i < N; ++i)
            {
                // Get wind estimation at cell i
                WindVectorXY w = gt_map[i];
                double module = sqrt(pow(w.x,2) + pow(w.y,2));
                double direction = atan2(w.y, w.x);
                
                // Cell center coordinates
                double cx = dimensions_meters[0] + (i % width + 0.5) * cell_size;
                double cy = dimensions_meters[2] + (i / width + 0.5) * cell_size;

                // Create Arrow marker
                Marker m;
                m.header.frame_id = "map";
                m.header.stamp = rclcpp::Time(0);
                m.ns = "gt_wind";
                m.id = static_cast<int>(i);
                m.type = Marker::ARROW;
                m.action = Marker::ADD;

                // Set pose
                m.pose.position.x = cx;
                m.pose.position.y = cy;
                m.pose.orientation = Utils::createQuaternionMsgFromYaw(direction);
                // shape
                m.scale.x = cell_size * module / max_module; // arrow length,
                m.scale.y = 0.03;      // arrow width
                m.scale.z = 0.05;      // arrow height
                // color -> must normalize to [0-199]
                Utils::get_arrow_color(module, max_module, m.color.r, m.color.g, m.color.b);
                m.color.a = 1.0;       // transparency
                // Add marker to Array
                wind_gt_array.markers.push_back(m);
            }

            if (verbose)
                RCLCPP_INFO(get_logger(), "[Cvalgt] Publishing Ground-Truth Wind Map with %zu markers.", wind_gt_array.markers.size());

            wind_gt_array_pub->publish(wind_gt_array);
        }
    }
    catch (const std::exception& e)
    {
        RCLCPP_ERROR(get_logger(), "[Cvalgt] Exception caught while reading CFD ground-truth wind map: '%s'", e.what());
    }
}


void Cvalgt::SimulateFixedWindObservations()
{
    // Insert fixed observations from the GT map (for testing purposes)

    Eigen::Vector2i dimensions = gmrf_map->map_size();
    double N = dimensions.x()*dimensions.y() - 1;
    std::vector<size_t> observation_indices;
    
    // Left vertical Inlet (10x6 scenarios)
    /*
    RCLCPP_WARN(get_logger(), "[Cvalgt] Using fixed observations at 10x6 Inlet.");
    for (size_t row = 11; row <= 14; ++row)
    {
        size_t i = row * dimensions.x() + 1;
        observation_indices.push_back(i);
    }
    
    */
    // Exp_C scenario: Set all the 6 inlets/outlets
    std::vector<int> activelets = {1, 2, 3, 4, 5 ,6};
    RCLCPP_WARN(get_logger(), "[Cvalgt] Using fixed observations at Exp_C Inlet.");
    
    // in/outlet 1 (top-left)
    if (std::find(activelets.begin(), activelets.end(), 1) != activelets.end())
    {
        for (size_t col = 14; col <= 17; ++col)
        {
            size_t i = (dimensions.y() - 2) * dimensions.x() + col;
            observation_indices.push_back(i);
        }
    }
   // in/outlet 2 (right-top)
   if (std::find(activelets.begin(), activelets.end(), 2) != activelets.end())
    {
        for (size_t row = (dimensions.y()-8); row <= (dimensions.y()-4); ++row)
        {
            size_t i = row * dimensions.x() + dimensions.x()-2;
            observation_indices.push_back(i);
        }
    }
    // in/outlet 3 (left-bottom)
    if (std::find(activelets.begin(), activelets.end(), 3) != activelets.end())
    {
        for (size_t row = 11; row <= 14; ++row)
        {
             size_t i = row * dimensions.x() + 1;
            observation_indices.push_back(i);
        }
    }
    // in/outlet 4 (right-middle)
    if (std::find(activelets.begin(), activelets.end(), 4) != activelets.end())
    {
        for (size_t row = 14; row <= 17; ++row)
        {
            size_t i = row * dimensions.x() + dimensions.x()-2;
            observation_indices.push_back(i);
        }
    }
    // in/outlet 5 (bottom-left)
    if (std::find(activelets.begin(), activelets.end(), 5) != activelets.end())
    {
        for (size_t col = 5; col <= 10; ++col)
        {
            size_t i = 1 * dimensions.x() + col;
            observation_indices.push_back(i);
        }
    }
    // in/outlet 6 (bottom-right)
    if (std::find(activelets.begin(), activelets.end(), 6) != activelets.end())
    {    
        for (size_t col = 29; col <= 35; ++col)
        {
            size_t i = 1 * dimensions.x() + col;
            observation_indices.push_back(i);
        }
    }
    

    // Create observations in the GMRF map from the GT map
    for (size_t idx : observation_indices)
    {
        // Ensure index is within bounds of the GT map
        if (idx >= gt_map.size())
        {
            RCLCPP_WARN(get_logger(), "[Cvalgt] Observation index %zu is out of bounds for GT map size %zu. Skipping this observation.", idx, gt_map.size());
            continue;
        }
       
        // Ensure cell is free
        if (!gmrf_map->is_cell_free(idx))
        {
            RCLCPP_WARN(get_logger(), "[Cvalgt] Observation index %zu corresponds to a non-free cell. Skipping this observation.", idx);
            continue;
        }
        
        // Read GT wind at that cell
        double wind_speed_x = gt_map[idx].x;
        double wind_speed_y = gt_map[idx].y;
        double module = sqrt(pow(wind_speed_x,2) + pow(wind_speed_y,2));
        double direction = atan2(wind_speed_y, wind_speed_x);

        // Cell center coordinates
        double x_pos_meters, y_pos_meters;
        x_pos_meters = gmrf_map->map_dimensions_meters()[0] + (idx % dimensions.x() + 0.5) * cell_size;
        y_pos_meters = gmrf_map->map_dimensions_meters()[2] + (idx / dimensions.x() + 0.5) * cell_size;

        // Insert observation in GRMF
        gmrf_map->insertObservation_GMRF(module, direction, observation_var_wind_speed, observation_var_wind_direction, x_pos_meters, y_pos_meters);        
    }
}


void Cvalgt::SimulateWindObservations(size_t N_obs, bool remove_old_observations)
{
    // 1. Setup Random Number Generator (RNG)
    // std::random_device provides a non-deterministic seed (best practice)
    // std::mt19937 is a fast, high-quality Mersenne Twister engine
    std::mt19937 gen(std::chrono::system_clock::now().time_since_epoch().count());
    
    // 2. Define Distribution
    // std::uniform_int_distribution ensures a uniform probability
    // for all integers in the range [0, N] (inclusive).
    Eigen::Vector2i dimensions = gmrf_map->map_size();
    double N = dimensions.x()*dimensions.y() - 1;
    std::uniform_int_distribution<> distrib(0, N);
    
    // Clear previous observations if requested
    if (remove_old_observations)
    {
        gmrf_map->clearObservations_GMRF();
    }

    // Avoid sampling the same cell multiple times
    std::vector<int> observed_cells;
    gmrf_map->getObservationsIdx(observed_cells);
    

    // Add N random observations from the GT map
    for (size_t i = 0; i < N_obs; ++i)
    {
        // 3. Generate and Return the Random Number
        int idx = distrib(gen);

        // Ensure we don't sample the same cell multiple times
        while (std::find(observed_cells.begin(), observed_cells.end(), idx) != observed_cells.end())
        {
            idx = distrib(gen);
        }

        // Ensure cell is free
        if (!gmrf_map->is_cell_free(idx))
        {
            --i; // If cell is not free, we don't count this iteration and try again
            continue;
        }        
        
        // Add to observed cells list
        observed_cells.push_back(idx);
        
        // Read GT wind at that cell
        double wind_speed_x = gt_map[idx].x;
        double wind_speed_y = gt_map[idx].y;
        double module = sqrt(pow(wind_speed_x,2) + pow(wind_speed_y,2));
        double direction = atan2(wind_speed_y, wind_speed_x);

        // Add gaussian noise to the observation
        std::normal_distribution<> noise_speed(0.0, sqrt(observation_var_wind_speed));
        std::normal_distribution<> noise_direction(0.0, sqrt(observation_var_wind_direction));

        // Cell center coordinates
        double x_pos_meters, y_pos_meters;
        x_pos_meters = gmrf_map->map_dimensions_meters()[0] + (idx % dimensions.x() + 0.5) * cell_size;
        y_pos_meters = gmrf_map->map_dimensions_meters()[2] + (idx / dimensions.x() + 0.5) * cell_size;

        // Insert observation in GRMF
        gmrf_map->insertObservation_GMRF(module + noise_speed(gen), direction + noise_direction(gen), observation_var_wind_speed, observation_var_wind_direction, x_pos_meters, y_pos_meters);
        if (verbose)
            RCLCPP_INFO(get_logger(), "[Cvalgt] Inserting observation %zu: speed=%.2f m/s, direction=%.2f rad, pos=(%.2f, %.2f) m", 
                        i, module, direction, x_pos_meters, y_pos_meters);
    }
}


void Cvalgt::publishMaps()
{
    if (!visualize_gmrf)
        return;
    
    visualization_msgs::msg::MarkerArray wind_array;
    visualization_msgs::msg::Marker wind_std_array;
    Utils::createWindMarkerArrayFromGMRF(*gmrf_map, "map", wind_array, wind_std_array);
    wind_array_pub->publish(wind_array);
    wind_std_array_pub->publish(wind_std_array);
}



void Cvalgt::update()
{
    // MAP estimation (with current lambdas)
    gmrf_map->MAP_estimation_GMRF(num_iterations_MAP);

    // Estimate uncertainty on final MAP estimation
    gmrf_map->computeUncertainty_GMRF();
}


std::vector<double> Cvalgt::compute_performance_metrics(const std::string& metric) const
{
    try
    {
        std::vector<double> residuals;
        Eigen::Vector2i dimensions = gmrf_map->map_size();
        size_t N = dimensions.x()*dimensions.y();
        size_t N2 = gt_map.size();
        double Nf = std::max(1.0, static_cast<double>(free_cell_indices.size()));
        if (N != N2) {
            RCLCPP_ERROR(get_logger(), "[gmrf-validation] ERROR: GMRF map size and GT map size do not match!");
            return {};
        }

        using visualization_msgs::msg::Marker;
        Marker metric_marker;
        if (visualize_gmrf)
        {
            // Visualization Metric is a single marker of type POINTS
            metric_marker.header.frame_id = "map";
            metric_marker.header.stamp = rclcpp::Time(0);
            metric_marker.ns = "gmrf_metric";
            metric_marker.id = 0;
            metric_marker.type = Marker::POINTS;
            metric_marker.action = Marker::ADD;
            // shape
            metric_marker.scale.x = cell_size;
            metric_marker.scale.y = cell_size;

            // Add one POINTS marker per cell
            metric_marker.points.reserve(N);
            metric_marker.colors.reserve(N);
        }

        // METRICS: AE, ME, RSE, NSP, NLPD
        double AE = 0.0;                // Angular Error (0-1)
        double ME = 0.0;                // Module Error or Speed difference (only magnitude) [0,inf]    
        double RSE = 0.0;               // Root Square Error (0-inf)
        double NSP = 0.0;               // Normalized Scalar Product [-1,1] -> [0,2] for optimization
        double NLPD = 0.0;              // Negative Log-Predictive Density (NLPD) (-inf, inf) -> [0, inf] for optimization

        double sum_AE = 0.0;                // Angular Error (0-1)
        double sum_ME = 0.0;                // To accumulate speed errors for average
        double sum_squared_error = 0.0;     // for RMSE (0-inf)
        double sum_gt_magnitudes = 0.0;     // To normalize the RMSE
        double sum_NSP = 0.0;               // Scalar product
        double sum_NLPD = 0.0;              // Negative Log-Predictive Density (NLPD) (-inf, inf) -> [0, inf] for optimization

        // Iterate over free cells!
        for (int i: free_cell_indices)
        {
            // Skip obstacles and boundary cells (if any) - inlets and outlets            
            //if (gmrf_map->is_cell_free(i) && !gmrf_map->is_cell_boundary(i))
            {
                // Get GMRF estimation at cell i (polars)
                WindVector est_wind = gmrf_map->getEstimation(i);
                double est_module = est_wind.module;
                double est_direction = est_wind.angle;
                double est_varModule = est_wind.varMod;
                double est_varDirection = est_wind.varAngle;
                double est_covModuleDirection = est_wind.covModAngle;
                // Cartesian (for plotting and metrics)
                double est_x = est_wind.x;
                double est_y = est_wind.y;
                double est_varX = est_wind.varX;
                double est_varY = est_wind.varY;
                double est_covXY = est_wind.covXY;
                
                // Get GT wind at cell i
                WindVectorXY gt_wind = gt_map[i];
                double gt_x = gt_wind.x;
                double gt_y = gt_wind.y;
                double gt_module = sqrt(pow(gt_wind.x,2) + pow(gt_wind.y,2));
                double gt_direction = atan2(gt_wind.y, gt_wind.x);
                
                // ===================================================
                // METRIC1 - AE: Angular Error (only direction) [0,1]
                // ===================================================
                try
                {
                    double AE = 0.0;
                    // Case 1: Both are zero (Perfect match of a dead-calm wind zone)
                    if (gt_module < 1e-6 && est_module < 1e-6) 
                    {
                        AE = 0.0;
                    }
                    // Case 2: GT has wind, but the model predicted zero wind
                    else if (est_module < 1e-6) 
                    {
                        // The model is failing to predict existing wind.                        
                        AE = 1.0; // Apply maximum penalty to actively repel the optimizer from degenerate zero-fields!
                    }
                    // Case 3: GT has no wind, but the model hallucinated wind
                    else if (gt_module < 1e-6) 
                    {
                        // The model is inventing wind where there is none.
                        AE = 1.0; // Maximum penalty 
                    }
                    // Case 4: Normal vectors (Both > 1e-6)
                    else
                    {
                        // 1. Calculate the dot product
                        double dot_product = (est_x * gt_x) + (est_y * gt_y);

                        // 2. Calculate Cosine Similarity (ranges from -1 to 1)
                        double cosine_similarity = dot_product / (est_module * gt_module);

                        // Clamp values to the valid range [-1, 1] for acos stability due to floating point errors
                        cosine_similarity = std::max(-1.0, std::min(1.0, cosine_similarity));

                        // 4. Calculate the angle in radians (ranges from 0 to pi)
                        double angle_radians = std::acos(cosine_similarity);

                        // 5. Normalize the angle from [0, pi] range to [0, 1] range
                        AE = angle_radians / M_PI;
                    }
                    
                    // 6. Accumulate (for average later)
                    sum_AE += AE;
                }
                catch(const std::exception& e)
                {
                    std::cerr << "AE exception: " << e.what() << '\n';
                }

                // ===================================================
                // METRIC2: Speed difference (only magnitude) [0,inf]
                // ===================================================
                try
                {
                    ME = std::abs( gt_module - est_module );
                    sum_ME += ME;
                }
                catch(const std::exception& e)
                {
                    std::cerr << "Speed difference exception: " << e.what() << '\n';
                }

                // ===================================================
                // METRIC3: RMSE in Cartesian (module+direction) [0,inf]
                // ===================================================
                try
                {
                    // Squared Euclidean distance
                    double diff_x = gt_x - est_x;
                    double diff_y = gt_y - est_y;
                    double squared_error = (diff_x * diff_x + diff_y * diff_y);
                    RSE = std::sqrt(squared_error);
                    sum_squared_error += squared_error;

                    // Calculate the magnitude of the ground truth vector: ||GT_i|| for normalization
                    sum_gt_magnitudes += gt_module;
                }
                catch(const std::exception& e)
                {
                    std::cerr << "RMSE exception: " << e.what() << '\n';
                }


                // ===================================================
                // METRIC4: NSP - Scalar Product (module+direction) [-1,1] --> [0,2]
                // ===================================================
                try
                {
                    double dot_product = est_x * gt_x + est_y * gt_y;
                    double gt_magnitude_squared = gt_x * gt_x + gt_y * gt_y;
                    double est_magnitude_squared = est_x * est_x + est_y * est_y;
                    double c = 0.001;  // small constant to avoid division by zero

                    // NSP computation, range [-1, 1]
                    // 1: perfect alignment, -1: opposite direction, 0: orthogonal
                    // |1|: perfect match in module, |<1|: error in module
                    NSP = (2*dot_product + c) / (gt_magnitude_squared + est_magnitude_squared + c);
                    
                    // Adjust range to [0,2] for optimization.
                    NSP = (1-NSP); // 0 means perfect alignment, 2 means opposite direction

                    // Accumulate for average later
                    sum_NSP += NSP;
                }
                catch(const std::exception& e)
                {
                    std::cerr << "NSP exception: " << e.what() << '\n';
                }                    

                // ===================================================
                // METRIC5: NLPD - Negative Log-Predictive Density (uncertainty aware) [-inf,inf]
                // ===================================================
                try
                {
                    // Cartesian NLPD (full covariance Gaussian in X and Y)
                    //=========================================
                    // Calculate residuals
                    double dx = gt_x - est_x;
                    double dy = gt_y - est_y;

                    // Calculate the determinant of the covariance matrix
                    // Adding a tiny epsilon to ensure numerical stability/invertibility
                    double det = (est_varX * est_varY) - (est_covXY * est_covXY);
                    double epsilon = 1e-6;
                    if (det < epsilon) det = epsilon; 

                    // Compute the Mahalanobis Distance squared (the quadratic form)
                    // Formula: (1/det) * [dx, dy] * [varY, -covXY; -covXY, varX] * [dx; dy]
                    double mahalanobis_sq = (1.0 / det) * (
                        dx * dx * est_varY + 
                        dy * dy * est_varX - 
                        2.0 * dx * dy * est_covXY
                    );

                    // Calculate NLPD (-inf,inf)
                    NLPD = 0.5 * (mahalanobis_sq + std::log(det) + 2.0 * std::log(2.0 * M_PI));
                    
                    // Adjust range for optimization (minimization towards 0)
                    double tau = 10.0; 
                    double scaled_nlpd = NLPD / tau;
                    
                    // PREVENT INFINITY OVERFLOW
                    if (scaled_nlpd > 700.0) scaled_nlpd = 700.0; 
                    
                    NLPD = std::exp(scaled_nlpd);  //This is necessary for optimizaton

                    // Accumulate for average later
                    sum_NLPD += NLPD;
                }
                catch(const std::exception& e)
                {
                    std::cerr << "NLPD exception: " << e.what() << '\n';
                }

                // ===================================================
                // RESIDUAL FOR OPTIMIZATION AND VISUALIZATION
                // ===================================================
                double current_viz_val = 0.0;               
                if (metric == "AE") {
                    residuals.push_back(AE);
                    current_viz_val = AE;
                }
                else if (metric == "ME") {
                    residuals.push_back(ME);
                    current_viz_val = ME;
                }
                else if (metric == "RSE") {
                    residuals.push_back(RSE);
                    current_viz_val = RSE;
                }
                else if (metric == "NSP") {
                    residuals.push_back(NSP);
                    current_viz_val = NSP;
                }
                else if (metric == "NLPD") {
                    residuals.push_back(NLPD);
                    current_viz_val = NLPD;
                }
                else if (metric == "ME&AE") {
                    double combined_error = 0.7 * AE + 0.3 * ME; 
                    residuals.push_back(combined_error);
                    current_viz_val = combined_error;
                }
                
                // ===============================
                // VISUALIZATION MARKERS IN RVIZ2
                // ===============================
                if (visualize_gmrf)
                {
                    // Cell center coordinates
                    double cx = 0.0, cy = 0.0;
                    gmrf_map->id2xy_public(i, cx, cy);
                    
                    // Add one point at the cell center
                    geometry_msgs::msg::Point p;
                    p.x = cx;
                    p.y = cy;
                    p.z = 0.0;
                    metric_marker.points.push_back(p);
                    
                    std_msgs::msg::ColorRGBA color;
                    // color -> must normalize to [0-199]
                    double max_value = 1.0; // default max value for normalization
                    if (metric == "AE")
                        max_value = 1.0; // Angular Error is already normalized to [0,1]
                    else if (metric == "ME")
                        max_value = 1.0; // Assuming max speed error of 1 m/s for normalization
                    else if (metric == "RSE")
                        max_value = 1.0; // Assuming max RMSE of 1 m/s for normalization
                    else if (metric == "NSP")
                        max_value = 2.0; // NSP ranges from 0 to 2 after adjustment
                    else if (metric == "NLPD")
                        max_value = 10.0;
                    else if (metric == "ME&AE")
                        max_value = 1.0; // Assuming max combined error of 1 for normalization
                    
                    // CAREFUL --> If we are asking for an average residual, set a per-cell metric
                    else if (metric == "AAE")
                    {
                        current_viz_val = AE;
                        max_value = 1.0;
                    }
                    else if (metric == "AME")
                    {
                        current_viz_val = ME;
                        max_value = 1.0;
                    }
                    else if (metric == "NRMSE")
                    {
                        current_viz_val = RSE;
                        max_value = 1.0;
                    }
                    else if (metric == "ANSP")
                    {
                        current_viz_val = NSP;
                        max_value = 2.0;
                    }
                    else if (metric == "ANLPD")
                    {
                        current_viz_val = NLPD;
                        max_value = 10.0;
                    }
                    else if (metric == "AVERAGES")
                    {
                        // For the average metrics, we can visualize the NSP as a representative metric
                        current_viz_val = NSP;
                        max_value = 2.0;
                    }

                    // Color
                    Utils::get_arrow_color(current_viz_val, max_value, color.r, color.g, color.b);
                    color.a = 1.0;       // transparency
                    metric_marker.colors.push_back(color);
                }
            }            
        } // end for each free-cell  

        if (visualize_gmrf) {
            metric_pub->publish(metric_marker);
        }

        // ===============================
        // ACCUMULATED METRICS (AVERAGES)
        // ===============================
        double AAE = sum_AE / Nf;       // Average Angular Error [0,1]
        double AME = sum_ME / Nf;       // Average Module Error (speed difference) [0,inf]
        double RMSE = std::sqrt(sum_squared_error / Nf);    // Root Mean Square Error (RMSE) [0,inf]
        double NRMSE = (sum_gt_magnitudes > 1e-6) ? (RMSE / ( sum_gt_magnitudes / Nf )) : RMSE; // Normalized RMSE (NRMSE) [0,inf]
        double ANSP = sum_NSP / Nf;                         // Normalized Scalar Product [0,2]
        double ANLPD = sum_NLPD / Nf;                       // Average Negative Log-Predictive Density (NLPD) [0,inf]

        if (verbose) {
            RCLCPP_INFO(this->get_logger(), "[gmrf-validation] Average Angular Error (rad) = %.2f", AAE);
            RCLCPP_INFO(this->get_logger(), "[gmrf-validation] Average Module Error (m/s) = %.2f", AME);
            RCLCPP_INFO(this->get_logger(), "[gmrf-validation] NRMSE (m/s)= %.2f", NRMSE);
            RCLCPP_INFO(this->get_logger(), "[gmrf-validation] Average Normalized Scalar Product = %.2f", ANSP);
            RCLCPP_INFO(this->get_logger(), "[gmrf-validation] Average exp(NLPD/10) = %.2f", ANLPD);
        }
        
        if (metric == "AAE") {
            residuals.push_back(AAE);
        }
        else if (metric == "AME") {
            residuals.push_back(AME);
        }
        else if (metric == "NRMSE") {
            residuals.push_back(NRMSE);
        }
        else if (metric == "ANSP") {
            residuals.push_back(ANSP);
        }
        else if (metric == "ANLPD") {
            residuals.push_back(ANLPD);
        }
        else if (metric == "AVERAGES") {    // used for logging
            residuals.push_back(AAE);
            residuals.push_back(AME);
            residuals.push_back(NRMSE);
            residuals.push_back(ANSP);
            residuals.push_back(ANLPD);
        }
        
        // return residuals for optimization (per-cell) or average metric (single value)
        return residuals;
        
    }
    catch (std::exception e)
    {
        std::cerr << "=============================================================" << std::endl;
        std::cerr << "[compute_performance_metrics] EXCEPTION: " << e.what() << std::endl;
        std::cerr << "=============================================================" << std::endl;
        return {1.0, 999.0, 0.0, 999.0};
    }
}


void Cvalgt::saveGMRFEstimationToCSV(const std::string& file_name)
{
    // Header
    Eigen::Vector2i dimensions = gmrf_map->map_size();
    size_t N = dimensions.x()*dimensions.y();
    std::vector<int> observed_cells;
    gmrf_map->getObservationsIdx(observed_cells);

    std::ofstream csv_file;
    csv_file.open(file_name);

    // Save header
    csv_file << "Dimensions_x," << dimensions.x() << "\n";
    csv_file << "Dimensions_y," << dimensions.y() << "\n";
    csv_file << "Cell_size," << cell_size << "\n";
    csv_file << "N_observations," << observed_cells.size() << "\n";
    csv_file << "Obs_idx,";
    for (size_t i = 0; i < observed_cells.size(); ++i)        
        csv_file << observed_cells[i] << (i < observed_cells.size() - 1 ? "," : "\n");
    csv_file << "cell_index,"
             << "gmrf_wind_x,"
             << "gmrf_wind_y,"
             << "gmrf_var_x,"
             << "gmrf_var_y,"
             << "gmrf_cov_xy,"
             << "gt_wind_x,"
             << "gt_wind_y\n";

    for (size_t i = 0; i < N; ++i)
    {
        // Skip occupied cells
        if (gmrf_map->is_cell_free(i))
        {
            // Get GMRF estimation at cell i
            WindVector est_wind = gmrf_map->getEstimation(i);
            double est_x = est_wind.x;
            double est_y = est_wind.y;
            double est_varX = est_wind.varX;
            double est_varY = est_wind.varY;
            double est_covXY = est_wind.covXY;

            // Get GT wind at cell i
            WindVectorXY gt_wind = gt_map[i];
            double gt_x = gt_wind.x;
            double gt_y = gt_wind.y;

            // Save to CSV
            csv_file << i << ","
                     << est_x << ","
                     << est_y << ","
                     << est_varX << ","
                     << est_varY << ","
                     << est_covXY << ","
                     << gt_x << ","
                     << gt_y << "\n";
        }
        else
        {
            // For occupied cells, save NaN
            csv_file << i << ","
                     << "NaN,NaN,NaN,NaN,NaN,NaN\n";
        }
    }

    csv_file.close();
    RCLCPP_INFO(get_logger(), "[gmrf-validation] GMRF estimation saved to CSV file: %s", file_name.c_str());
}


void Cvalgt::update_lambdas(double lambda_advection, double lambda_mass, double lambda_diffusion, double lambda_obstacles)
{
    gmrf_map->update_lambdas(lambda_advection, lambda_mass, lambda_diffusion, lambda_obstacles);
}


void Cvalgt::read_lambdas(double &lambda_advection, double &lambda_mass, double &lambda_diffusion, double &lambda_obstacles)
{
    gmrf_map->read_lambdas(lambda_advection, lambda_mass, lambda_diffusion, lambda_obstacles);
}

void Cvalgt::clearEstimation()
{
    gmrf_map->clearEstimation();
}


struct Stats {
    double mean;
    double variance;
    double std_dev;
};


Stats calculate_stats(const std::vector<double>& data) 
{
    if (data.empty()) return {0.0, 0.0, 0.0};

    // 1. Calcular la media de forma eficiente
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    double mean = sum / data.size();

    // 2. Calcular la suma de cuadrados de las diferencias
    double sq_diff_sum = std::accumulate(data.begin(), data.end(), 0.0,
        [mean](double acc, double x) {
            return acc + (x - mean) * (x - mean);
        });

    // 3. Calcular varianza y desviación
    double variance = sq_diff_sum / data.size();
    double std_dev = std::sqrt(variance);

    return {mean, variance, std_dev};
}


// Generates an LHS matrix of size [num_samples x num_dims]
std::vector<std::vector<double>> generateLHS(int num_samples, int num_dims, double min_val, double max_val, int seed = 42) {
    std::vector<std::vector<double>> lhs(num_samples, std::vector<double>(num_dims));
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    double interval_size = (max_val - min_val) / num_samples;

    for (int d = 0; d < num_dims; ++d) {
        // Create an ordered list of intervals and shuffle them
        std::vector<int> intervals(num_samples);
        std::iota(intervals.begin(), intervals.end(), 0);
        std::shuffle(intervals.begin(), intervals.end(), generator);

        for (int s = 0; s < num_samples; ++s) {
            // Pick a random point within the assigned interval
            double random_offset = dist(generator) * interval_size;
            lhs[s][d] = min_val + (intervals[s] * interval_size) + random_offset;
        }
    }
    return lhs;
}

//-----------------------------------------------------------------------------
//                                    MAIN
//----------------------------------------------------------------------------
int main(int argc, char** argv)
{
    // Initialize ROS2
    rclcpp::init(argc, argv);

    // Create and Init the GMRF-wind-validation node
    auto my_gmrf_map = std::make_shared<Cvalgt>();

    //=========================================
    // SELECT EXPERIMENT
    // 1- Optimize Lambda parameters (scenario, wind, Nobs)
    // 2- Simple GMRF estimation with fixed Lambdas
    // 3- Manual analysis (Dynamic reconfiguration) of the effect of each Lambda on the GMRF estimation
    //=========================================
    
    switch (my_gmrf_map->experiment_number)
    {
    case 1: //Optimize Lambda parameters (Fixed Observations -> at Inlets)
    {        
        // CERES configuration
        // ===================================================================
        // PHASE 1: ANSP (structural optimization of the GMRF)
        // ===================================================================
        ceres::Solver::Options options1;
        options1.linear_solver_type = ceres::DENSE_QR;
        options1.max_num_iterations = 100;
        options1.use_nonmonotonic_steps = true; // Using a non-monotonic steps can help "jump" over small ridges
        options1.minimizer_progress_to_stdout = false;
        options1.function_tolerance = 1e-3;
        options1.gradient_tolerance = 1e-3;
        options1.parameter_tolerance = 1e-4;

        // ===================================================================
        // PHASE 2: NLPD (log_k) uncertainty-aware optimization of the GMRF
        // ===================================================================
        ceres::Solver::Options options2;
        options2.linear_solver_type = ceres::DENSE_QR;
        options2.max_num_iterations = 100;
        options2.use_nonmonotonic_steps = true; // Using a non-monotonic steps can help "jump" over small ridges
        options2.minimizer_progress_to_stdout = false;
        options2.function_tolerance = 1e-6;
        options2.parameter_tolerance = 1e-5;

        // ===================================================================
        // PRE-COMPUTE LATIN HYPERCUBE SAMPLES (coarse grid search for good initializations of the optimization)
        // ===================================================================
        int num_repetitions = 100;
        // We explore bounds [-3.0, 3.0] to get a wide variety of initializations.
        std::vector<std::vector<double>> lhs_samples = generateLHS(num_repetitions, 3, -3.0, 3.0, 12345);

        
        //=========================================
        // LOOP - NUM OBSERVATIONS
        //=========================================
        for (int i = 6; i <= 6 ; i += 1) 
        {
            std::cerr << "================== Optimizing for N_obs = " << i << " ==================\n";

            // Save to CSV file
            std::string filename_lambda = "Lambda_values_" + std::to_string(i) + "_obs_" + my_gmrf_map->optimization_metric + ".csv";
            std::ofstream file_lambda(filename_lambda);
            if (file_lambda.is_open()) 
            {
                // Header
                file_lambda << "Repetition,AAE,AME,NRMSE,ANSP,ANLPD,Lambda_Advection,Lambda_Mass,Lambda_Diffusion,Lambda_Obstacles\n";
                file_lambda.close();
            } else {
                std::cerr << "Error: No se pudo abrir el archivo para escribir." << std::endl;
            }

            //=========================================
            // LOOP - M TIMES (different initializations)
            //=========================================
            for (int repeat = 0; repeat < num_repetitions; ++repeat)
            {
                std::cerr << "---- Repetition " << repeat << " ----\n";
                // 1. Clear previous estimation (important to avoid bias in the optimization)
                my_gmrf_map->clearEstimation();

                // 2. Simulate fixed observations (inlets)
                my_gmrf_map->SimulateFixedWindObservations();

                // 3. Initialize parameters for this specific run
                std::string initi_mode = "LHS"; // "LHS" for Latin Hypercube Sampling, "random", "fixed"
                double alpha[4];
                double log_k[1]; 

                if (initi_mode == "random")
                {
                    // Randomize initial alpha and log_k for each repetition
                    std::default_random_engine generator(repeat); // Seed with repetition number for reproducibility
                    std::normal_distribution<double> distribution(0.0, 0.25); // mean=0, stddev=0.25
                    for (int j = 0; j < 4; ++j)
                        alpha[j] = distribution(generator);
                    alpha[2] = 0.0; // Keep alpha_diffusion fixed at 0 for all repetitions to maintain the same reference scale for the other lambdas
                    log_k[0] = std::log(1/my_gmrf_map->observation_var_wind_speed);
                }
                else if (initi_mode == "LHS")
                {
                    // Use Latin Hypercube Sampling for initialization
                    alpha[0] = lhs_samples[repeat][0];
                    alpha[1] = lhs_samples[repeat][1];
                    alpha[2] = 0.0;                     // anchor
                    alpha[3] = lhs_samples[repeat][2]; 
                    log_k[0] = std::log(1/my_gmrf_map->observation_var_wind_speed); //controls uncertainty
                }
                else
                {
                    // Fixed initialization (same for all repetitions)
                    alpha[0] = 1.0; // alpha_advection
                    alpha[1] = 1.0; // alpha_mass
                    alpha[2] = 0.0; // alpha_diffusion (reference)
                    alpha[3] = 1.0; // alpha_obstacles
                    log_k[0] = std::log(1/my_gmrf_map->observation_var_wind_speed);
                }

                //-------------------------------------------------
                // 4. PHASE-1: ERROR IN THE STRUCTURE OF THE GMRF
                //-------------------------------------------------
                ceres::Problem problem_phase1;

                // Get the number residuals
                int num_residuals = 0;
                if (my_gmrf_map->optimization_metric == "AE" || my_gmrf_map->optimization_metric == "ME" || 
                    my_gmrf_map->optimization_metric == "RSE" || my_gmrf_map->optimization_metric == "NSP" || 
                    my_gmrf_map->optimization_metric == "NLPD" || my_gmrf_map->optimization_metric == "ME&AE") 
                {
                    // One residual per valid cell
                    num_residuals = my_gmrf_map->free_cell_indices.size();
                } 
                else if (my_gmrf_map->optimization_metric == "AVERAGES") 
                {
                    num_residuals = 5; 
                } 
                else 
                {
                    num_residuals = 1; // AAE, AME, NRMSE, etc.
                }
                
                // Create a DYNAMIC cost function
                ceres::DynamicNumericDiffCostFunction<CostFunctor>* cost_function1 =
                    new ceres::DynamicNumericDiffCostFunction<CostFunctor>(
                        new CostFunctor(my_gmrf_map.get(), my_gmrf_map->optimization_metric)
                    );

                // Define the parameter blocks (4 alphas, 1 log_k)
                cost_function1->AddParameterBlock(4); 
                cost_function1->AddParameterBlock(1);

                // Define the number of residuals
                cost_function1->SetNumResiduals(num_residuals);

                // Add it to the problem
                problem_phase1.AddResidualBlock(cost_function1, nullptr, alpha, log_k);
                
                // Set fixed and varable parameters for this phase
                problem_phase1.SetParameterBlockVariable(alpha);
                problem_phase1.SetParameterBlockConstant(log_k);  // Keep k constant if not optimizing NLPD 

                // Dock one alpha parameter to avoid scale ambiguity (only optimize 3 alphas)
                std::vector<int> constant_alpha_indices;
                constant_alpha_indices.push_back(2); // Fix alpha_diffusion to 0 (lambda_diffusion = k*e^0 = k) and optimize the relative importance of the other lambdas with respect to diffusion
#if CERES_VERSION_MAJOR >= 2 && CERES_VERSION_MINOR >= 2
                ceres::SubsetManifold* alpha_parameterization =
                    new ceres::SubsetManifold(4, constant_alpha_indices);
                problem_phase1.SetManifold(alpha, alpha_parameterization);
#else
                ceres::SubsetParameterization* alpha_parameterization =
                    new ceres::SubsetParameterization(4, constant_alpha_indices);
                problem_phase1.SetParameterization(alpha, alpha_parameterization);
#endif

                // Set reasonable bounds to avoid numerical issues and extreme values
                for (int i : {0, 1, 3}) {
                    problem_phase1.SetParameterLowerBound(alpha, i, -5.0);
                    problem_phase1.SetParameterUpperBound(alpha, i, 5.0);
                }                

                // Solve Phase 1 (optimization)
                ceres::Solver::Summary summary1;
                ceres::Solve(options1, &problem_phase1, &summary1);
                std::cerr << summary1.BriefReport() << "\n";

                if (false)
                {
                    //--------------------
                    // 5. PHASE-2: NLPD
                    //--------------------
                    ceres::Problem problem_phase2;
                    ceres::CostFunction* cost_function2 =
                        new ceres::NumericDiffCostFunction<CostFunctor, ceres::CENTRAL, 1, 4, 1>(
                            new CostFunctor(my_gmrf_map.get(), "anlpd")
                        );
                    problem_phase2.AddResidualBlock(cost_function2, nullptr, alpha, log_k);
                    // Set fixed and varable parameters for this phase
                    problem_phase2.SetParameterBlockConstant(alpha);    // Keep alphas constant in this phase to focus on optimizing k
                    problem_phase2.SetParameterBlockVariable(log_k);    // Ensure log_k is variable in this phase
                    problem_phase2.SetParameterLowerBound(log_k, 0, -5.0);  // Min log_k = -5 --> k = 0.0067 (very small lambdas) 
                    problem_phase2.SetParameterUpperBound(log_k, 0, 5.0);   // Max log_k = 5 --> k = 148.4 (very large lambdas)

                    ceres::Solver::Summary summary2;
                    ceres::Solve(options2, &problem_phase2, &summary2);
                    std::cerr << summary2.BriefReport() << "\n";
                    
                    if (!summary1.IsSolutionUsable() || !summary2.IsSolutionUsable())
                    {
                        std::cerr << "[WARNING] Optimización inestable en la repetición " 
                                << repeat << ". Descartando resultados para evitar outliers.\n";
                        continue; 
                    }
                }
            

                // 6. Optimized Results
                // ---------------------
                double final_k = std::exp(log_k[0]);
                double sum_exp = std::exp(alpha[0]) + std::exp(alpha[1]) + std::exp(alpha[2]) + std::exp(alpha[3]);
                double lambda[4];
                for (int i = 0; i < 4; ++i){
                    lambda[i] = final_k * (std::exp(alpha[i]) / sum_exp);
                }
                
                // 7. Make sure uncertainty is computed one last time with final values
                my_gmrf_map->update_lambdas(lambda[0], lambda[1], lambda[2], lambda[3]);
                my_gmrf_map->update();

                // 7. Performance metrics for logging (optimal values)
                std::vector<double> metrics = my_gmrf_map->compute_performance_metrics("AVERAGES");
                
                // 8. Append results to the CSV file
                std::ofstream file_append(filename_lambda, std::ios::app);
                if (file_append.is_open()) 
                {
                    file_append << repeat << "," 
                                << metrics[0] << ","    // AAE
                                << metrics[1] << ","    // AME
                                << metrics[2] << ","    // NRMSE
                                << metrics[3] << ","    // ANSP
                                << metrics[4] << ","    // ANLPD
                                << lambda[0] << "," 
                                << lambda[1] << ","
                                << lambda[2] << ","
                                << lambda[3] << "\n";
                    file_append.close();
                } else {
                    std::cerr << "Error: No se pudo abrir el archivo para escribir." << std::endl;
                }
            } // end loop M times.
        } // end loop N observations
        break;
    }
    case 2:     // Manual tunning of lambda parameters (dynamic reconfigure)
    {
        // Simulate Observations   
        my_gmrf_map->SimulateFixedWindObservations();
        
        // Start a background thread to handle ROS callbacks
        rclcpp::executors::MultiThreadedExecutor executor;
        executor.add_node(my_gmrf_map);
        std::thread spin_thread([&executor]() {
            executor.spin(); // This will handle rqt_reconfigure instantly!
        });

        rclcpp::Rate rate(1);
        while (rclcpp::ok())
        {
            // Update Lambda values (from parameter server)
            double lambda_advection = my_gmrf_map->get_parameter("GMRF_lambdaPrior_advection").as_double();
            double lambda_mass_conservation = my_gmrf_map->get_parameter("GMRF_lambdaPrior_mass_conservation").as_double();
            double lambda_diffusion = my_gmrf_map->get_parameter("GMRF_lambdaPrior_diffusion").as_double();
            double lambda_obstacles = my_gmrf_map->get_parameter("GMRF_lambdaPrior_obstacles").as_double();
            my_gmrf_map->update_lambdas(lambda_advection, lambda_mass_conservation, lambda_diffusion, lambda_obstacles);

            // Estimate MAP + Uncertainty
            my_gmrf_map->update();
            my_gmrf_map->publishMaps();

            // Compute performance metrics
            std::vector<double> metrics = my_gmrf_map->compute_performance_metrics("AVERAGES");
            
            // Display metrics in the console
            RCLCPP_INFO(my_gmrf_map->get_logger(), "[gmrf-validation] AAE: %.2f rad, AME: %.2f rad, RMSE: %.2f m/s, ANSP: %.2f, ANLPD: %.2f", 
                        metrics[0], metrics[1], metrics[2], metrics[3], metrics[4]);

            // rclcpp::spin_some(my_gmrf_map); --> Spin is handled in the background thread
            rate.sleep();

            // Save GMRF estimation to CSV file (for debugging/visualization purposes)
            if (true) {
                std::string filename_csv = "gmrf_estimation.csv";
                my_gmrf_map->saveGMRFEstimationToCSV(filename_csv);
                //RCLCPP_INFO(my_gmrf_map->get_logger(), "[gmrf-validation] GMRF estimation saved to CSV file: %s", filename_csv.c_str());
            }
        }

        // Clean up the thread when exiting
        rclcpp::shutdown();
        if (spin_thread.joinable()) {
            spin_thread.join();
        }
        break;
    }
    case 3:     //Increasing number of observations
    {
        // Start a background thread to handle ROS callbacks
        rclcpp::executors::MultiThreadedExecutor executor;
        executor.add_node(my_gmrf_map);
        std::thread spin_thread([&executor]() {
            executor.spin(); // This will handle rqt_reconfigure instantly!
        });

        // Set Lambda values (from parameter server)
        double lambda_advection = my_gmrf_map->get_parameter("GMRF_lambdaPrior_advection").as_double();
        double lambda_mass_conservation = my_gmrf_map->get_parameter("GMRF_lambdaPrior_mass_conservation").as_double();
        double lambda_diffusion = my_gmrf_map->get_parameter("GMRF_lambdaPrior_diffusion").as_double();
        double lambda_obstacles = my_gmrf_map->get_parameter("GMRF_lambdaPrior_obstacles").as_double();
        my_gmrf_map->update_lambdas(lambda_advection, lambda_mass_conservation, lambda_diffusion, lambda_obstacles);

        // Experiment folder
        std::string experiment_folder = "experiment2";
        std::filesystem::create_directories(experiment_folder);

        int num_repetitions = 100;
        for (int repeat = 0; repeat < num_repetitions; ++repeat)
        {
            // Clear GMRF estimation to avoid bias from previous runs
            my_gmrf_map->clearEstimation();
            bool clear_observations = true;

            // Create repetition folder
            std::string repetition_folder = experiment_folder + "/repetition_" + std::to_string(repeat);
            std::filesystem::create_directories(repetition_folder);
            
            int num_observations = 0;
            int max_observations = 200;
            for (num_observations = 0; num_observations <= max_observations; num_observations += 5)
            {
                RCLCPP_INFO(my_gmrf_map->get_logger(), "================== Simulating %d Observations ==================", num_observations);
                
                // Add new Observations
                my_gmrf_map->SimulateWindObservations(5, clear_observations);
                clear_observations = false; // Only clear observations for the first batch
                
                // Estimate MAP + Uncertainty
                my_gmrf_map->update();            

                // Save GMRF estimation to CSV file
                std::string filename_csv = repetition_folder + "/gmrf_estimation_expC_" + std::to_string(num_observations) + "_obs.csv";
                my_gmrf_map->saveGMRFEstimationToCSV(filename_csv);
                
                // Compute performance metrics
                std::vector<double> metrics = my_gmrf_map->compute_performance_metrics("AVERAGES");
                
                // Display metrics in the console
                RCLCPP_INFO(my_gmrf_map->get_logger(), "[gmrf-validation] AAE: %.2f rad, AME: %.2f rad, RMSE: %.2f m/s, ANSP: %.2f, ANLPD: %.2f", 
                            metrics[0], metrics[1], metrics[2], metrics[3], metrics[4]);
            }
        }
        
        // Clean up the thread when exiting
        rclcpp::shutdown();
        if (spin_thread.joinable()) {
            spin_thread.join();
        }
        break;
    }

    case 4:     //Computational Time
    {
        // Start a background thread to handle ROS callbacks
        rclcpp::executors::MultiThreadedExecutor executor;
        executor.add_node(my_gmrf_map);
        std::thread spin_thread([&executor]() {
            executor.spin(); // This will handle rqt_reconfigure instantly!
        });

        // Set Lambda values (from parameter server)
        double lambda_advection = my_gmrf_map->get_parameter("GMRF_lambdaPrior_advection").as_double();
        double lambda_mass_conservation = my_gmrf_map->get_parameter("GMRF_lambdaPrior_mass_conservation").as_double();
        double lambda_diffusion = my_gmrf_map->get_parameter("GMRF_lambdaPrior_diffusion").as_double();
        double lambda_obstacles = my_gmrf_map->get_parameter("GMRF_lambdaPrior_obstacles").as_double();
        my_gmrf_map->update_lambdas(lambda_advection, lambda_mass_conservation, lambda_diffusion, lambda_obstacles);
       
        int num_repetitions = 100;
        for (int repeat = 0; repeat < num_repetitions; ++repeat)
        {
            // Clear GMRF estimation to avoid bias from previous runs
            my_gmrf_map->clearEstimation();
            
            // Add new Observations
            my_gmrf_map->SimulateWindObservations(60, true);
                
            // Estimate MAP + Uncertainty
            my_gmrf_map->update();
        }
        
        // Clean up the thread when exiting
        rclcpp::shutdown();
        if (spin_thread.joinable()) {
            spin_thread.join();
        }
        break;
    }

    default:
        break;
    }
       
    
    return 0;
}
