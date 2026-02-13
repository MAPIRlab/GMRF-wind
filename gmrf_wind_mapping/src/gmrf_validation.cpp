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

using namespace std::placeholders;


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
    GMRF_lambdaPrior_reg = declare_parameter<double>("GMRF_lambdaPrior_reg", 1.0);
    GMRF_lambdaPrior_flux_conservation = declare_parameter<double>("GMRF_lambdaPrior_flux_conservation", 1.0);
    GMRF_lambdaPrior_obstacles = declare_parameter<double>("GMRF_lambdaPrior_obstacles", 1.0);

    // Experiment
     experiment_number = declare_parameter<int>("experiment_number", 0);
    

    // Publishers
    //----------------------------------
    wind_array_pub = create_publisher<visualization_msgs::msg::MarkerArray>("wind_array_pub", 1);
    wind_std_array_pub = create_publisher<visualization_msgs::msg::Marker>("wind_std_array_pub", 1);
    
    rclcpp::QoS qos_profile(1); // Keep the history depth low (1 is often sufficient for static data)
    qos_profile.durability(RMW_QOS_POLICY_DURABILITY_TRANSIENT_LOCAL);
    wind_gt_array_pub = create_publisher<visualization_msgs::msg::MarkerArray>("wind_gt_array_pub", qos_profile);
    metric_pub = create_publisher<visualization_msgs::msg::Marker>("metric_pub", 1);
    
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

        for (int y = 0; y < height; y++)
            for (int x = 0; x < width; x++)
            {
                double value = mapImage.at<uint8_t>(y, x) / 255.;
                int8_t& occupancy = occupancy_map.data.at(width * (height - y - 1) + x);
                if (value > free_thresh)
                    occupancy = 0;
                else if (value < occupied_thresh)
                    occupancy = 100;
                else
                    occupancy = -1;
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
    occMap.resolution = occupancy_map.info.resolution;      // cell size (m)
    occMap.width = occupancy_map.info.width;        // number of cells in X direction
    occMap.height = occupancy_map.info.height;       // number of cells in Y direction
    occMap.origin_x = occupancy_map.info.origin.position.x;        // world coordinates of the origin of the map
    occMap.origin_y = occupancy_map.info.origin.position.y;        // world coordinates of the origin of the map
    
    // Create the GMRF-Map and initialize its Prior Factors
    gmrf_map = std::make_unique<CGMRF_map>(occMap, 
                                        cell_size, 
                                        GMRF_lambdaPrior_reg, 
                                        GMRF_lambdaPrior_flux_conservation,
                                        GMRF_lambdaPrior_obstacles,
                                        verbose,
                                        true // estimateTiming
                                        );
    RCLCPP_INFO(get_logger(), "[GMRF-validation] GMRF Initialized");
    module_init = true;
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
                if (module < 0.01) m.scale.x = cell_size/5;    // arrow length,
                else m.scale.x = cell_size;      // arrow length,
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


void Cvalgt::SimulateWindObservations(size_t N_obs)
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
    
    // Clear previous observations
    gmrf_map->clearObservations_GMRF();
    std::vector<size_t> observed_cells; // To keep track of which cells have been observed

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

        // Insert observation
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
    gmrf_map->MAP_estimation_GMRF();

    // Estimate uncertainty on final MAP estimation
    gmrf_map->computeUncertainty_GMRF();
}


std::array<double,4> Cvalgt::compute_performance_metrics() const
{
    try
    {
        RCLCPP_INFO(this->get_logger(), "[metrics] Computing performance metrics...");
        // Compare current GMRF estimation with GT map
        Eigen::Vector2i dimensions = gmrf_map->map_size();
        size_t N = dimensions.x()*dimensions.y();
        size_t N2 = gt_map.size();
        //RCLCPP_INFO(get_logger(), "[gmrf-validation] Map size: GMRF=%zu cells, GT=%zu cells", N, N2);
        if (N != N2)
        {
            RCLCPP_ERROR(get_logger(), "[gmrf-validation] ERROR: GMRF map size and GT map size do not match!");
            return {0.0, 0.0, 0.0};
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

        // METRICS
        double sum_AE = 0.0;                // Angular Error (0-1)
        double sum_squared_error = 0.0;     // for RMSE (0-inf)
        double sum_gt_magnitudes = 0.0;     // To normalize the RMSE
        double sum_normalized_scalar_products = 0.0;   // Scalar product
        double NLPD = 0.0;                  // Negative Log-Predictive Density (NLPD) (-inf, inf)
        double sum_NLPD = 0.0;              // Negative Log-Predictive Density (NLPD) (-inf, inf)

        RCLCPP_INFO(this->get_logger(), "[metrics] 4EachCell...");
        for (size_t i = 0; i < N; ++i)
        {
            // Skip occupied cells
            if (gmrf_map->is_cell_free(i))
            {
                //RCLCPP_INFO(this->get_logger(), "[gmrf-validation] LOOP i=%zu/%zu", i, N);

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
                // METRIC1: Angular Error (only direction) [0,1]
                // ===================================================
                try
                {
                    // 1. Calculate the dot product
                    double dot_product = est_x * gt_x + est_y * gt_y;

                    // 2. Calculate Cosine Similarity (ranges from -1 to 1)
                    double cosine_similarity = dot_product / (est_module * gt_module + 1e-6);

                    // Clamp values to the valid range [-1, 1] for acos stability due to floating point errors
                    cosine_similarity = std::max(-1.0, std::min(1.0, cosine_similarity));

                    // 4. Calculate the angle in radians (ranges from 0 to pi)
                    double angle_radians = std::acos(cosine_similarity);

                    // 5. Normalize the angle from [0, pi] range to [0, 1] range
                    double AE = angle_radians / M_PI;

                    // 6. Accumulate
                    sum_AE += AE;
                }
                catch(const std::exception& e)
                {
                    std::cerr << "AAE exception: " << e.what() << '\n';
                }

                // ===================================================
                // METRIC2: RMSE in Cartesian [0,inf]
                // ===================================================
                try
                {
                    // Squared Euclidean distance
                    double diff_x = gt_x - est_x;
                    double diff_y = gt_y - est_y;
                    double squared_error = (diff_x * diff_x + diff_y * diff_y);
                    sum_squared_error += squared_error;

                    // Calculate the magnitude of the ground truth vector: ||GT_i|| for normalization
                    sum_gt_magnitudes += std::sqrt(gt_x * gt_x + gt_y * gt_y);
                }
                catch(const std::exception& e)
                {
                    std::cerr << "RMSE exception: " << e.what() << '\n';
                }


                // ===================================================
                // METRIC3: NSP - Normalized Scalar Product [-1,1]
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
                    double normalized_dot_product = (2*dot_product + c) / (gt_magnitude_squared + est_magnitude_squared + c);
                    sum_normalized_scalar_products += normalized_dot_product;
                }
                catch(const std::exception& e)
                {
                    std::cerr << "NSP exception: " << e.what() << '\n';
                }                    

                // ===================================================
                // METRIC4: NLPD - Negative Log-Predictive Density (uncertainty aware) [-inf,inf]
                // ===================================================
                try
                {
                    // 4.1 Cartesian NLPD (full covariance Gaussian in X and Y)
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

                    // Calculate NLPD
                    double nlpd_cart = 0.5 * (mahalanobis_sq + std::log(det) + 2.0 * std::log(2.0 * M_PI));


                    // 4.2 weighthed Polar NLPD computation
                    //=========================================
                    // Residuals
                    double dMod = gt_module - est_module;
                    double dAng = gt_direction - est_direction;

                    // Normalize angle difference to [-PI, PI] to avoid jumps at the boundary
                    while (dAng > M_PI) dAng -= 2.0 * M_PI;
                    while (dAng < -M_PI) dAng += 2.0 * M_PI;

                    // Covariance Matrix (using Jacobian to propagate from Cartesian to Polar)
                    double var_r = std::max(1e-6, est_varModule);
                    double var_theta = std::max(1e-6, est_varDirection);
                    double cov_r_theta = est_covModuleDirection;

                    // Calculate Determinant and Mahalanobis for Polar NLPD
                    double detP = (var_r * var_theta) - (cov_r_theta * cov_r_theta);
                    if (detP < 1e-6) detP = 1e-6;

                    double mahalanobis_sq_polar = (1.0 / detP) * (
                        dMod * dMod * var_theta + 
                        dAng * dAng * var_r - 
                        2.0 * dMod * dAng * cov_r_theta
                    );

                    double nlpd_polar = 0.5 * (mahalanobis_sq_polar + std::log(detP) + 2.0 * std::log(2.0 * M_PI));


                    // NLPD to consider (cartesian or polar)
                    //=========================================
                    NLPD = nlpd_cart;
                    
                    // Accumulate
                    sum_NLPD += NLPD;
                }
                catch(const std::exception& e)
                {
                    std::cerr << "NLPD exception: " << e.what() << '\n';
                }


                // ===============================
                // VISUALIZATION METRIC
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
                    Utils::get_arrow_color(NLPD, 1, color.r, color.g, color.b);
                    color.a = 1.0;       // transparency
                    metric_marker.colors.push_back(color);
                }
            }
        } // end for each cell

        if (visualize_gmrf)
            metric_pub->publish(metric_marker);


        // Average Angular Error
        double AAE = sum_AE / N;
        RCLCPP_INFO(this->get_logger(), "[gmrf-validation] Average Angular Error (rad) = %.2f", AAE);

        // RMSE & NRMSE
        double RMSE = std::sqrt(sum_squared_error / N);
        double NRMSE = RMSE / ( sum_gt_magnitudes / N);
        RCLCPP_INFO(this->get_logger(), "[gmrf-validation] RMSE (m/s)= %.2f", RMSE);

        // Normalized Scalar Product
        double ANSP = sum_normalized_scalar_products / N;
        RCLCPP_INFO(this->get_logger(), "[gmrf-validation] Average Normalized Scalar Product = %.2f", ANSP);

        // Average NLPD
        double ANLPD = sum_NLPD / N;
        RCLCPP_INFO(this->get_logger(), "[gmrf-validation] Average NLPD = %.2f", ANLPD);

        return {AAE, RMSE, ANSP, ANLPD};
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

    std::ofstream csv_file;
    csv_file.open(file_name);

    // Save header
    csv_file << "Dimensions_x," << dimensions.x() << "\n";
    csv_file << "Dimensions_y," << dimensions.y() << "\n";
    csv_file << "Cell_size," << cell_size << "\n";
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


void Cvalgt::update_lambdas(double lambda_reg, double lambda_flux, double lambda_obstacles)
{
    gmrf_map->update_lambdas(lambda_reg, lambda_flux, lambda_obstacles);
}

void Cvalgt::read_lambdas(double &lambda_reg, double &lambda_flux, double &lambda_obstacles)    
{
    gmrf_map->read_lambdas(lambda_reg, lambda_flux, lambda_obstacles);
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
    //=========================================
    
    switch (my_gmrf_map->experiment_number)
    {
    case 1: //Optimize Lambda parameters with NLPD (scenario, wind, Nobs)
    {
        // Initial values for the 3 Lambda parameters (Reg, Flux, Obstacles)
        double parameters[3] = {0.1, 0.1, 0.1};
        double ref_parameters[3] = {1, 10, 10};    // We use these values as baseline
        for (int j = 0; j < 3; ++j)
            parameters[j] = static_cast<double>(rand()) / RAND_MAX * 10; // Random value between 0 and 10
        my_gmrf_map->update_lambdas(parameters[0], parameters[1], parameters[2]);
        // CERES OPTIMIZER
        ceres::Problem problem;

        // <CostFunctor, Number of Residuals, Number of Parameters>
        ceres::CostFunction* cost_function =
            new ceres::NumericDiffCostFunction<CostFunctor, ceres::CENTRAL, 1, 3>(
                new CostFunctor(my_gmrf_map.get())
            );

        problem.AddResidualBlock(cost_function, nullptr, parameters);

        // Configure and Run the solver
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::DENSE_QR;
        options.max_num_iterations = 100;
        options.use_nonmonotonic_steps = true; // Using a non-monotonic steps can help "jump" over small ridges
        options.minimizer_progress_to_stdout = false;
        // Ensure lambdas are positive
        problem.SetParameterLowerBound(parameters, 0, 1e-6); // Lambda_Reg >= 0
        problem.SetParameterLowerBound(parameters, 1, 1e-6); // Lambda_Flux >= 0
        problem.SetParameterLowerBound(parameters, 2, 1e-6); // Lambda_Obstacles >= 0

        // DATA COLLECTION
        std::string filename_results = "gmrf_opt_metrics_lambda_Nobs.csv";
        std::ofstream file_results(filename_results);
        if (file_results.is_open()) 
        {
            // Header
            file_results << "N_Obs,AAE,RMSE,ANSP,NLPD,Lambda_Reg_mean,Lambda_Reg_std,Lambda_Flux_mean,Lambda_Flux_std,Lambda_Obstacles_mean,Lambda_Obstacles_std,GMRF_estimation_filename\n";
            file_results.close();
        } else {
            std::cerr << "Error: No se pudo abrir el archivo para escribir." << std::endl;
        }

        //=========================================
        // LOOP - NUM OBSERVATIONS
        //=========================================
        Eigen::Vector2i dimensions = my_gmrf_map->gmrf_map->map_size();
        int max_obs = dimensions.x() * dimensions.y();
        max_obs = 100; // Limit to 100 observations for testing
        for (int i = 50; i <= max_obs; i += 10) 
        {
            std::cerr << "================== Optimizing for N_obs = " << i << " ==================\n";
            
            // Repeat M times to get average performance
            std::vector<double> lambda1;
            std::vector<double> lambda2;
            std::vector<double> lambda3;

            // Save to CSV file
            std::string filename_lambda = "Lambda_values_obs_" + std::to_string(i) + ".csv";
            std::ofstream file_lambda(filename_lambda);
            if (file_lambda.is_open()) 
            {
                // Header
                file_lambda << "Repetition,AAE,RMSE,ANSP,NLPD,Lambda_Reg,Lambda_Flux,Lambda_Obstacles\n";
                file_lambda.close();
            } else {
                std::cerr << "Error: No se pudo abrir el archivo para escribir." << std::endl;
            }

            for (int repeat = 0; repeat < 20; ++repeat)
            {
                // Reset parameters on each iteration to avoid local minima (slower)
                for (int j = 0; j < 3; ++j)
                    parameters[j] = static_cast<double>(rand()) / RAND_MAX * 10; // Random value between 0 and 10
                my_gmrf_map->update_lambdas(parameters[0], parameters[1], parameters[2]);

                // Simulate N (random) observations (also clears previous observations)
                my_gmrf_map->SimulateWindObservations(i);

                // Optimize Lambdas for this scenario and observations
                ceres::Solver::Summary summary;
                ceres::Solve(options, &problem, &summary);

                // Store opt lambda values
                lambda1.push_back(parameters[0]);
                lambda2.push_back(parameters[1]);
                lambda3.push_back(parameters[2]);

                // Update GMRF with optimal lambda values
                my_gmrf_map->update_lambdas(parameters[0], parameters[1], parameters[2]);

                // Perform final estimation with optimized lambdas
                my_gmrf_map->update();

                // Compute performance metrics
                std::array<double, 4> metrics = my_gmrf_map->compute_performance_metrics();
                double AAE_opt = metrics[0];  // AAE is the first element
                double RMSE_opt = metrics[1]; // RMSE is the second element
                double ANSP_opt = metrics[2]; // ANSP is the third element
                double NLPD_opt = metrics[3]; // NLPD is the fourth element
                
                // Compute estimation with baseline parameters for comparison
                my_gmrf_map->update_lambdas(ref_parameters[0], ref_parameters[1], ref_parameters[2]);
                my_gmrf_map->update();
                std::array<double, 4> metrics_ref = my_gmrf_map->compute_performance_metrics();
                double AAE_ref = metrics_ref[0];  // AAE is the first element
                double RMSE_ref = metrics_ref[1]; // RMSE is the second element
                double ANSP_ref = metrics_ref[2]; // ANSP is the third element
                double NLPD_ref = metrics_ref[3]; // NLPD is the fourth element

                // Append results to the CSV file
                std::ofstream file_append(filename_lambda, std::ios::app);
                if (file_append.is_open()) 
                {
                    file_append << repeat << "," 
                                << AAE_opt << "," 
                                << RMSE_opt << "," 
                                << ANSP_opt << "," 
                                << NLPD_opt << "," 
                                << parameters[0] << "," 
                                << parameters[1] << ","
                                << parameters[2] << "\n";
                    // Baseline
                    file_append << repeat << "," 
                                << AAE_ref << "," 
                                << RMSE_ref << "," 
                                << ANSP_ref << "," 
                                << NLPD_ref << "," 
                                << ref_parameters[0] << "," 
                                << ref_parameters[1] << ","
                                << ref_parameters[2] << "\n";
                    file_append.close();
                } else {
                    std::cerr << "Error: No se pudo abrir el archivo para escribir." << std::endl;
                }
            } // end repeat

            // Compute mean/std of lambda values
            Stats stats1 = calculate_stats(lambda1);
            Stats stats2 = calculate_stats(lambda2);
            Stats stats3 = calculate_stats(lambda3);
                    
            // Update GMRF with mean lambda values
            parameters[0] = stats1.mean;
            parameters[1] = stats2.mean;
            parameters[2] = stats3.mean;            
            my_gmrf_map->update_lambdas(parameters[0], parameters[1], parameters[2]);

            // Perform final estimation with optimized lambdas
            my_gmrf_map->update();

            // Compute performance metrics
            std::array<double, 4> metrics = my_gmrf_map->compute_performance_metrics();
            double AAE_opt = metrics[0];  // AAE is the first element
            double RMSE_opt = metrics[1]; // RMSE is the second element
            double ANSP_opt = metrics[2]; // ANSP is the third element
            double NLPD_opt = metrics[3]; // NLPD is the fourth element

            // Append results to the CSV file
            std::ofstream file_append(filename_results, std::ios::app);
            if (file_append.is_open()) 
            {
                // Save estimation to file (for later visualization)
                std::string file_name = "gmrf_opt_estimation_obs_" + std::to_string(i) + ".csv";
                file_append << i << "," 
                            << AAE_opt << "," 
                            << RMSE_opt << "," 
                            << ANSP_opt << "," 
                            << NLPD_opt << "," 
                            << stats1.mean << "," << stats1.std_dev << ","
                            << stats2.mean << "," << stats2.std_dev << ","
                            << stats3.mean << "," << stats3.std_dev << ","
                            << file_name << "\n";
                file_append.close();

                // Save GMRF estimation to CSV
                my_gmrf_map->saveGMRFEstimationToCSV(file_name);
            } else {
                std::cerr << "Error: No se pudo abrir el archivo para escribir." << std::endl;
            }
        }
        break;
    }

    case 2:
    {
        // Simple GMRF estimation with fixed Lambdas
        // =========================================

        // Initial values for the 3 Lambda parameters (Reg, Flux, Obstacles)
        double parameters[3] = {my_gmrf_map->GMRF_lambdaPrior_reg, 
                                my_gmrf_map->GMRF_lambdaPrior_flux_conservation, 
                                my_gmrf_map->GMRF_lambdaPrior_obstacles};

        // Simulate 50 noisy random observations (wont change during the test)
        my_gmrf_map->SimulateWindObservations(50);

        // Update GMRF with new lambda values
        my_gmrf_map->update_lambdas(parameters[0], parameters[1], parameters[2]);
        
        // Perform MAP + Uncertainty estimation
        my_gmrf_map->update();

        // Save GMRF estimation to CSV
        std::string file_name = "gmrf_estimation_lambda_" + 
                                std::to_string(parameters[0]) + "_" + 
                                std::to_string(parameters[1]) + "_" + 
                                std::to_string(parameters[2]) + ".csv";
        my_gmrf_map->saveGMRFEstimationToCSV(file_name);

        break;

        // Iterate over each of the 3 parameter positions
        // Define the set of values you want to test
        std::vector<double> test_values = {0.01, 0.1, 1.0, 10.0, 100.0, 1000.0};

        for (int i = 0; i < 3; ++i) 
        {
            double original_val = parameters[i]; // Store original to reset later

            for (double val : test_values) 
            {
                // Update the specific parameter
                parameters[i] = val;

                // Update GMRF with new lambda values
                my_gmrf_map->update_lambdas(parameters[0], parameters[1], parameters[2]);
                
                // Perform MAP + Uncertainty estimation
                my_gmrf_map->update();

                // Save GMRF estimation to CSV
                std::string file_name = "gmrf_estimation_lambda_" + 
                                        std::to_string(parameters[0]) + "_" + 
                                        std::to_string(parameters[1]) + "_" + 
                                        std::to_string(parameters[2]) + ".csv";
                my_gmrf_map->saveGMRFEstimationToCSV(file_name);
            }

            // Reset parameter to base value before moving to the next index
            parameters[i] = original_val;
        }
        
        break;
    }
    default:
        break;
    }
       
    
    return 0;
}
