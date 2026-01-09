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

    // Lambda/weights for the different priors and observation factors
    GMRF_lambdaPrior_reg = declare_parameter<double>("GMRF_lambdaPrior_reg", 1.0);
    GMRF_lambdaPrior_flux_conservation = declare_parameter<double>("GMRF_lambdaPrior_flux_conservation", 1.0);
    GMRF_lambdaPrior_obstacles = declare_parameter<double>("GMRF_lambdaPrior_obstacles", 1.0);
    GMRF_lambdaObs = declare_parameter<double>("GMRF_lambdaObs", 1.0);

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

    // 4. Simulate K fixed observations (for testing)
    //SimulateWindObservations();

    // 5. Initialize metrics file
    metrics_filename = "gmrf_evaluation_metrics.csv";
    std::ofstream metrics_file;
    metrics_file.open(metrics_filename, std::ios_base::app); // append mode
    metrics_file << "GMRF_lambdaPrior_reg,"     // regularization
                 << "GMRF_lambdaPrior_flux,"    // mass conservation
                 << "GMRF_lambdaPrior_obs,"     // obstacles
                 << "GMRF_lambdaPrior_mea,"     // measurements
                 << "RMSE,"
                 << "AAE,"
                 << "NLPD\n";
    metrics_file.close();
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
                                        GMRF_lambdaObs,
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

    // Add N random observations from the GT map
    for (size_t i = 0; i < N_obs; ++i)
    {
        // 3. Generate and Return the Random Number
        int idx = distrib(gen);
        
        // Read GT wind at that cell
        double wind_speed_x = gt_map[idx].x;
        double wind_speed_y = gt_map[idx].y;
        double module = sqrt(pow(wind_speed_x,2) + pow(wind_speed_y,2));
        double direction = atan2(wind_speed_y, wind_speed_x);

        // Cell center coordinates
        double x_pos_meters, y_pos_meters;
        x_pos_meters = gmrf_map->map_dimensions_meters()[0] + (idx % dimensions.x() + 0.5) * cell_size;
        y_pos_meters = gmrf_map->map_dimensions_meters()[2] + (idx / dimensions.x() + 0.5) * cell_size;

        // Insert observation
        gmrf_map->insertObservation_GMRF(module, direction, x_pos_meters, y_pos_meters, GMRF_lambdaObs);
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
    /*
    // Update Lambda parameters (read parameter server)
    GMRF_lambdaPrior_reg = get_parameter("GMRF_lambdaPrior_reg").as_double();
    GMRF_lambdaPrior_flux_conservation = get_parameter("GMRF_lambdaPrior_flux_conservation").as_double();
    GMRF_lambdaPrior_obstacles = get_parameter("GMRF_lambdaPrior_obstacles").as_double();
    GMRF_lambdaObs = get_parameter("GMRF_lambdaObs").as_double();
    GMRF_lambdaObsLoss = get_parameter("GMRF_lambdaObsLoss").as_double();
    gmrf_map->update_lambdas(GMRF_lambdaPrior_reg, GMRF_lambdaPrior_flux_conservation, GMRF_lambdaPrior_obstacles, GMRF_lambdaObs);
    */     

    // Optimize Lambda hyperparameters (Log Marginal Likelihood maximization)
    // optimize_LML(bool performLMLOptimization, int maxIterations, double learningRate, double LML_threshold)
    // gmrf_map->optimize_LML(true, 1, 0.001, 0.1);

    // MAP estimation (with current lambdas)
    gmrf_map->MAP_estimation_GMRF();

    // Estimate uncertainty on final MAP estimation
    gmrf_map->computeUncertainty_GMRF();
}


std::array<double,3> Cvalgt::compute_performance_metrics() const
{
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
    double sum_NLPD = 0.0;              // Negative Log-Predictive Density (NLPD) (-inf, inf)

    for (size_t i = 0; i < N; ++i)
    {
        // Skip occupied cells
        if (gmrf_map->is_cell_free(i))
        {
            //RCLCPP_INFO(this->get_logger(), "[gmrf-validation] LOOP i=%zu/%zu", i, N);

            // Get GMRF estimation at cell i
            WindVector est_wind = gmrf_map->getEstimation(i);
            double est_module = est_wind.module;
            double est_direction = est_wind.direction;
            double est_x = est_wind.module * cos(est_wind.direction);
            double est_y = est_wind.module * sin(est_wind.direction);
            double est_std = est_wind.stdDev;
            double est_stdX = est_wind.stdDevX;
            double est_stdY = est_wind.stdDevY;

            // Get GT wind at cell i
            WindVectorXY gt_wind = gt_map[i];
            double gt_x = gt_wind.x;
            double gt_y = gt_wind.y;
            double gt_module = sqrt(pow(gt_wind.x,2) + pow(gt_wind.y,2));
            double gt_direction = atan2(gt_wind.y, gt_wind.x);
            
            // ===================================================
            // METRIC1: Angular Error (only direction) [0,1]
            // ===================================================
            // 1. Calculate the dot product
            double dot_product = est_x * gt_x + est_y * gt_y;

            // 2. Calculate the magnitude (L2 Norm) of each vector
            double mag_v1 = std::sqrt(est_x*est_x + est_y*est_y);
            double mag_v2 = std::sqrt(gt_x*gt_x + gt_y*gt_y);

            // 3. Calculate Cosine Similarity (ranges from -1 to 1)
            double cosine_similarity = dot_product / (mag_v1 * mag_v2 + 1e-6);

            // Clamp values to the valid range [-1, 1] for acos stability due to floating point errors
            cosine_similarity = std::max(-1.0, std::min(1.0, cosine_similarity));

            // 4. Calculate the angle in radians (ranges from 0 to pi)
            double angle_radians = std::acos(cosine_similarity);

            // 5. Normalize the angle from [0, pi] range to [0, 1] range
            double AE = angle_radians / M_PI;

            // 6. Accumulate
            sum_AE += AE;

            // ===================================================
            // METRIC2: RMSE (only module) [0,inf]
            // ===================================================
            double diff_x = gt_x - est_x;
            double diff_y = gt_y - est_y;
            double squared_error = (diff_x * diff_x + diff_y * diff_y);
            sum_squared_error += squared_error;

            // Calculate the magnitude of the ground truth vector: ||GT_i|| for normalization
            sum_gt_magnitudes += std::sqrt(gt_x * gt_x + gt_y * gt_y);
            

            // ===================================================
            // METRIC3: NLPD - Negative Log-Predictive Density (uncertainty aware) [-inf,inf]
            // ===================================================
            double mahalanobis_distance_squared = (pow(gt_x - est_x, 2) / (est_stdX * est_stdX + 1e-6)) +
                                                  (pow(gt_y - est_y, 2) / (est_stdY * est_stdY + 1e-6));

            double NLPD = 0.5 * mahalanobis_distance_squared + std::log(est_stdX) + std::log(est_stdY) + std::log(2 * M_PI);
            
            // Accumulate
            sum_NLPD += NLPD;

            // ===============================
            // VISUALIZE METRIC
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
    }
    metric_pub->publish(metric_marker);


    // RMSE & NRMSE
    double RMSE = std::sqrt(sum_squared_error / N);
    double NRMSE = RMSE / ( sum_gt_magnitudes / N);
    RCLCPP_INFO(this->get_logger(), "[gmrf-validation] RMSE (m/s)= %.2f", RMSE);

    // Average Angular Error
    double AAE = sum_AE / N;
    RCLCPP_INFO(this->get_logger(), "[gmrf-validation] Average Angular Error (rad) = %.2f", AAE);

    // Average NLPD
    double ANLPD = sum_NLPD / N;
    RCLCPP_INFO(this->get_logger(), "[gmrf-validation] Average NLPD = %.2f", ANLPD);

    return {AAE, RMSE, ANLPD};
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
             << "gmrf_std_x,"
             << "gmrf_std_y,"
             << "gt_wind_x,"
             << "gt_wind_y\n";

    for (size_t i = 0; i < N; ++i)
    {
        // Skip occupied cells
        if (gmrf_map->is_cell_free(i))
        {
            // Get GMRF estimation at cell i
            WindVector est_wind = gmrf_map->getEstimation(i);
            double est_x = est_wind.module * cos(est_wind.direction);
            double est_y = est_wind.module * sin(est_wind.direction);
            double est_stdX = est_wind.stdDevX;
            double est_stdY = est_wind.stdDevY;

            // Get GT wind at cell i
            WindVectorXY gt_wind = gt_map[i];
            double gt_x = gt_wind.x;
            double gt_y = gt_wind.y;

            // Save to CSV
            csv_file << i << ","
                     << est_x << ","
                     << est_y << ","
                     << est_stdX << ","
                     << est_stdY << ","
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


void Cvalgt::update_lambdas(double lambda_reg, double lambda_flux, double lambda_obstacles, double lambda_obs)
{
    gmrf_map->update_lambdas(lambda_reg, lambda_flux, lambda_obstacles, lambda_obs);
}

void Cvalgt::read_lambdas(double &lambda_reg, double &lambda_flux, double &lambda_obstacles, double &lambda_obs)    
{
    gmrf_map->read_lambdas(lambda_reg, lambda_flux, lambda_obstacles, lambda_obs);
}



void saveMetricsToCSV(const std::string& filename, 
               const std::vector<int>& v_n_obs, 
               const std::vector<double>& v_AAE,
               const std::vector<double>& v_RMSE,
               const std::vector<double>& v_NLPD,
               const std::vector<double>& l_reg, 
               const std::vector<double>& l_flux, 
               const std::vector<double>& l_obst,
               const std::vector<double>& l_obs)
{
    std::ofstream file(filename);

    if (file.is_open()) {
        // 1. Escribir el encabezado (Header)
        file << "N_Obs, AAE, RMSE, NLPD, Lambda_Reg, Lambda_Flux, Lambda_Obstacles, Lambda_Observations\n";

        // 2. Escribir los datos fila por fila
        for (size_t j = 0; j < v_n_obs.size(); ++j) {
            file << v_n_obs[j] << "," 
                 << v_AAE[j] << "," 
                 << v_RMSE[j] << "," 
                 << v_NLPD[j] << "," 
                 << l_reg[j] << ","
                 << l_flux[j] << "," 
                 << l_obst[j] << ","
                 << l_obs[j] << "\n";
        }

        file.close();
        std::cout << "Datos guardados exitosamente en " << filename << std::endl;
    } else {
        std::cerr << "Error: No se pudo abrir el archivo para escribir." << std::endl;
    }
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

    // Initial values for the 4 Lambda parameters (Reg, Flux, Obstacles, Observations)
    double parameters[4] = {0.1, 0.1, 0.1, 0.1};

    // CERES OPTIMIZER
    ceres::Problem problem;

    // <CostFunctor, Number of Residuals, Number of Parameters>
    ceres::CostFunction* cost_function =
        new ceres::NumericDiffCostFunction<CostFunctor, ceres::CENTRAL, 1, 4>(
            new CostFunctor(my_gmrf_map.get())
        );

    problem.AddResidualBlock(cost_function, nullptr, parameters);

    // 3. Configure and Run the solver
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;


    // LOOP
    //=========================================
    // Vectores para almacenar los resultados
    std::vector<int> v_n_obs;
    std::vector<double> v_AAE;
    std::vector<double> v_RMSE;
    std::vector<double> v_NLPD;
    std::vector<double> l_reg;
    std::vector<double> l_flux;
    std::vector<double> l_obst;
    std::vector<double> l_obs;
    
    for (int i = 1; i <= 500; i += 10) 
    {
        // Simulate N (random) observations
        my_gmrf_map->SimulateWindObservations(i);

        // Optimize Lambdas
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);

        // Output results
        std::cerr << summary.FullReport() << "\n";
        std::cerr << "Final parameters: "
                    << parameters[0] << ", " << parameters[1] << ", "
                    << parameters[2] << ", " << parameters[3] << "\n";

        // Compute performance metrics with current parameters
        std::array<double, 3> metrics = my_gmrf_map->compute_performance_metrics();
        double AAE_opt = metrics[0];  // AAE is the first element
        double RMSE_opt = metrics[1]; // RMSE is the second element
        double NLPD_opt = metrics[2]; // NLPD is the third element

        // Guardar resultados
        v_n_obs.push_back(i);
        v_AAE.push_back(AAE_opt);
        v_RMSE.push_back(RMSE_opt);
        v_NLPD.push_back(NLPD_opt);
        l_reg.push_back(parameters[0]);
        l_flux.push_back(parameters[1]);
        l_obst.push_back(parameters[2]);
        l_obs.push_back(parameters[3]);

        // Update Estimation + Uncertainty
        //my_gmrf_map->update();
        
        // Publish Map as markers (RVIZ2)
        //my_gmrf_map->publishMaps();

        // Save estimation to file (for visualization later)
        std::string file_name = "gmrf_opt_estimation_obs_" + std::to_string(i) + ".csv";
        my_gmrf_map->saveGMRFEstimationToCSV(file_name);
    }

    std::string nombreArchivo = "gmrf_nlpd_optimization_n_obs3.csv";
    saveMetricsToCSV(nombreArchivo, v_n_obs, v_AAE, v_RMSE, v_NLPD, l_reg, l_flux, l_obst, l_obs);
    return 0;
}
