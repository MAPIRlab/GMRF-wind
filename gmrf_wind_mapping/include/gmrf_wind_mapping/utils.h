#pragma once

#include <filesystem>
#include <geometry_msgs/msg/quaternion.hpp>
#include <nav_msgs/msg/occupancy_grid.hpp>
#include <opencv2/highgui.hpp>
#include <tf2_geometry_msgs/tf2_geometry_msgs.hpp>
#include <visualization_msgs/msg/marker_array.hpp>
#include <visualization_msgs/msg/marker.hpp>
#include <rclcpp/rclcpp.hpp>

namespace Utils
{

    // Get yaw from a quaternion
    inline double getYaw(const geometry_msgs::msg::Quaternion& quat)
    {
        tf2::Quaternion tfquat;
        tf2::fromMsg(quat, tfquat);

        tf2::Matrix3x3 m(tfquat);
        double roll, pitch, yaw;
        m.getRPY(roll, pitch, yaw);
        return yaw;
    }


    // Create a quaternion message from a yaw angle
    inline geometry_msgs::msg::Quaternion createQuaternionMsgFromYaw(double yaw)
    {
        return tf2::toMsg(tf2::Quaternion(tf2::Vector3(0, 0, 1), yaw));
    }

    
    // Parse a grayscale image as an occupancy gridmap (To be improved)
    inline nav_msgs::msg::OccupancyGrid parseMapImage(const std::string& path)
    {
        if (!std::filesystem::exists(path))
        {
            RCLCPP_ERROR(rclcpp::get_logger("GMRF"), "Tried to parse map image at path %s, but it does not exist", path.c_str());
            raise(SIGTRAP);
        }

        cv::Mat mapImage = cv::imread(path, cv::IMREAD_GRAYSCALE);
        cv::flip(mapImage, mapImage, 0);
        size_t width = mapImage.size().width;
        size_t height = mapImage.size().height;

        nav_msgs::msg::OccupancyGrid occupancyGrid;
        occupancyGrid.data.resize(width * height);
        for (int i = 0; i < width * height; i++)
            occupancyGrid.data[i] = (int8_t)(100 - std::clamp((int)mapImage.data[i], 0, 100));

        return occupancyGrid;
    }

    void get_arrow_color(double module, double max_module, float& color_r, float& color_g, float& color_b)
    {
        // Predefined colormap (jet-like) with 200 entries
        // Declare the arrays as static. Initialization only happens once.
        // Use const to ensure they are read-only after initialization.
        static const float temp_color_r[200] = {
            0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
            0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
            0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
            0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34,
            0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58, 0.60, 0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.74, 0.76, 0.78, 0.80,
            0.82, 0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
            1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
            1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.98, 0.96, 0.94, 0.92, 0.90, 0.88, 0.86, 0.84, 0.82,
            0.80, 0.78, 0.76, 0.74, 0.72, 0.70, 0.68, 0.66, 0.64, 0.62, 0.60, 0.58, 0.56, 0.54, 0.52, 0.50};
        static const float temp_color_g[200] = {
            0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
            0.00, 0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.42,
            0.44, 0.46, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58, 0.60, 0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.74, 0.76, 0.78, 0.80, 0.82, 0.84, 0.86, 0.88,
            0.90, 0.92, 0.94, 0.96, 0.98, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
            1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
            1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.98, 0.96, 0.94, 0.92, 0.90, 0.88, 0.86, 0.84, 0.82, 0.80, 0.78, 0.76, 0.74,
            0.72, 0.70, 0.68, 0.66, 0.64, 0.62, 0.60, 0.58, 0.56, 0.54, 0.52, 0.50, 0.48, 0.46, 0.44, 0.42, 0.40, 0.38, 0.36, 0.34, 0.32, 0.30, 0.28,
            0.26, 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.08, 0.06, 0.04, 0.02, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
            0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
        static const float temp_color_b[200] = {
            0.52, 0.54, 0.56, 0.58, 0.60, 0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.74, 0.76, 0.78, 0.80, 0.82, 0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96,
            0.98, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
            1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
            1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.98, 0.96, 0.94, 0.92, 0.90, 0.88, 0.86, 0.84, 0.82, 0.80, 0.78, 0.76, 0.74, 0.72, 0.70, 0.68, 0.66,
            0.64, 0.62, 0.60, 0.58, 0.56, 0.54, 0.52, 0.50, 0.48, 0.46, 0.44, 0.42, 0.40, 0.38, 0.36, 0.34, 0.32, 0.30, 0.28, 0.26, 0.24, 0.22, 0.20,
            0.18, 0.16, 0.14, 0.12, 0.10, 0.08, 0.06, 0.04, 0.02, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
            0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
            0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
            0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};

        // Get color index
        size_t idx_color = 199 * (module / max_module);
        color_r = temp_color_r[idx_color];
        color_g = temp_color_g[idx_color];
        color_b = temp_color_b[idx_color];
    }

    void createWindMarkerArrayFromGMRF(CGMRF_map& my_map, 
                                        std::string frame_id, 
                                        visualization_msgs::msg::MarkerArray& wind_array, 
                                        visualization_msgs::msg::Marker& wind_std_array)
    {
        using visualization_msgs::msg::Marker;
        using visualization_msgs::msg::MarkerArray;

        // Clear previous markers
        wind_array.markers.clear();
        wind_std_array.points.clear();

        // Get GMRF dimensions
        Eigen::Vector2i dims = my_map.map_size();
        int width = dims.x();
        int height = dims.y();
        if (width <= 0 || height <= 0)
            return;
        // Number of cells
        size_t N = static_cast<size_t>(width) * static_cast<size_t>(height);

        // Get cell size
        double cell_size = my_map.map_resolution();
        if (cell_size <= 0.0)
            return;

        // Get max wind vector in the map (to normalize the plot)
        double max_module = 0.0;
        double max_var = 0.0;
        for (size_t i = 0; i < N; ++i)
        {
            WindVector w = my_map.getEstimation(static_cast<int>(i));
            // module 
            if (w.module > max_module)
                max_module = w.module;
            // Uncertainty (det of Sigma). Can also consider the Trace of Sigma or the max eigenvalue
            double detSigma = sqrt(w.varX * w.varY - w.covXY * w.covXY);
            if (detSigma > max_var)
                max_var = detSigma;
        }

        // Uncertainty is a single marker of type POINTS (heatMap of uncertainty)
        wind_std_array.header.frame_id = frame_id;
        wind_std_array.header.stamp = rclcpp::Time(0);
        wind_std_array.ns = "gmrf_wind_stddev";
        wind_std_array.id = 0;
        wind_std_array.type = Marker::POINTS;
        wind_std_array.action = Marker::ADD;
        // shape
        wind_std_array.scale.x = cell_size;
        wind_std_array.scale.y = cell_size;

        // And one ARROW/POINTS marker per cell for wind_vector/var
        wind_array.markers.reserve(N);
        wind_std_array.points.reserve(N);
        wind_std_array.colors.reserve(N);
        for (size_t i = 0; i < N; ++i)
        {
            // Get wind estimation at cell i
            WindVector w = my_map.getEstimation(static_cast<int>(i));
            Eigen::Vector2d vec = w.asEigen();                              //vec.x(); vec.y()

            // Cell center coordinates
            double cx = 0.0, cy = 0.0;
            my_map.id2xy_public(i, cx, cy);

            // Point marker for Uncertainty
            {
                // Add one point at the cell center
                geometry_msgs::msg::Point p;
                p.x = cx;
                p.y = cy;
                p.z = 0.0;
                wind_std_array.points.push_back(p);
                
                std_msgs::msg::ColorRGBA color;
                // color -> must normalize to [0-199]
                // Uncertainty (det of Sigma). Can also consider the Trace of Sigma or the max eigenvalue
                double detSigma = sqrt(w.varX * w.varY - w.covXY * w.covXY);
                get_arrow_color(detSigma, max_var, color.r, color.g, color.b);
                color.a = 1.0;       // transparency
                wind_std_array.colors.push_back(color);
            }

            // Arrow marker for wind vector
            if (w.module < 0.01 || max_module < 0.01)
                continue;   // skip near-zero vectors
            
            // Create Arrow marker
            Marker m;
            m.header.frame_id = frame_id;
            m.header.stamp = rclcpp::Time(0);
            m.ns = "gmrf_wind";
            m.id = static_cast<int>(i);
            m.type = Marker::ARROW;
            m.action = Marker::ADD;

            // Set pose
            m.pose.position.x = cx;
            m.pose.position.y = cy;
            m.pose.orientation = createQuaternionMsgFromYaw(atan2(vec.y(), vec.x()));
            // shape
            m.scale.x = cell_size; // arrow length,
            m.scale.y = 0.03;      // arrow width
            m.scale.z = 0.05;      // arrow height
            // color -> must normalize to [0-199]
            get_arrow_color(w.module, max_module, m.color.r, m.color.g, m.color.b);
            m.color.a = 1.0;       // transparency

            // Add marker to Array
            wind_array.markers.push_back(m);
        }
    }
} // namespace Utils