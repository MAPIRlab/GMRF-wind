#pragma once

#include <filesystem>
#include <geometry_msgs/msg/quaternion.hpp>
#include <nav_msgs/msg/occupancy_grid.hpp>
#include <opencv2/highgui.hpp>
#include <tf2_geometry_msgs/tf2_geometry_msgs.hpp>

namespace Utils
{

    inline double getYaw(const geometry_msgs::msg::Quaternion& quat)
    {
        tf2::Quaternion tfquat;
        tf2::fromMsg(quat, tfquat);

        tf2::Matrix3x3 m(tfquat);
        double roll, pitch, yaw;
        m.getRPY(roll, pitch, yaw);
        return yaw;
    }

    inline geometry_msgs::msg::Quaternion createQuaternionMsgFromYaw(double yaw)
    {
        return tf2::toMsg(tf2::Quaternion(tf2::Vector3(0, 0, 1), yaw));
    }

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
} // namespace Utils