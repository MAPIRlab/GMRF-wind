#pragma once

#include <geometry_msgs/msg/quaternion.hpp>
#include <tf2_geometry_msgs/tf2_geometry_msgs.hpp>

namespace Utils
{

    static double getYaw(const geometry_msgs::msg::Quaternion& quat)
    {
        tf2::Quaternion tfquat;
        tf2::fromMsg(quat, tfquat);

        tf2::Matrix3x3 m(tfquat);
        double roll, pitch, yaw;
        m.getRPY(roll, pitch, yaw);
        return yaw;
    }

    static geometry_msgs::msg::Quaternion createQuaternionMsgFromYaw(double yaw)
    {
        return tf2::toMsg(tf2::Quaternion(tf2::Vector3(0, 0, 1), yaw));
    }

}