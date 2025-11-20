import os
from pathlib import Path
from ament_index_python.packages import get_package_share_directory
from launch import LaunchDescription
from launch.conditions import IfCondition
from launch.conditions import UnlessCondition
from launch.actions import DeclareLaunchArgument
from launch.actions import IncludeLaunchDescription
from launch.launch_description_sources import PythonLaunchDescriptionSource
from launch.substitutions import ThisLaunchFileDir
from launch.actions import ExecuteProcess
from launch.substitutions import LaunchConfiguration, PythonExpression
from launch_ros.actions import Node

def generate_launch_description():

    return LaunchDescription([

            Node(
                package='gmrf_wind_mapping',
                executable='gmrf_wind_mapping_node',
                name='gmrf_wind_mapping_node',
                output='screen',
                parameters=[
                    {"frame_id": "map"},                              # Frame where to plot the map, usually (map)
                    {"sensor_topic": "/anemometer"},                  # Topic where the anemometer measurements are published
                    {"map_file": ""},                                 # Path to a pre-recorded Occupancy GridMap file (grayscale Image). If empty, will listen to map_topic
                    {"map_topic": "map"},                             # Topic where the Occupancy GridMap is published
                    {"exec_freq": 10.0},                              # Frequency (Hz) to execute the GMRF update step
                    {"cell_size": 0.25},                              # Size of each cell in the GMRF grid (meters)
                    {"verbose": False},                               # Verbose mode for debugging
                    # Lambda weights for the different priors & Obs in the GMRF model
                    {"GMRF_lambdaPrior_reg": 1.0},                      # Regularization -> neighbour cells have similar wind vectors
                    {"GMRF_lambdaPrior_mass_conservation": 10000.0},    # Mass conservation law -> divergence of the wind field is zero
                    {"GMRF_lambdaPrior_obstacles": 10.0},               # Obstacles --> cells close to obstacles has only tangencial wind
                    {"GMRF_lambdaObs": 10.0},                         # The initial weight (Lambda) of each observation
                    {"GMRF_lambdaObsLoss": 0.0}                       # The loss of information (Lambda) of the observations with each iteration
                ]
            ),
    ])