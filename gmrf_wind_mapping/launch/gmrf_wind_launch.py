import os
from launch import LaunchDescription
from launch.actions import (DeclareLaunchArgument, SetLaunchConfiguration, IncludeLaunchDescription,
                            SetEnvironmentVariable, OpaqueFunction, GroupAction, Shutdown, ExecuteProcess)
from launch.launch_description_sources import FrontendLaunchDescriptionSource, PythonLaunchDescriptionSource
from launch.substitutions import LaunchConfiguration, PathJoinSubstitution, FindExecutable
from launch_ros.actions import Node, PushRosNamespace
from ament_index_python.packages import get_package_share_directory
from launch.frontend.parse_substitution import parse_substitution

# ===========================

def launch_setup(context, *args, **kwargs):
    
    # CONFIG
    use_sim_time = False
    pkg_dir = get_package_share_directory("gmrf_wind_mapping")
    test_env_dir = get_package_share_directory("test_env")

    gmrf_wind = Node(
        package="gmrf_wind_mapping",
        executable="gmrf_wind_mapping_node",
        name="gmrf",
        parameters=[
            {"frame_id": "map"},                              # Frame where to plot the map, usually (map)
            {"sensor_topic": "/anemometer"},                  # Topic where the anemometer measurements are published
            #{"map_yaml_file": os.path.join(pkg_dir, "launch", "demo_map.yaml")},                                 # Path to a pre-recorded Occupancy GridMap file (grayscale Image). If empty, will listen to map_topic
            {"map_yaml_file": os.path.join(test_env_dir, "scenarios", "10x6_empty_room", "occupancy.yaml")},                                 # Path to a pre-recorded Occupancy GridMap file (grayscale Image). If empty, will listen to map_topic
            {"map_topic": "map"},                             # Topic where the Occupancy GridMap is published
            {"exec_freq": 10.0},                              # Frequency (Hz) to execute the GMRF update step
            {"cell_size": 0.3},                              # Size of each cell in the GMRF grid (meters)
            {"verbose": False},                               # Verbose mode for debugging
            {"visualize_gmrf": True},                         # Visualize the GMRF wind field in RViz
            # Lambda weights for the different priors & Obs in the GMRF model
            {"observation_var_wind_speed": 0.00025},           # Variance of the wind speed measurement (m/s)^2
            {"observation_var_wind_direction": 0.00025},       # Variance of the wind direction measurement (rad)^2
            {"GMRF_lambdaPrior_advection": 10.0},              # Advection constraint -> neighboring cells should have similar wind values in the direction of the wind
            {"GMRF_lambdaPrior_mass_conservation": 10.0},      # Mass conservation law -> divergence of the wind field is zero
            {"GMRF_lambdaPrior_diffusion": 1.0},               # Diffusion constraint -> neighboring cells should have similar wind values in all directions
            {"GMRF_lambdaPrior_vorticity": 1.0},               # Vorticity constraint -> curl of the wind field is zero
            {"GMRF_lambdaPrior_obstacles": 1.0},               # Obstacles --> cells close to obstacles has only tangencial wind
        ]
    )

    rviz = Node(
        package = "rviz2",
        executable = "rviz2",
        name = "rviz2",
        # prefix="xterm -e",
        arguments=[
            "-d" + os.path.join(pkg_dir, "launch", "gmrf.rviz")
        ],
    )
    
    
    actions = []    
    actions.append(rviz)
    actions.append(gmrf_wind)
    return actions


def generate_launch_description():

    launch_description = [
        # Set env var to print messages to stdout immediately
        SetEnvironmentVariable("RCUTILS_LOGGING_BUFFERED_STREAM", "1"),
        SetEnvironmentVariable("RCUTILS_COLORIZED_OUTPUT", "1"),
        SetLaunchConfiguration(
            name="pkg_dir",
            value=[get_package_share_directory("gmrf_wind_mapping")],
        ),
    ]
    
    launch_description.append(OpaqueFunction(function=launch_setup))

    return LaunchDescription(launch_description)