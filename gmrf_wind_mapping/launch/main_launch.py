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

    map_server = Node(
        package='nav2_map_server',
        executable='map_server',
        name='map_server',
        output='screen',
        parameters=[
            {'use_sim_time': use_sim_time},
            {'yaml_filename' : os.path.join(pkg_dir, 'launch', 'demo_map.yaml')}
        ]
    )

    gmrf_wind = Node(
        package="gmrf_wind_mapping",
        executable="gmrf_wind_mapping_node",
        name="gmrf",
        parameters=[
            {"frame_id": "map"},                              # Frame where to plot the map, usually (map)
            {"sensor_topic": "/anemometer"},                  # Topic where the anemometer measurements are published
            {"map_file": ""},                                 # Path to a pre-recorded Occupancy GridMap file (grayscale Image). If empty, will listen to map_topic
            {"map_topic": "map"},                             # Topic where the Occupancy GridMap is published
            {"exec_freq": 10.0},                              # Frequency (Hz) to execute the GMRF update step
            {"cell_size": 0.25},                              # Size of each cell in the GMRF grid (meters)
            {"verbose": False},                               # Verbose mode for debugging
            {"visualize_gmrf": True},                         # Visualize the GMRF wind field in RViz
            # Lambda weights for the different priors & Obs in the GMRF model
            {"GMRF_lambdaPrior_reg": 1.0},                      # Regularization -> neighbour cells have similar wind vectors
            {"GMRF_lambdaPrior_mass_conservation": 1000.0},    # Mass conservation law -> divergence of the wind field is zero
            {"GMRF_lambdaPrior_obstacles": 10.0},               # Obstacles --> cells close to obstacles has only tangencial wind
            {"GMRF_lambdaObs": 10000.0},                         # The initial weight (Lambda) of each observation
            {"GMRF_lambdaObsLoss": 0.0}                       # The loss of information (Lambda) of the observations with each iteration
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

    # LIFECYCLE MANAGER
    lifecicle = Node(
        package='nav2_lifecycle_manager',
        executable='lifecycle_manager',
        name='lifecycle_manager_localization',
        output='screen',
        parameters=[
            {'use_sim_time': use_sim_time},
            {'autostart': True},
            {'node_names': ['map_server']}
        ]
    )
    
    actions = []    
    actions.append(rviz)
    actions.append(gmrf_wind)    
    actions.append(map_server)
    actions.append(lifecicle)
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