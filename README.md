# HW2: Control your Robot

To install any missing dependencies:

````bash
./src/quick_build.sh
````

Then source the overlay in every terminal needed for the behaviour demontration

````bash
source install/setup.bash
````


## Kinematic Control

To start the RViz2 visualization of the manipulator:

````bash
ros2 launch iiwa_bringup iiwa.launch.py 
````
This command will bring up the robot in RViz2 and load the position command interface and position controller. The parameters `command_interface` and `robot_controller` can be set in order to choose which interface to choose between velocity and effort.

Velocity interface and controller:
````bash
ros2 launch iiwa_bringup iiwa.launch.py command_interface:="velocity" robot_controller:="velocity_controller"
````

Effort interface and controller:
````bash
ros2 launch iiwa_bringup iiwa.launch.py command_interface:="effort" robot_controller:="effort_controller"
````

In another terminal, the ros2_kdl_package can be started via launch file, loading the arguments specified in the args.yaml file (in /ros2_kdl_pacakge/config) via the command:

````bash
 ros2 launch ros2_kdl_package ros2_kdl.launch.py
````

#### Veclotity Control in the Null Space

In order to start the newly implemented velocity controller in the null space, the user can either set the `ctrl` parameter in the args.yaml to `velocity_ctrl_null` or start the node via:

````bash
 ros2 run ros2_kdl_package ros2_kdl_node --ros-args -p cmd_interface:="velocity" -p ctrl:="velocity_ctrl_null" -p end_position:="[0.2, 0.0, 0.3]"
````

Any other parameter, such as the gain `Kp` and `lambda`, can be further specified in the command line or the yaml file.

#### Action server mode

The action server mode can be set by modifying the `action_server_mode` parameter in the yaml file to true or via command line:

````bash
 ros2 run ros2_kdl_package ros2_kdl_node --ros-args -p action_server_mode:="true" -p end_position:="[0.2, 0.0, 0.3]"
````
The node will wait for an action client request. The action client can be started in another terminal via:

````bash
ros2 run ros2_kdl_package ros2_kdl_action_client 
````

The position and orientation error will be printed to screen.

## Vision Based Control

The gazebo world containing the aruco tag and the iiwa manipulator can be started by exporting the environment variable:

````bash
export GZ_SIM_RESOURCE_PATH=$GZ_SIM_RESOURCE_PATH:~/ros2_ws/src/ros2_iiwa/iiwa_description/gazebo/models/arucotag
````
To start just the gazebo world containing the aruco tag:
````bash
ros2 launch ros_gz_sim gz_sim.launch.py gz_args:="~/ros2_ws/src/ros2_iiwa/iiwa_description/gazebo/worlds/empty.world"
````
 Otherwise to start the simulation comprising the robot and the controllers(in this case with velocity interface and controller):

````bash
ros2 launch iiwa_bringup iiwa.launch.py use_sim:=true gz_args:="-r ~/ros2_ws/src/ros2_iiwa/iiwa_description/gazebo/worlds/empty.world" command_interface:="velocity" robot_controller:="velocity_controller"
````

In another terminal the aruco single node can be started to recognize the tag:

````bash
ros2 run aruco_ros single --ros-args \
  -p marker_id:=201 \
  -p marker_size:=0.1 \
  -p camera_frame:=link_7 \
  -p marker_frame:=aruco_marker_frame \
  -p reference_frame:=link_7 \
  -r /image:=/videocamera \
  -r /camera_info:=/camera_info
````

Via `ros2 topic echo` command the `PoseStamped` message with the transformation can be read:

````bash
ros2 topic echo /aruco_single/pose --once
````
Then to start the vision based control (Note: this control has not the expected behaviour):
````bash
 ros2 run ros2_kdl_package ros2_kdl_node --ros-args -p cmd_interface:="velocity" -p ctrl:="vision_ctrl"
````
To send a service call, addressing the point 2.c, to modify the marker position the command is:

````bash
ros2 service call /world/default/set_pose ros_gz_interfaces/srv/SetEntityPose "{
  entity: { name: 'arucotag', id: 0, type: 2 },
  pose: { position: {x: 0.7, y: 0.2, z: 2.1}, orientation: {w: 1.0, x: 0.0, y: 0.0, z: 0.0} }
}"
````
(Note to send a service call it is mandatory to run the ros2_kdl.launch.py, because the node defining the bridge is created inside it)
