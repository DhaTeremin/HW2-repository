#!/bin/bash
set -e  # stop on error

# Update packages
sudo apt update

# Install ROS dependencies
rosdep install --ignore-src --from-paths . -y -r

# Build the workspace
colcon build

# Source the workspace
source install/setup.bash

echo "Build and setup completed."
