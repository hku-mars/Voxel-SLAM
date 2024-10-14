# Voxel-SLAM: A Complete, Accurate, and Versatile LiDAR-Inertial SLAM System

**Voxel-SLAM** is a complete, accurate, and versatile LiDAR-inertial SLAM system that fully utilizes short-term, mid-term, long-term, and multi-map data associations. It includes five modules: initialization, odometry, local mapping, loop closure, and global mapping. The initialization can provide accurate states and local map in a static or dynamic initial state. The odometry estimates current states and detect potential system divergence. The local mapping refine the states and local map within the sliding window by a LiDAR-inertial BA. The loop closure can detect in multiple sessions. The global mapping refine the global map with an efficient hierarchical global BA. The system overview is:

<div align="center">
    <a href="https://youtu.be/Cg9W01aIUzE" target="_blank">
    <img src="./figure/systemoverview.png" width = 60% >
</div>

### Related Video

The video of **Voxel-SLAM** is available on [YouTube](https://youtu.be/Cg9W01aIUzE).

### Related works

Related paper is available on [**arxiv**](https://arxiv.org/abs/2410.08935).

### Competitions

Voxel-SLAM has been served as a subsystem to participate in [ICRA HILTI 2023 SLAM Challenge](https://hilti-challenge.com/leader-board-2023.html) (**2nd** place on the LiDAR single-session) and [ICCV 2023 SLAM Challenge](https://superodometry.com/iccv23_challenge_LiI) (**1st** place on the LiDAR inertial track).

### Codes

1. Voxel-SLAM: The codes will be published once the paper is accepted.

2. VoxelSLAMPointCloud2: A customized plugin for RViz. It has the same usage to original "PointCloud2" in RViz, but it can clear the point cloud map automatically when receiving an empty point cloud, with any **Decay Time** of the plugin. 

   (1) Put the "VoxelSLAMPointCloud2" within the same "src" folder to your project and make.

   (2) “roslaunch” your program with RViz. 

   (3) Click the "Add" button to add the "VoxelSLAMPointCloud2".

