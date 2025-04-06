#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <execution>
#include <fstream>
#include <mutex>
#include <pcl/common/io.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <ros/ros.h>
#include <sstream>
#include <stdio.h>
#include <string>
#include <unordered_map>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>

#define HASH_P 116101
#define MAX_N 10000000000

typedef struct ConfigSetting {
  /* for point cloud pre-preocess*/
  int useful_corner_num_ = 30;

  /* for key points*/
  float plane_merge_normal_thre_;
  float plane_merge_dis_thre_;
  float plane_detection_thre_ = 0.01;
  float voxel_size_ = 1.0;
  int voxel_init_num_ = 10;
  int proj_plane_num_ = 3;
  float proj_image_resolution_ = 0.5;
  float proj_image_high_inc_ = 0.5;
  float proj_dis_min_ = 0;
  float proj_dis_max_ = 5;
  float summary_min_thre_ = 10;
  int line_filter_enable_ = 0;
  int touch_filter_enable_ = 0;

  /* for STD */
  float descriptor_near_num_ = 10;
  float descriptor_min_len_ = 1;
  float descriptor_max_len_ = 10;
  float non_max_suppression_radius_ = 3.0;
  float std_side_resolution_ = 0.2;

  /* for place recognition*/
  int skip_near_num_ = 20;
  int candidate_num_ = 50;
  float rough_dis_threshold_ = 0.03;
  float similarity_threshold_ = 0.7;
  float icp_threshold_ = 0.5;
  float normal_threshold_ = 0.1;
  float dis_threshold_ = 0.3;

} ConfigSetting;

typedef struct BinaryDescriptor {
  std::vector<bool> occupy_array_;
  unsigned char summary_;
  Eigen::Vector3d location_;
} BinaryDescriptor;

typedef struct BinaryDescriptorF {
  // std::vector<bool> occupy_array_;
  // bool occupy_array_[49];
  unsigned char summary_;
  Eigen::Vector3f location_;
} BinaryDescriptorF;

// 1kb,12.8
typedef struct STD {
  Eigen::Vector3d triangle_;
  Eigen::Vector3d angle_;
  Eigen::Vector3d center_;
  // unsigned short frame_number_;
  int frame_number_;
  // std::vector<unsigned short> score_frame_;
  // std::vector<Eigen::Matrix3d> position_list_;
  BinaryDescriptor binary_A_;
  BinaryDescriptor binary_B_;
  BinaryDescriptor binary_C_;
} STD;

typedef struct BTCPlane {
  pcl::PointXYZINormal p_center_;
  Eigen::Vector3d center_;
  Eigen::Vector3d normal_;
  Eigen::Matrix3d covariance_;
  float radius_ = 0;
  float min_eigen_value_ = 1;
  float d_ = 0;
  int id_ = 0;
  int sub_plane_num_ = 0;
  int points_size_ = 0;
  bool is_plane_ = false;
} BTCPlane;

typedef struct STDMatchList {
  std::vector<std::pair<STD, STD>> match_list_;
  std::pair<int, int> match_id_;
  int match_frame_;
  double mean_dis_;
} STDMatchList;

class BTCVOXEL_LOC {
public:
  int64_t x, y, z;

  BTCVOXEL_LOC(int64_t vx = 0, int64_t vy = 0, int64_t vz = 0)
      : x(vx), y(vy), z(vz) {}

  bool operator==(const BTCVOXEL_LOC &other) const {
    return (x == other.x && y == other.y && z == other.z);
  }
};

// Hash value
namespace std {
template <> struct hash<BTCVOXEL_LOC> {
  int64_t operator()(const BTCVOXEL_LOC &s) const {
    using std::hash;
    using std::size_t;
    return ((((s.z) * HASH_P) % MAX_N + (s.y)) * HASH_P) % MAX_N + (s.x);
  }
};
} // namespace std

class STD_LOC {
public:
  int64_t x, y, z, a, b, c;

  STD_LOC(int64_t vx = 0, int64_t vy = 0, int64_t vz = 0, int64_t va = 0,
          int64_t vb = 0, int64_t vc = 0)
      : x(vx), y(vy), z(vz), a(va), b(vb), c(vc) {}

  bool operator==(const STD_LOC &other) const {
    return (x == other.x && y == other.y && z == other.z);
    // return (x == other.x && y == other.y && z == other.z && a == other.a &&
    //         b == other.b && c == other.c);
  }
};

namespace std {
template <> struct hash<STD_LOC> {
  int64_t operator()(const STD_LOC &s) const {
    using std::hash;
    using std::size_t;
    return ((((s.z) * HASH_P) % MAX_N + (s.y)) * HASH_P) % MAX_N + (s.x);
  }
};
} // namespace std

class BTCOctoTree {
public:
  ConfigSetting config_setting_;
  std::vector<Eigen::Vector3d> voxel_points_;
  BTCPlane *plane_ptr_;
  int layer_;
  int octo_state_; // 0 is end of tree, 1 is not
  int merge_num_ = 0;
  bool is_project_ = false;
  std::vector<Eigen::Vector3d> project_normal;
  bool is_publish_ = false;
  BTCOctoTree *leaves_[8];
  double voxel_center_[3]; // x, y, z
  float quater_length_;
  bool init_octo_;

  // for plot
  bool is_check_connect_[6];
  bool connect_[6];
  BTCOctoTree *connect_tree_[6];

  BTCOctoTree(const ConfigSetting &config_setting)
      : config_setting_(config_setting) {
    voxel_points_.clear();
    octo_state_ = 0;
    layer_ = 0;
    init_octo_ = false;
    for (int i = 0; i < 8; i++) {
      leaves_[i] = nullptr;
    }
    // for plot
    for (int i = 0; i < 6; i++) {
      is_check_connect_[i] = false;
      connect_[i] = false;
      connect_tree_[i] = nullptr;
    }
    plane_ptr_ = new BTCPlane;
  }
  void init_plane();
  void init_octo_tree();

  ~BTCOctoTree()
  {
    delete plane_ptr_;
  }

};

void load_config_setting(std::string &config_file,
                         ConfigSetting &config_setting);

double binary_similarity(const BinaryDescriptor &b1,
                         const BinaryDescriptor &b2);

bool binary_greater_sort(BinaryDescriptor a, BinaryDescriptor b);
bool plane_greater_sort(BTCPlane *plane1, BTCPlane *plane2);

// double
// calc_triangle_dis(const std::vector<std::pair<STD, STD>> &match_std_list);

void read_parameters(ros::NodeHandle &nh, ConfigSetting &config_setting, int isHighFly);

Eigen::Vector3d normal2vec(const pcl::PointXYZINormal &pi);

template <typename T> Eigen::Vector3d point2vec(const T &pi) {
  Eigen::Vector3d vec(pi.x, pi.y, pi.z);
  return vec;
}

double time_inc(std::chrono::_V2::system_clock::time_point &t_end,
                std::chrono::_V2::system_clock::time_point &t_begin);


class STDescManager {
public:
  STDescManager() = default;

  ConfigSetting config_setting_;

  unsigned int current_frame_id_;

  STDescManager(ConfigSetting &config_setting)
      : config_setting_(config_setting) {
    current_frame_id_ = 0;
  };

  // std::vector<BinaryDescriptor> vec_binary;

  // hash table, save all descriptors
  std::unordered_map<STD_LOC, std::vector<STD>> data_base_;

  // save all key clouds, optional
  // std::vector<pcl::PointCloud<pcl::PointXYZI>::Ptr> key_cloud_vec_;

  // save all binary descriptors of key frame
  // std::vector<std::vector<BinaryDescriptor>> history_binary_list_;

  // save all planes of key frame, required
  std::vector<pcl::PointCloud<pcl::PointXYZINormal>::Ptr> plane_cloud_vec_;

  /*Three main processing functions*/

  // generate STDescs from a point cloud
  void GenerateSTDescs(pcl::PointCloud<pcl::PointXYZI>::Ptr &input_cloud,
                       std::vector<STD> &stds_vec, int id);

  // search result <candidate_id, plane icp score>. -1 for no loop
  void SearchLoop(std::vector<STD> &stds_vec,
                  std::pair<int, double> &loop_result,
                  std::pair<Eigen::Vector3d, Eigen::Matrix3d> &loop_transform,
                  std::vector<std::pair<STD, STD>> &loop_std_pair, pcl::PointCloud<pcl::PointXYZINormal>::Ptr pl_cur);

  // add descriptors to database
  void AddSTDescs(const std::vector<STD> &stds_vec);

  // Geometrical optimization by plane-to-plane ico
  void PlaneGeomrtricIcp(
      const pcl::PointCloud<pcl::PointXYZINormal>::Ptr &source_cloud,
      const pcl::PointCloud<pcl::PointXYZINormal>::Ptr &target_cloud,
      std::pair<Eigen::Vector3d, Eigen::Matrix3d> &transform);

private:
  /*Following are sub-processing functions*/

  // voxelization and plane detection
  void init_voxel_map(const pcl::PointCloud<pcl::PointXYZI>::Ptr &input_cloud,
                      std::unordered_map<BTCVOXEL_LOC, BTCOctoTree *> &voxel_map);

  // acquire planes from voxel_map
  void get_plane(const std::unordered_map<BTCVOXEL_LOC, BTCOctoTree *> &voxel_map,
                 pcl::PointCloud<pcl::PointXYZINormal>::Ptr &plane_cloud);

  void get_project_plane(std::unordered_map<BTCVOXEL_LOC, BTCOctoTree *> &feat_map,
                         std::vector<BTCPlane *> &project_plane_list);

  void merge_plane(std::vector<BTCPlane *> &origin_list,
                   std::vector<BTCPlane *> &merge_plane_list);

  // extract corner points from pre-build voxel map and clouds

  void binary_extractor(const std::vector<BTCPlane *> proj_plane_list,
                        const pcl::PointCloud<pcl::PointXYZI>::Ptr &input_cloud,
                        std::vector<BinaryDescriptor> &binary_descriptor_list);

  void extract_binary(const Eigen::Vector3d &project_center,
                      const Eigen::Vector3d &project_normal,
                      const pcl::PointCloud<pcl::PointXYZI>::Ptr &input_cloud,
                      std::vector<BinaryDescriptor> &binary_list);

  // non maximum suppression, to control the number of corners
  void non_maxi_suppression(std::vector<BinaryDescriptor> &binary_list);

  // build STDescs from corner points.
  void generate_std(const std::vector<BinaryDescriptor> &binary_list,
                    const int &frame_id, std::vector<STD> &std_list);

  // Select a specified number of candidate frames according to the number of
  // STDesc rough matches
  void candidate_selector(std::vector<STD> &stds_vec,
                          std::vector<STDMatchList> &candidate_matcher_vec);

  // Get the best candidate frame by geometry check
  void
  candidate_verify(STDMatchList &candidate_matcher, double &verify_score,
                   std::pair<Eigen::Vector3d, Eigen::Matrix3d> &relative_pose,
                   std::vector<std::pair<STD, STD>> &sucess_match_vec, pcl::PointCloud<pcl::PointXYZINormal>::Ptr pl_cur);

  // Get the transform between a matched std pair
  void triangle_solver(std::pair<STD, STD> &std_pair, Eigen::Vector3d &t,
                       Eigen::Matrix3d &rot);

  // Geometrical verification by plane-to-plane icp threshold
  double plane_geometric_verify(
      const pcl::PointCloud<pcl::PointXYZINormal>::Ptr &source_cloud,
      const pcl::PointCloud<pcl::PointXYZINormal>::Ptr &target_cloud,
      const std::pair<Eigen::Vector3d, Eigen::Matrix3d> &transform);
};