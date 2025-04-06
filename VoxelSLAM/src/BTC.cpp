#include "BTC.h"

void read_parameters(ros::NodeHandle &nh, ConfigSetting &config_setting, int isHighFly)
{
  if(!isHighFly)
  {
    config_setting.useful_corner_num_ = 100;
    config_setting.plane_merge_normal_thre_ = 0.1; 
    config_setting.plane_merge_dis_thre_ = 0.3;
    config_setting.plane_detection_thre_ = 0.01;
    config_setting.voxel_size_ = 1; // adjust
    config_setting.voxel_init_num_ = 10;
    config_setting.proj_plane_num_ = 2; // adjust
    config_setting.proj_image_resolution_ = 0.5;
    config_setting.proj_image_high_inc_ = 0.1;
    config_setting.proj_dis_min_ = 0;
    config_setting.proj_dis_max_ = 5;
    config_setting.summary_min_thre_ = 10;
    config_setting.line_filter_enable_ = 1;
    config_setting.touch_filter_enable_ = 0;

    config_setting.descriptor_near_num_ = 15; // adjust
    config_setting.descriptor_min_len_ = 2;
    config_setting.descriptor_max_len_ = 50; // adjust
    config_setting.non_max_suppression_radius_ = 2;
    config_setting.std_side_resolution_ = 0.2;

    config_setting.skip_near_num_ = 30;
    config_setting.candidate_num_ = 20;
    config_setting.rough_dis_threshold_ = 0.01;
    config_setting.similarity_threshold_ = 0.7;
    config_setting.icp_threshold_ = 0.15;
    config_setting.normal_threshold_ = 0.2;
    config_setting.dis_threshold_ = 0.5;
  }
  else
  {
    // avia fly
    config_setting.useful_corner_num_ = 200;
    config_setting.plane_merge_normal_thre_ = 0.3;  // 0.1
    config_setting.plane_merge_dis_thre_ = 0.6;     // 0.3
    config_setting.plane_detection_thre_ = 0.05;    // 0.01
    config_setting.voxel_size_ = 2;                 // adjust
    config_setting.voxel_init_num_ = 10;
    config_setting.proj_plane_num_ = 1;  // adjust
    config_setting.proj_image_resolution_ = 0.5;
    config_setting.proj_image_high_inc_ = 0.2;  // adjust
    config_setting.proj_dis_min_ = 0;
    config_setting.proj_dis_max_ = 10;     // adjust
    config_setting.summary_min_thre_ = 6;  // adjust
    config_setting.line_filter_enable_ = 0;
    config_setting.touch_filter_enable_ = 0;

    config_setting.descriptor_near_num_ = 15;        // adjust
    config_setting.descriptor_min_len_ = 3;          // adjust
    config_setting.descriptor_max_len_ = 50;         // adjust
    config_setting.non_max_suppression_radius_ = 3;  // adjust
    config_setting.std_side_resolution_ = 0.2;

    config_setting.skip_near_num_ = 30;
    config_setting.candidate_num_ = 100;
    config_setting.rough_dis_threshold_ = 0.01;
    config_setting.similarity_threshold_ = 0.5;  // adjust
    config_setting.icp_threshold_ = 0.15;
    config_setting.normal_threshold_ = 0.2;
    config_setting.dis_threshold_ = 0.5;
  }
}

double binary_similarity(const BinaryDescriptor &b1,
                         const BinaryDescriptor &b2) {
  double dis = 0;
  for (size_t i = 0; i < b1.occupy_array_.size(); i++) {
    // to be debug hanming distance
    if (b1.occupy_array_[i] == true && b2.occupy_array_[i] == true) {
      dis += 1;
    }
  }
  return 2 * dis / (b1.summary_ + b2.summary_);
}

bool binary_greater_sort(BinaryDescriptor a, BinaryDescriptor b) {
  return (a.summary_ > b.summary_);
}

bool plane_greater_sort(BTCPlane *plane1, BTCPlane *plane2) {
  return plane1->points_size_ > plane2->points_size_;
}

void BTCOctoTree::init_octo_tree() {
  if (voxel_points_.size() > config_setting_.voxel_init_num_) {
    init_plane();
  }
}

void BTCOctoTree::init_plane() {
  plane_ptr_->covariance_ = Eigen::Matrix3d::Zero();
  plane_ptr_->center_ = Eigen::Vector3d::Zero();
  plane_ptr_->normal_ = Eigen::Vector3d::Zero();
  plane_ptr_->points_size_ = voxel_points_.size();
  plane_ptr_->radius_ = 0;
  for (auto &pi : voxel_points_) {
    plane_ptr_->covariance_ += pi * pi.transpose();
    plane_ptr_->center_ += pi;
  }
  plane_ptr_->center_ = plane_ptr_->center_ / plane_ptr_->points_size_;
  plane_ptr_->covariance_ =
      plane_ptr_->covariance_ / plane_ptr_->points_size_ -
      plane_ptr_->center_ * plane_ptr_->center_.transpose();
  Eigen::EigenSolver<Eigen::Matrix3d> es(plane_ptr_->covariance_);
  Eigen::Matrix3cd evecs = es.eigenvectors();
  Eigen::Vector3cd evals = es.eigenvalues();
  Eigen::Vector3d evalsReal;
  evalsReal = evals.real();
  Eigen::Matrix3d::Index evalsMin, evalsMax;
  evalsReal.rowwise().sum().minCoeff(&evalsMin);
  evalsReal.rowwise().sum().maxCoeff(&evalsMax);
  int evalsMid = 3 - evalsMin - evalsMax;
  if (evalsReal(evalsMin) < config_setting_.plane_detection_thre_) 
  {
    plane_ptr_->normal_ << evecs.real()(0, evalsMin), evecs.real()(1, evalsMin),
        evecs.real()(2, evalsMin);
    plane_ptr_->min_eigen_value_ = evalsReal(evalsMin);
    plane_ptr_->radius_ = sqrt(evalsReal(evalsMax));
    plane_ptr_->is_plane_ = true;

    plane_ptr_->d_ = -(plane_ptr_->normal_(0) * plane_ptr_->center_(0) +
                       plane_ptr_->normal_(1) * plane_ptr_->center_(1) +
                       plane_ptr_->normal_(2) * plane_ptr_->center_(2));
    plane_ptr_->p_center_.x = plane_ptr_->center_(0);
    plane_ptr_->p_center_.y = plane_ptr_->center_(1);
    plane_ptr_->p_center_.z = plane_ptr_->center_(2);
    plane_ptr_->p_center_.normal_x = plane_ptr_->normal_(0);
    plane_ptr_->p_center_.normal_y = plane_ptr_->normal_(1);
    plane_ptr_->p_center_.normal_z = plane_ptr_->normal_(2);
  } else {
    plane_ptr_->is_plane_ = false;
  }
}

// double
// calc_triangle_dis(const std::vector<std::pair<STD, STD>> &match_std_list) {
//   double mean_triangle_dis = 0;
//   for (auto var : match_std_list) {
//     mean_triangle_dis += (var.first.triangle_ - var.second.triangle_).norm() /
//                          var.first.triangle_.norm();
//   }
//   if (match_std_list.size() > 0) {
//     mean_triangle_dis = mean_triangle_dis / match_std_list.size();
//   } else {
//     mean_triangle_dis = -1;
//   }
//   return mean_triangle_dis;
// }

void STDescManager::GenerateSTDescs(
    pcl::PointCloud<pcl::PointXYZI>::Ptr &input_cloud,
    std::vector<STD> &stds_vec, int id) 
{ // step1, voxelization and plane dection
  std::unordered_map<BTCVOXEL_LOC, BTCOctoTree *> voxel_map;
  init_voxel_map(input_cloud, voxel_map);
  pcl::PointCloud<pcl::PointXYZINormal>::Ptr plane_cloud(
      new pcl::PointCloud<pcl::PointXYZINormal>);
  get_plane(voxel_map, plane_cloud);
  // std::cout << "[Description] planes size:" << plane_cloud->size() << std::endl;
  plane_cloud->header.seq = id;
  plane_cloud_vec_.push_back(plane_cloud);

  // step3, extraction binary descriptors
  std::vector<BTCPlane *> proj_plane_list;
  std::vector<BTCPlane *> merge_plane_list;
  get_project_plane(voxel_map, proj_plane_list);
  if (proj_plane_list.size() == 0) 
  {
    BTCPlane *single_plane = new BTCPlane;
    single_plane->normal_ << 0, 0, 1;
    single_plane->center_ << input_cloud->points[0].x, input_cloud->points[0].y,
        input_cloud->points[0].z;
    merge_plane_list.push_back(single_plane);
  } else {
    sort(proj_plane_list.begin(), proj_plane_list.end(), plane_greater_sort);
    merge_plane(proj_plane_list, merge_plane_list);
    sort(merge_plane_list.begin(), merge_plane_list.end(), plane_greater_sort);
  }
  std::vector<BinaryDescriptor> binary_list;
  binary_extractor(merge_plane_list, input_cloud, binary_list);
  // history_binary_list_.push_back(binary_list);
  // corner_cloud_vec_.push_back(corner_points);
  // std::cout << "[Description] corners size:" << binary_list.size() << std::endl;

  // vec_binary = binary_list;

  // step4, generate stable triangle descriptors
  stds_vec.clear();
  generate_std(binary_list, current_frame_id_, stds_vec);
  // std::cout << "[Description] stds size:" << stds_vec.size() << std::endl;

  // step5, clear memory
  for (auto iter = voxel_map.begin(); iter != voxel_map.end(); iter++) {
    delete (iter->second);
  }
  return;
}

void STDescManager::SearchLoop(
    std::vector<STD> &stds_vec, std::pair<int, double> &loop_result,
    std::pair<Eigen::Vector3d, Eigen::Matrix3d> &loop_transform,
    std::vector<std::pair<STD, STD>> &loop_std_pair, pcl::PointCloud<pcl::PointXYZINormal>::Ptr pl_cur) 
{
  if (stds_vec.size() == 0) {
    // ROS_ERROR_STREAM("No STDescs!");
    loop_result = std::pair<int, double>(-1, 0);
    return;
  }
  // step1, select candidates, default number 50
  std::vector<STDMatchList> candidate_matcher_vec;
  candidate_selector(stds_vec, candidate_matcher_vec);

  // step2, select best candidates from rough candidates
  double best_score = 0;
  int best_candidate_id = -1;
  int triggle_candidate = -1;
  std::pair<Eigen::Vector3d, Eigen::Matrix3d> best_transform;
  std::vector<std::pair<STD, STD>> best_sucess_match_vec;
  for (size_t i = 0; i < candidate_matcher_vec.size(); i++) 
  {
    double verify_score = -1;
    std::pair<Eigen::Vector3d, Eigen::Matrix3d> relative_pose;
    std::vector<std::pair<STD, STD>> sucess_match_vec;
    candidate_verify(candidate_matcher_vec[i], verify_score, relative_pose,
                     sucess_match_vec, pl_cur);
    // std::cout << "[Retreival] try frame:"
    //           << candidate_matcher_vec[i].match_id_.second
    //           << ", rough size:" << candidate_matcher_vec[i].match_list_.size()
    //           << ", score:" << verify_score << std::endl;
    if (verify_score > best_score) {
      best_score = verify_score;
      best_candidate_id = candidate_matcher_vec[i].match_id_.second;
      best_transform = relative_pose;
      best_sucess_match_vec = sucess_match_vec;
      triggle_candidate = i;
      // std::cout << "[Retreival] best candidate:" << best_candidate_id
      //           << ", score:" << best_score << std::endl;
    }
  }

  if (best_score > config_setting_.icp_threshold_) {
    loop_result = std::pair<int, double>(best_candidate_id, best_score);
    loop_transform = best_transform;
    loop_std_pair = best_sucess_match_vec;
    return;
  } else {
    loop_result = std::pair<int, double>(-1, 0);
    return;
  }
}

void STDescManager::AddSTDescs(const std::vector<STD> &stds_vec) {
  // update frame id
  current_frame_id_++;
  for (auto &single_std : stds_vec) {
    // calculate the position of single std
    STD_LOC position;
    position.x = (int)(single_std.triangle_[0] + 0.5);
    position.y = (int)(single_std.triangle_[1] + 0.5);
    position.z = (int)(single_std.triangle_[2] + 0.5);
    auto iter = data_base_.find(position);
    if (iter != data_base_.end()) {
      iter->second.push_back(single_std);
    } else {
      std::vector<STD> descriptor_vec;
      descriptor_vec.push_back(single_std);
      data_base_[position] = descriptor_vec;
    }
  }
  return;
}

void STDescManager::init_voxel_map(
    const pcl::PointCloud<pcl::PointXYZI>::Ptr &input_cloud,
    std::unordered_map<BTCVOXEL_LOC, BTCOctoTree *> &voxel_map) 
{
  uint plsize = input_cloud->size();
  for (uint i = 0; i < plsize; i++) {
    Eigen::Vector3d p_c(input_cloud->points[i].x, input_cloud->points[i].y,
                        input_cloud->points[i].z);
    double loc_xyz[3];
    for (int j = 0; j < 3; j++) {
      loc_xyz[j] = p_c[j] / config_setting_.voxel_size_;
      if (loc_xyz[j] < 0) {
        loc_xyz[j] -= 1.0;
      }
    }
    BTCVOXEL_LOC position((int64_t)loc_xyz[0], (int64_t)loc_xyz[1],
                       (int64_t)loc_xyz[2]);
    auto iter = voxel_map.find(position);
    if (iter != voxel_map.end()) {
      voxel_map[position]->voxel_points_.push_back(p_c);
    } else {
      BTCOctoTree *octo_tree = new BTCOctoTree(config_setting_);
      voxel_map[position] = octo_tree;
      voxel_map[position]->voxel_points_.push_back(p_c);
    }
  }
  std::vector<std::unordered_map<BTCVOXEL_LOC, BTCOctoTree *>::iterator> iter_list;
  std::vector<size_t> index;
  size_t i = 0;
  for (auto iter = voxel_map.begin(); iter != voxel_map.end(); ++iter) {
    index.push_back(i);
    i++;
    iter_list.push_back(iter);
    // iter->second->init_octo_tree();
  }
  // std::for_each(
  //     std::execution::par_unseq, index.begin(), index.end(),
  //     [&](const size_t &i) { iter_list[i]->second->init_octo_tree(); });

  for(const size_t &i: index)
    iter_list[i]->second->init_octo_tree();

}

void STDescManager::get_plane(
    const std::unordered_map<BTCVOXEL_LOC, BTCOctoTree *> &voxel_map,
    pcl::PointCloud<pcl::PointXYZINormal>::Ptr &plane_cloud) {
  for (auto iter = voxel_map.begin(); iter != voxel_map.end(); iter++) {
    if (iter->second->plane_ptr_->is_plane_) {
      pcl::PointXYZINormal pi;
      pi.x = iter->second->plane_ptr_->center_[0];
      pi.y = iter->second->plane_ptr_->center_[1];
      pi.z = iter->second->plane_ptr_->center_[2];
      pi.normal_x = iter->second->plane_ptr_->normal_[0];
      pi.normal_y = iter->second->plane_ptr_->normal_[1];
      pi.normal_z = iter->second->plane_ptr_->normal_[2];
      plane_cloud->push_back(pi);
    }
  }
}

void STDescManager::get_project_plane(
    std::unordered_map<BTCVOXEL_LOC, BTCOctoTree *> &voxel_map,
    std::vector<BTCPlane *> &project_plane_list) 
{
  std::vector<BTCPlane *> origin_list;
  for (auto iter = voxel_map.begin(); iter != voxel_map.end(); iter++) {
    if (iter->second->plane_ptr_->is_plane_) {
      origin_list.push_back(iter->second->plane_ptr_);
    }
  }
  for (size_t i = 0; i < origin_list.size(); i++)
    origin_list[i]->id_ = 0;
  int current_id = 1;
  for (auto iter = origin_list.end() - 1; iter != origin_list.begin(); iter--) 
  {
    for (auto iter2 = origin_list.begin(); iter2 != iter; iter2++) 
    {
      Eigen::Vector3d normal_diff = (*iter)->normal_ - (*iter2)->normal_;
      Eigen::Vector3d normal_add = (*iter)->normal_ + (*iter2)->normal_;
      double dis1 =
          fabs((*iter)->normal_(0) * (*iter2)->center_(0) +
               (*iter)->normal_(1) * (*iter2)->center_(1) +
               (*iter)->normal_(2) * (*iter2)->center_(2) + (*iter)->d_);
      double dis2 =
          fabs((*iter2)->normal_(0) * (*iter)->center_(0) +
               (*iter2)->normal_(1) * (*iter)->center_(1) +
               (*iter2)->normal_(2) * (*iter)->center_(2) + (*iter2)->d_);
      if (normal_diff.norm() < config_setting_.plane_merge_normal_thre_ ||
          normal_add.norm() < config_setting_.plane_merge_normal_thre_)
        if (dis1 < config_setting_.plane_merge_dis_thre_ &&
            dis2 < config_setting_.plane_merge_dis_thre_) {
          if ((*iter)->id_ == 0 && (*iter2)->id_ == 0) 
          {
            (*iter)->id_ = current_id;
            (*iter2)->id_ = current_id;
            current_id++;
          } else if ((*iter)->id_ == 0 && (*iter2)->id_ != 0)
            (*iter)->id_ = (*iter2)->id_;
          else if ((*iter)->id_ != 0 && (*iter2)->id_ == 0)
            (*iter2)->id_ = (*iter)->id_;
        }
    }
  }
  std::vector<BTCPlane *> merge_list;
  std::vector<int> merge_flag;

  for (size_t i = 0; i < origin_list.size(); i++) 
  {
    auto it =
        std::find(merge_flag.begin(), merge_flag.end(), origin_list[i]->id_);
    if (it != merge_flag.end())
      continue;
    if (origin_list[i]->id_ == 0) {
      continue;
    }
    BTCPlane *merge_plane = new BTCPlane;
    (*merge_plane) = (*origin_list[i]);
    bool is_merge = false;
    for (size_t j = 0; j < origin_list.size(); j++) 
    {
      if (i == j)
        continue;
      if (origin_list[j]->id_ == origin_list[i]->id_) 
      {
        is_merge = true;
        Eigen::Matrix3d P_PT1 =
            (merge_plane->covariance_ +
             merge_plane->center_ * merge_plane->center_.transpose()) *
            merge_plane->points_size_;
        Eigen::Matrix3d P_PT2 =
            (origin_list[j]->covariance_ +
             origin_list[j]->center_ * origin_list[j]->center_.transpose()) *
            origin_list[j]->points_size_;
        Eigen::Vector3d merge_center =
            (merge_plane->center_ * merge_plane->points_size_ +
             origin_list[j]->center_ * origin_list[j]->points_size_) /
            (merge_plane->points_size_ + origin_list[j]->points_size_);
        Eigen::Matrix3d merge_covariance =
            (P_PT1 + P_PT2) /
                (merge_plane->points_size_ + origin_list[j]->points_size_) -
            merge_center * merge_center.transpose();
        merge_plane->covariance_ = merge_covariance;
        merge_plane->center_ = merge_center;
        merge_plane->points_size_ =
            merge_plane->points_size_ + origin_list[j]->points_size_;
        merge_plane->sub_plane_num_++;
        // for (size_t k = 0; k < origin_list[j]->cloud.size(); k++) {
        //   merge_plane->cloud.points.push_back(origin_list[j]->cloud.points[k]);
        // }
        Eigen::EigenSolver<Eigen::Matrix3d> es(merge_plane->covariance_);
        Eigen::Matrix3cd evecs = es.eigenvectors();
        Eigen::Vector3cd evals = es.eigenvalues();
        Eigen::Vector3d evalsReal;
        evalsReal = evals.real();
        Eigen::Matrix3f::Index evalsMin, evalsMax;
        evalsReal.rowwise().sum().minCoeff(&evalsMin);
        evalsReal.rowwise().sum().maxCoeff(&evalsMax);
        Eigen::Vector3d evecMin = evecs.real().col(evalsMin);
        merge_plane->normal_ << evecs.real()(0, evalsMin),
            evecs.real()(1, evalsMin), evecs.real()(2, evalsMin);
        merge_plane->radius_ = sqrt(evalsReal(evalsMax));
        merge_plane->d_ = -(merge_plane->normal_(0) * merge_plane->center_(0) +
                            merge_plane->normal_(1) * merge_plane->center_(1) +
                            merge_plane->normal_(2) * merge_plane->center_(2));
        merge_plane->p_center_.x = merge_plane->center_(0);
        merge_plane->p_center_.y = merge_plane->center_(1);
        merge_plane->p_center_.z = merge_plane->center_(2);
        merge_plane->p_center_.normal_x = merge_plane->normal_(0);
        merge_plane->p_center_.normal_y = merge_plane->normal_(1);
        merge_plane->p_center_.normal_z = merge_plane->normal_(2);
      }
    }
    if (is_merge) {
      merge_flag.push_back(merge_plane->id_);
      merge_list.push_back(merge_plane);
    }
  }
  project_plane_list = merge_list;
}

void STDescManager::merge_plane(std::vector<BTCPlane *> &origin_list,
                                std::vector<BTCPlane *> &merge_plane_list) {
  if (origin_list.size() == 1) {
    merge_plane_list = origin_list;
    return;
  }
  for (size_t i = 0; i < origin_list.size(); i++)
    origin_list[i]->id_ = 0;
  int current_id = 1;
  for (auto iter = origin_list.end() - 1; iter != origin_list.begin(); iter--) 
  {
    for (auto iter2 = origin_list.begin(); iter2 != iter; iter2++) 
    {
      Eigen::Vector3d normal_diff = (*iter)->normal_ - (*iter2)->normal_;
      Eigen::Vector3d normal_add = (*iter)->normal_ + (*iter2)->normal_;
      double dis1 =
          fabs((*iter)->normal_(0) * (*iter2)->center_(0) +
               (*iter)->normal_(1) * (*iter2)->center_(1) +
               (*iter)->normal_(2) * (*iter2)->center_(2) + (*iter)->d_);
      double dis2 =
          fabs((*iter2)->normal_(0) * (*iter)->center_(0) +
               (*iter2)->normal_(1) * (*iter)->center_(1) +
               (*iter2)->normal_(2) * (*iter)->center_(2) + (*iter2)->d_);
      if (normal_diff.norm() < config_setting_.plane_merge_normal_thre_ ||
          normal_add.norm() < config_setting_.plane_merge_normal_thre_)
        if (dis1 < config_setting_.plane_merge_dis_thre_ &&
            dis2 < config_setting_.plane_merge_dis_thre_) {
          if ((*iter)->id_ == 0 && (*iter2)->id_ == 0) {
            (*iter)->id_ = current_id;
            (*iter2)->id_ = current_id;
            current_id++;
          } else if ((*iter)->id_ == 0 && (*iter2)->id_ != 0)
            (*iter)->id_ = (*iter2)->id_;
          else if ((*iter)->id_ != 0 && (*iter2)->id_ == 0)
            (*iter2)->id_ = (*iter)->id_;
        }
    }
  }
  std::vector<int> merge_flag;

  for (size_t i = 0; i < origin_list.size(); i++) {
    auto it =
        std::find(merge_flag.begin(), merge_flag.end(), origin_list[i]->id_);
    if (it != merge_flag.end())
      continue;
    if (origin_list[i]->id_ == 0) {
      merge_plane_list.push_back(origin_list[i]);
      continue;
    }
    BTCPlane *merge_plane = new BTCPlane;
    (*merge_plane) = (*origin_list[i]);
    bool is_merge = false;
    for (size_t j = 0; j < origin_list.size(); j++) {
      if (i == j)
        continue;
      if (origin_list[j]->id_ == origin_list[i]->id_) {
        is_merge = true;
        Eigen::Matrix3d P_PT1 =
            (merge_plane->covariance_ +
             merge_plane->center_ * merge_plane->center_.transpose()) *
            merge_plane->points_size_;
        Eigen::Matrix3d P_PT2 =
            (origin_list[j]->covariance_ +
             origin_list[j]->center_ * origin_list[j]->center_.transpose()) *
            origin_list[j]->points_size_;
        Eigen::Vector3d merge_center =
            (merge_plane->center_ * merge_plane->points_size_ +
             origin_list[j]->center_ * origin_list[j]->points_size_) /
            (merge_plane->points_size_ + origin_list[j]->points_size_);
        Eigen::Matrix3d merge_covariance =
            (P_PT1 + P_PT2) /
                (merge_plane->points_size_ + origin_list[j]->points_size_) -
            merge_center * merge_center.transpose();
        merge_plane->covariance_ = merge_covariance;
        merge_plane->center_ = merge_center;
        merge_plane->points_size_ =
            merge_plane->points_size_ + origin_list[j]->points_size_;
        merge_plane->sub_plane_num_ += origin_list[j]->sub_plane_num_;
        // for (size_t k = 0; k < origin_list[j]->cloud.size(); k++) {
        //   merge_plane->cloud.points.push_back(origin_list[j]->cloud.points[k]);
        // }
        Eigen::EigenSolver<Eigen::Matrix3d> es(merge_plane->covariance_);
        Eigen::Matrix3cd evecs = es.eigenvectors();
        Eigen::Vector3cd evals = es.eigenvalues();
        Eigen::Vector3d evalsReal;
        evalsReal = evals.real();
        Eigen::Matrix3f::Index evalsMin, evalsMax;
        evalsReal.rowwise().sum().minCoeff(&evalsMin);
        evalsReal.rowwise().sum().maxCoeff(&evalsMax);
        Eigen::Vector3d evecMin = evecs.real().col(evalsMin);
        merge_plane->normal_ << evecs.real()(0, evalsMin),
            evecs.real()(1, evalsMin), evecs.real()(2, evalsMin);
        merge_plane->radius_ = sqrt(evalsReal(evalsMax));
        merge_plane->d_ = -(merge_plane->normal_(0) * merge_plane->center_(0) +
                            merge_plane->normal_(1) * merge_plane->center_(1) +
                            merge_plane->normal_(2) * merge_plane->center_(2));
        merge_plane->p_center_.x = merge_plane->center_(0);
        merge_plane->p_center_.y = merge_plane->center_(1);
        merge_plane->p_center_.z = merge_plane->center_(2);
        merge_plane->p_center_.normal_x = merge_plane->normal_(0);
        merge_plane->p_center_.normal_y = merge_plane->normal_(1);
        merge_plane->p_center_.normal_z = merge_plane->normal_(2);
      }
    }
    if (is_merge) {
      merge_flag.push_back(merge_plane->id_);
      merge_plane_list.push_back(merge_plane);
    }
  }
}

void STDescManager::binary_extractor(
    const std::vector<BTCPlane *> proj_plane_list,
    const pcl::PointCloud<pcl::PointXYZI>::Ptr &input_cloud,
    std::vector<BinaryDescriptor> &binary_descriptor_list) 
{
  binary_descriptor_list.clear();
  std::vector<BinaryDescriptor> temp_binary_list;
  Eigen::Vector3d last_normal(0, 0, 0);
  int useful_proj_num = 0;
  for (int i = 0; i < proj_plane_list.size(); i++) 
  {
    std::vector<BinaryDescriptor> prepare_binary_list;
    Eigen::Vector3d proj_center = proj_plane_list[i]->center_;
    Eigen::Vector3d proj_normal = proj_plane_list[i]->normal_;
    if ((proj_normal - last_normal).norm() < 0.3 ||
        (proj_normal + last_normal).norm() > 0.3) {
      last_normal = proj_normal;
      // std::cout << "proj normal:" << proj_normal.transpose() << std::endl;
      useful_proj_num++;
      extract_binary(proj_center, proj_normal, input_cloud,
                     prepare_binary_list);
      for (auto &bi : prepare_binary_list) {
        temp_binary_list.push_back(bi);
      }
      if (useful_proj_num == config_setting_.proj_plane_num_) {
        break;
      }
    }
  }
  non_maxi_suppression(temp_binary_list);
  if (config_setting_.useful_corner_num_ > temp_binary_list.size()) {
    binary_descriptor_list = temp_binary_list;
  } else {
    std::sort(temp_binary_list.begin(), temp_binary_list.end(),
              binary_greater_sort);
    for (size_t i = 0; i < config_setting_.useful_corner_num_; i++) {
      binary_descriptor_list.push_back(temp_binary_list[i]);
    }
  }
  return;
}

void STDescManager::extract_binary(
    const Eigen::Vector3d &project_center,
    const Eigen::Vector3d &project_normal,
    const pcl::PointCloud<pcl::PointXYZI>::Ptr &input_cloud,
    std::vector<BinaryDescriptor> &binary_list) 
{
  binary_list.clear();
  double binary_min_dis = config_setting_.summary_min_thre_;
  double resolution = config_setting_.proj_image_resolution_;
  double dis_threshold_min = config_setting_.proj_dis_min_;
  double dis_threshold_max = config_setting_.proj_dis_max_;
  double high_inc = config_setting_.proj_image_high_inc_;
  bool line_filter_enable = config_setting_.line_filter_enable_;
  double A = project_normal[0];
  double B = project_normal[1];
  double C = project_normal[2];
  double D =
      -(A * project_center[0] + B * project_center[1] + C * project_center[2]);
  std::vector<Eigen::Vector3d> projection_points;
  Eigen::Vector3d x_axis(1, 1, 0);
  // Eigen::Vector3d x_axis(1, 0, 0);
  if (C != 0) {
    x_axis[2] = -(A + B) / C;
  } else if (B != 0) {
    x_axis[1] = -A / B;
  } else {
    x_axis[0] = 0;
    x_axis[1] = 1;
  }
  x_axis.normalize();
  Eigen::Vector3d y_axis = project_normal.cross(x_axis);
  y_axis.normalize();
  double ax = x_axis[0];
  double bx = x_axis[1];
  double cx = x_axis[2];
  double dx = -(ax * project_center[0] + bx * project_center[1] +
                cx * project_center[2]);
  double ay = y_axis[0];
  double by = y_axis[1];
  double cy = y_axis[2];
  double dy = -(ay * project_center[0] + by * project_center[1] +
                cy * project_center[2]);
  std::vector<Eigen::Vector2d> point_list_2d;
  pcl::PointCloud<pcl::PointXYZ> point_list_3d;
  std::vector<double> dis_list_2d;
  for (size_t i = 0; i < input_cloud->size(); i++) 
  {
    double x = input_cloud->points[i].x;
    double y = input_cloud->points[i].y;
    double z = input_cloud->points[i].z;
    double dis = fabs(x * A + y * B + z * C + D);
    pcl::PointXYZ pi;
    if (dis < dis_threshold_min || dis > dis_threshold_max) {
      continue;
    } else {
      if (dis > dis_threshold_min && dis <= dis_threshold_max) {
        pi.x = x;
        pi.y = y;
        pi.z = z;
      }
    }
    Eigen::Vector3d cur_project;

    cur_project[0] = (-A * (B * y + C * z + D) + x * (B * B + C * C)) /
                     (A * A + B * B + C * C);
    cur_project[1] = (-B * (A * x + C * z + D) + y * (A * A + C * C)) /
                     (A * A + B * B + C * C);
    cur_project[2] = (-C * (A * x + B * y + D) + z * (A * A + B * B)) /
                     (A * A + B * B + C * C);
    pcl::PointXYZ p;
    p.x = cur_project[0];
    p.y = cur_project[1];
    p.z = cur_project[2];
    double project_x =
        cur_project[0] * ay + cur_project[1] * by + cur_project[2] * cy + dy;
    double project_y =
        cur_project[0] * ax + cur_project[1] * bx + cur_project[2] * cx + dx;
    Eigen::Vector2d p_2d(project_x, project_y);
    point_list_2d.push_back(p_2d);
    dis_list_2d.push_back(dis);
    point_list_3d.points.push_back(pi);
  }
  double min_x = 10;
  double max_x = -10;
  double min_y = 10;
  double max_y = -10;
  if (point_list_2d.size() <= 5) {
    return;
  }
  for (auto &pi : point_list_2d) {
    if (pi[0] < min_x) {
      min_x = pi[0];
    }
    if (pi[0] > max_x) {
      max_x = pi[0];
    }
    if (pi[1] < min_y) {
      min_y = pi[1];
    }
    if (pi[1] > max_y) {
      max_y = pi[1];
    }
  }
  // segment project cloud
  int segmen_base_num = 5;
  double segmen_len = segmen_base_num * resolution;
  int x_segment_num = (max_x - min_x) / segmen_len + 1;
  int y_segment_num = (max_y - min_y) / segmen_len + 1;
  int x_axis_len = (int)((max_x - min_x) / resolution + segmen_base_num);
  int y_axis_len = (int)((max_y - min_y) / resolution + segmen_base_num);

  std::vector<double> **dis_container = new std::vector<double> *[x_axis_len];
  BinaryDescriptor **binary_container = new BinaryDescriptor *[x_axis_len];
  for (int i = 0; i < x_axis_len; i++) {
    dis_container[i] = new std::vector<double>[y_axis_len];
    binary_container[i] = new BinaryDescriptor[y_axis_len];
  }
  double **img_count = new double *[x_axis_len];
  for (int i = 0; i < x_axis_len; i++) {
    img_count[i] = new double[y_axis_len];
  }
  double **dis_array = new double *[x_axis_len];
  for (int i = 0; i < x_axis_len; i++) {
    dis_array[i] = new double[y_axis_len];
  }
  double **mean_x_list = new double *[x_axis_len];
  for (int i = 0; i < x_axis_len; i++) {
    mean_x_list[i] = new double[y_axis_len];
  }
  double **mean_y_list = new double *[x_axis_len];
  for (int i = 0; i < x_axis_len; i++) {
    mean_y_list[i] = new double[y_axis_len];
  }
  for (int x = 0; x < x_axis_len; x++) {
    for (int y = 0; y < y_axis_len; y++) {
      img_count[x][y] = 0;
      mean_x_list[x][y] = 0;
      mean_y_list[x][y] = 0;
      dis_array[x][y] = 0;
      std::vector<double> single_dis_container;
      dis_container[x][y] = single_dis_container;
    }
  }

  for (size_t i = 0; i < point_list_2d.size(); i++) {
    int x_index = (int)((point_list_2d[i][0] - min_x) / resolution);
    int y_index = (int)((point_list_2d[i][1] - min_y) / resolution);
    mean_x_list[x_index][y_index] += point_list_2d[i][0];
    mean_y_list[x_index][y_index] += point_list_2d[i][1];
    img_count[x_index][y_index]++;
    dis_container[x_index][y_index].push_back(dis_list_2d[i]);
  }

  for (int x = 0; x < x_axis_len; x++) {
    for (int y = 0; y < y_axis_len; y++) {
      // calc segment dis array
      if (img_count[x][y] > 0) {
        int cut_num = (dis_threshold_max - dis_threshold_min) / high_inc;
        std::vector<bool> occup_list;
        std::vector<double> cnt_list;
        BinaryDescriptor single_binary;
        for (size_t i = 0; i < cut_num; i++) {
          cnt_list.push_back(0);
          occup_list.push_back(false);
        }
        for (size_t j = 0; j < dis_container[x][y].size(); j++) {
          int cnt_index =
              (dis_container[x][y][j] - dis_threshold_min) / high_inc;
          cnt_list[cnt_index]++;
        }
        double segmnt_dis = 0;
        for (size_t i = 0; i < cut_num; i++) {
          if (cnt_list[i] >= 1) {
            segmnt_dis++;
            occup_list[i] = true;
          }
        }
        dis_array[x][y] = segmnt_dis;
        single_binary.occupy_array_ = occup_list;
        single_binary.summary_ = segmnt_dis;
        binary_container[x][y] = single_binary;
      }
    }
  }

  // filter by distance
  std::vector<double> max_dis_list;
  std::vector<int> max_dis_x_index_list;
  std::vector<int> max_dis_y_index_list;

  for (int x_segment_index = 0; x_segment_index < x_segment_num;
       x_segment_index++) {
    for (int y_segment_index = 0; y_segment_index < y_segment_num;
         y_segment_index++) {
      double max_dis = 0;
      int max_dis_x_index = -10;
      int max_dis_y_index = -10;
      for (int x_index = x_segment_index * segmen_base_num;
           x_index < (x_segment_index + 1) * segmen_base_num; x_index++) {
        for (int y_index = y_segment_index * segmen_base_num;
             y_index < (y_segment_index + 1) * segmen_base_num; y_index++) {
          if (dis_array[x_index][y_index] > max_dis) {
            max_dis = dis_array[x_index][y_index];
            max_dis_x_index = x_index;
            max_dis_y_index = y_index;
          }
        }
      }
      if (max_dis >= binary_min_dis) {
        bool is_touch = true;
        if (config_setting_.touch_filter_enable_) {
          is_touch = binary_container[max_dis_x_index][max_dis_y_index]
                         .occupy_array_[0] ||
                     binary_container[max_dis_x_index][max_dis_y_index]
                         .occupy_array_[1] ||
                     binary_container[max_dis_x_index][max_dis_y_index]
                         .occupy_array_[2] ||
                     binary_container[max_dis_x_index][max_dis_y_index]
                         .occupy_array_[3];
        }

        if (is_touch) {
          max_dis_list.push_back(max_dis);
          max_dis_x_index_list.push_back(max_dis_x_index);
          max_dis_y_index_list.push_back(max_dis_y_index);
        }
      }
    }
  }
  // calc line or not
  std::vector<Eigen::Vector2i> direction_list;
  Eigen::Vector2i d(0, 1);
  direction_list.push_back(d);
  d << 1, 0;
  direction_list.push_back(d);
  d << 1, 1;
  direction_list.push_back(d);
  d << 1, -1;
  direction_list.push_back(d);
  for (size_t i = 0; i < max_dis_list.size(); i++) {
    Eigen::Vector2i p(max_dis_x_index_list[i], max_dis_y_index_list[i]);
    if (p[0] <= 0 || p[0] >= x_axis_len - 1 || p[1] <= 0 ||
        p[1] >= y_axis_len - 1) {
      continue;
    }
    bool is_add = true;

    if (line_filter_enable) {
      for (int j = 0; j < 4; j++) {
        Eigen::Vector2i p(max_dis_x_index_list[i], max_dis_y_index_list[i]);
        if (p[0] <= 0 || p[0] >= x_axis_len - 1 || p[1] <= 0 ||
            p[1] >= y_axis_len - 1) {
          continue;
        }
        Eigen::Vector2i p1 = p + direction_list[j];
        Eigen::Vector2i p2 = p - direction_list[j];
        double threshold = dis_array[p[0]][p[1]] - 3;
        if (dis_array[p1[0]][p1[1]] >= threshold) {
          if (dis_array[p2[0]][p2[1]] >= 0.5 * dis_array[p[0]][p[1]]) {
            is_add = false;
          }
        }
        if (dis_array[p2[0]][p2[1]] >= threshold) {
          if (dis_array[p1[0]][p1[1]] >= 0.5 * dis_array[p[0]][p[1]]) {
            is_add = false;
          }
        }
        if (dis_array[p1[0]][p1[1]] >= threshold) {
          if (dis_array[p2[0]][p2[1]] >= threshold) {
            is_add = false;
          }
        }
        if (dis_array[p2[0]][p2[1]] >= threshold) {
          if (dis_array[p1[0]][p1[1]] >= threshold) {
            is_add = false;
          }
        }
      }
    }
    if (is_add) {
      double px =
          mean_x_list[max_dis_x_index_list[i]][max_dis_y_index_list[i]] /
          img_count[max_dis_x_index_list[i]][max_dis_y_index_list[i]];
      double py =
          mean_y_list[max_dis_x_index_list[i]][max_dis_y_index_list[i]] /
          img_count[max_dis_x_index_list[i]][max_dis_y_index_list[i]];
      Eigen::Vector3d coord = py * x_axis + px * y_axis + project_center;
      pcl::PointXYZ pi;
      pi.x = coord[0];
      pi.y = coord[1];
      pi.z = coord[2];
      BinaryDescriptor single_binary =
          binary_container[max_dis_x_index_list[i]][max_dis_y_index_list[i]];
      single_binary.location_ = coord;
      binary_list.push_back(single_binary);
    }
  }
  for (int i = 0; i < x_axis_len; i++) {
    delete[] binary_container[i];
    delete[] dis_container[i];
    delete[] img_count[i];
    delete[] dis_array[i];
    delete[] mean_x_list[i];
    delete[] mean_y_list[i];
  }
  delete[] binary_container;
  delete[] dis_container;
  delete[] img_count;
  delete[] dis_array;
  delete[] mean_x_list;
  delete[] mean_y_list;
}

void STDescManager::non_maxi_suppression(
    std::vector<BinaryDescriptor> &binary_list) {
  pcl::PointCloud<pcl::PointXYZ>::Ptr prepare_key_cloud(
      new pcl::PointCloud<pcl::PointXYZ>);
  pcl::KdTreeFLANN<pcl::PointXYZ> kd_tree;
  std::vector<int> pre_count_list;
  std::vector<bool> is_add_list;
  for (auto &var : binary_list) {
    pcl::PointXYZ pi;
    pi.x = var.location_[0];
    pi.y = var.location_[1];
    pi.z = var.location_[2];
    prepare_key_cloud->push_back(pi);
    pre_count_list.push_back(var.summary_);
    is_add_list.push_back(true);
  }
  if(prepare_key_cloud->size() == 0) return;
  kd_tree.setInputCloud(prepare_key_cloud);
  std::vector<int> pointIdxRadiusSearch;
  std::vector<float> pointRadiusSquaredDistance;
  double radius = config_setting_.non_max_suppression_radius_;
  for (size_t i = 0; i < prepare_key_cloud->size(); i++) {
    pcl::PointXYZ searchPoint = prepare_key_cloud->points[i];
    if (kd_tree.radiusSearch(searchPoint, radius, pointIdxRadiusSearch,
                             pointRadiusSquaredDistance) > 0) {
      Eigen::Vector3d pi(searchPoint.x, searchPoint.y, searchPoint.z);
      for (size_t j = 0; j < pointIdxRadiusSearch.size(); ++j) {
        // Eigen::Vector3d pj(
        //     prepare_key_cloud->points[pointIdxRadiusSearch[j]].x,
        //     prepare_key_cloud->points[pointIdxRadiusSearch[j]].y,
        //     prepare_key_cloud->points[pointIdxRadiusSearch[j]].z);
        if (pointIdxRadiusSearch[j] == i) {
          continue;
        }
        if (pre_count_list[i] <= pre_count_list[pointIdxRadiusSearch[j]]) {
          is_add_list[i] = false;
        }
      }
    }
  }
  std::vector<BinaryDescriptor> pass_binary_list;
  for (size_t i = 0; i < is_add_list.size(); i++) {
    if (is_add_list[i]) {
      pass_binary_list.push_back(binary_list[i]);
    }
  }
  binary_list.clear();
  for (auto &var : pass_binary_list) {
    binary_list.push_back(var);
  }
  return;
}

void STDescManager::generate_std(
    const std::vector<BinaryDescriptor> &binary_list, const int &frame_id,
    std::vector<STD> &std_list) {
  double scale = 1.0 / config_setting_.std_side_resolution_;
  std::unordered_map<BTCVOXEL_LOC, bool> feat_map;
  pcl::PointCloud<pcl::PointXYZ> key_cloud;
  for (auto &var : binary_list) {
    pcl::PointXYZ pi;
    pi.x = var.location_[0];
    pi.y = var.location_[1];
    pi.z = var.location_[2];
    key_cloud.push_back(pi);
  }
  pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr kd_tree(
      new pcl::KdTreeFLANN<pcl::PointXYZ>);
  kd_tree->setInputCloud(key_cloud.makeShared());
  int K = config_setting_.descriptor_near_num_;
  std::vector<int> pointIdxNKNSearch(K);
  std::vector<float> pointNKNSquaredDistance(K);
  for (size_t i = 0; i < key_cloud.size(); i++) 
  {
    pcl::PointXYZ searchPoint = key_cloud.points[i];
    if (kd_tree->nearestKSearch(searchPoint, K, pointIdxNKNSearch,
                                pointNKNSquaredDistance) > 0) 
    {
      for (int m = 1; m < K - 1; m++) 
      {
        for (int n = m + 1; n < K; n++) 
        {
          pcl::PointXYZ p1 = searchPoint;
          pcl::PointXYZ p2 = key_cloud.points[pointIdxNKNSearch[m]];
          pcl::PointXYZ p3 = key_cloud.points[pointIdxNKNSearch[n]];
          double a = sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) +
                          pow(p1.z - p2.z, 2));
          double b = sqrt(pow(p1.x - p3.x, 2) + pow(p1.y - p3.y, 2) +
                          pow(p1.z - p3.z, 2));
          double c = sqrt(pow(p3.x - p2.x, 2) + pow(p3.y - p2.y, 2) +
                          pow(p3.z - p2.z, 2));
          if (a > config_setting_.descriptor_max_len_ ||
              b > config_setting_.descriptor_max_len_ ||
              c > config_setting_.descriptor_max_len_ ||
              a < config_setting_.descriptor_min_len_ ||
              b < config_setting_.descriptor_min_len_ ||
              c < config_setting_.descriptor_min_len_) {
            continue;
          }
          double temp;
          Eigen::Vector3d A, B, C;
          Eigen::Vector3i l1, l2, l3;
          Eigen::Vector3i l_temp;
          l1 << 1, 2, 0;
          l2 << 1, 0, 3;
          l3 << 0, 2, 3;
          if (a > b) {
            temp = a;
            a = b;
            b = temp;
            l_temp = l1;
            l1 = l2;
            l2 = l_temp;
          }
          if (b > c) {
            temp = b;
            b = c;
            c = temp;
            l_temp = l2;
            l2 = l3;
            l3 = l_temp;
          }
          if (a > b) {
            temp = a;
            a = b;
            b = temp;
            l_temp = l1;
            l1 = l2;
            l2 = l_temp;
          }
          if (fabs(c - (a + b)) < 0.2) {
            continue;
          }

          pcl::PointXYZ d_p;
          d_p.x = a * 1000;
          d_p.y = b * 1000;
          d_p.z = c * 1000;
          BTCVOXEL_LOC position((int64_t)d_p.x, (int64_t)d_p.y, (int64_t)d_p.z);
          auto iter = feat_map.find(position);
          Eigen::Vector3d normal_1, normal_2, normal_3;
          BinaryDescriptor binary_A;
          BinaryDescriptor binary_B;
          BinaryDescriptor binary_C;
          if (iter == feat_map.end()) {
            if (l1[0] == l2[0]) {
              A << p1.x, p1.y, p1.z;
              binary_A = binary_list[i];
            } else if (l1[1] == l2[1]) {
              A << p2.x, p2.y, p2.z;
              binary_A = binary_list[pointIdxNKNSearch[m]];
            } else {
              A << p3.x, p3.y, p3.z;
              binary_A = binary_list[pointIdxNKNSearch[n]];
            }
            if (l1[0] == l3[0]) {
              B << p1.x, p1.y, p1.z;
              binary_B = binary_list[i];
            } else if (l1[1] == l3[1]) {
              B << p2.x, p2.y, p2.z;
              binary_B = binary_list[pointIdxNKNSearch[m]];
            } else {
              B << p3.x, p3.y, p3.z;
              binary_B = binary_list[pointIdxNKNSearch[n]];
            }
            if (l2[0] == l3[0]) {
              C << p1.x, p1.y, p1.z;
              binary_C = binary_list[i];
            } else if (l2[1] == l3[1]) {
              C << p2.x, p2.y, p2.z;
              binary_C = binary_list[pointIdxNKNSearch[m]];
            } else {
              C << p3.x, p3.y, p3.z;
              binary_C = binary_list[pointIdxNKNSearch[n]];
            }
            STD single_descriptor;
            single_descriptor.binary_A_ = binary_A;
            single_descriptor.binary_B_ = binary_B;
            single_descriptor.binary_C_ = binary_C;
            single_descriptor.center_ = (A + B + C) / 3;
            single_descriptor.triangle_ << scale * a, scale * b, scale * c;
            single_descriptor.angle_[0] = fabs(5 * normal_1.dot(normal_2));
            single_descriptor.angle_[1] = fabs(5 * normal_1.dot(normal_3));
            single_descriptor.angle_[2] = fabs(5 * normal_3.dot(normal_2));
            // single_descriptor.angle << 0, 0, 0;
            single_descriptor.frame_number_ = frame_id;
            // single_descriptor.score_frame_.push_back(frame_number);
            // Eigen::Matrix3d triangle_positon;
            // triangle_positon.block<3, 1>(0, 0) = A;
            // triangle_positon.block<3, 1>(0, 1) = B;
            // triangle_positon.block<3, 1>(0, 2) = C;
            // single_descriptor.position_list_.push_back(triangle_positon);
            // single_descriptor.triangle_scale_ = scale;
            feat_map[position] = true;
            std_list.push_back(single_descriptor);
          }
        }
      }
    }
  }
}

void STDescManager::candidate_selector(
    std::vector<STD> &current_STD_list,
    std::vector<STDMatchList> &candidate_matcher_vec) 
{
  std::vector<double> match_array(plane_cloud_vec_.size(), 0);
  std::vector<int> match_list_index;
  std::vector<Eigen::Vector3i> voxel_round;
  for (int x = -1; x <= 1; x++) {
    for (int y = -1; y <= 1; y++) {
      for (int z = -1; z <= 1; z++) {
        Eigen::Vector3i voxel_inc(x, y, z);
        voxel_round.push_back(voxel_inc);
      }
    }
  }
  std::vector<bool> useful_match(current_STD_list.size(), false);
  std::vector<std::vector<size_t>> useful_match_index(current_STD_list.size());
  std::vector<std::vector<STD_LOC>> useful_match_position(
      current_STD_list.size());
  std::vector<size_t> index(current_STD_list.size());
  // for (size_t i = 0; i < index.size(); ++i) {
  //   index[i] = i;
  //   useful_match[i] = false;
  // }
  // std::mutex mylock;
  // auto t0 = std::chrono::high_resolution_clock::now();

  // std::for_each(
  //     std::execution::par_unseq, index.begin(), index.end(),
  //     [&](const size_t &i) 
  // for(const size_t &i: index)
  for (size_t i = 0; i < current_STD_list.size(); ++i) 
  {
    STD &descriptor = current_STD_list[i];
    STD_LOC position;
    // int best_index = 0;
    // STD_LOC best_position;
    double dis_threshold =
        descriptor.triangle_.norm() *
        config_setting_.rough_dis_threshold_; // old 0.005
    for (auto &voxel_inc : voxel_round) 
    {
      position.x = (int)(descriptor.triangle_[0] + voxel_inc[0]);
      position.y = (int)(descriptor.triangle_[1] + voxel_inc[1]);
      position.z = (int)(descriptor.triangle_[2] + voxel_inc[2]);
      Eigen::Vector3d voxel_center((double)position.x + 0.5,
                                    (double)position.y + 0.5,
                                    (double)position.z + 0.5);
      if ((descriptor.triangle_ - voxel_center).norm() < 1.5) 
      {
        auto iter = data_base_.find(position);
        if (iter != data_base_.end()) 
        {
          bool is_push_position = false;
          for (size_t j = 0; j < iter->second.size(); j++) 
          {
            if ((descriptor.frame_number_ -
                  iter->second[j].frame_number_) >
                config_setting_.skip_near_num_) 
            {
              double dis =
                  (descriptor.triangle_ - iter->second[j].triangle_)
                      .norm();
              if (dis < dis_threshold) {
                double similarity =
                    (binary_similarity(descriptor.binary_A_,
                                        iter->second[j].binary_A_) +
                      binary_similarity(descriptor.binary_B_,
                                        iter->second[j].binary_B_) +
                      binary_similarity(descriptor.binary_C_,
                                        iter->second[j].binary_C_)) /
                    3;
                if (similarity > config_setting_.similarity_threshold_) {
                  useful_match[i] = true;
                  useful_match_position[i].push_back(position);
                  useful_match_index[i].push_back(j);
                }
              }
            }
          }
        }
      }
    }
  }
  // );
  std::vector<Eigen::Vector2i, Eigen::aligned_allocator<Eigen::Vector2i>>
      index_recorder;
  
  for (size_t i = 0; i < useful_match.size(); i++) 
  {
    if (useful_match[i]) 
    {
      for (size_t j = 0; j < useful_match_index[i].size(); j++) 
      {
        match_array[data_base_[useful_match_position[i][j]]
                              [useful_match_index[i][j]]
                                  .frame_number_] += 1;
        Eigen::Vector2i match_index(i, j);
        index_recorder.push_back(match_index);
        // match_list.push_back(single_match_pair);
        match_list_index.push_back(
            data_base_[useful_match_position[i][j]][useful_match_index[i][j]]
                .frame_number_);
      }
    }
  }

  for (int cnt = 0; cnt < config_setting_.candidate_num_; cnt++) 
  {
    // double max_vote = 1;
    // int max_vote_index = -1;
    // for (int i = 0; i < 20000; i++) 
    // for(int i=0; i<match_array.size(); i++)
    // {
    //   if (match_array[i] > max_vote) {
    //     max_vote = match_array[i];
    //     max_vote_index = i;
    //   }
    // }

    double max_vote_index = std::max_element(match_array.begin(), match_array.end()) - match_array.begin();
    int max_vote = match_array[max_vote_index];

    STDMatchList match_triangle_list;
    if (max_vote_index >= 0 && max_vote >= 5) 
    {
      match_array[max_vote_index] = 0;
      match_triangle_list.match_frame_ = max_vote_index;
      match_triangle_list.match_id_.first = current_frame_id_;
      match_triangle_list.match_id_.second = max_vote_index;
      double mean_dis = 0;
      for (size_t i = 0; i < index_recorder.size(); i++) 
      {
        if (match_list_index[i] == max_vote_index) 
        {
          std::pair<STD, STD> single_match_pair;
          single_match_pair.first = current_STD_list[index_recorder[i][0]];
          single_match_pair.second =
              data_base_[useful_match_position[index_recorder[i][0]]
                                              [index_recorder[i][1]]]
                        [useful_match_index[index_recorder[i][0]]
                                           [index_recorder[i][1]]];
          match_triangle_list.match_list_.push_back(single_match_pair);
        }
      }
      candidate_matcher_vec.push_back(match_triangle_list);
    }
    else
      break;
    
  }
}

void STDescManager::candidate_verify(
    STDMatchList &candidate_matcher, double &verify_score,
    std::pair<Eigen::Vector3d, Eigen::Matrix3d> &relative_pose,
    std::vector<std::pair<STD, STD>> &sucess_match_list, pcl::PointCloud<pcl::PointXYZINormal>::Ptr pl_cur) 
{
  sucess_match_list.clear();
  double dis_threshold = 3;
  // std::time_t solve_time = 0;
  // std::time_t verify_time = 0;
  int skip_len = (int)(candidate_matcher.match_list_.size() / 50) + 1;
  int use_size = candidate_matcher.match_list_.size() / skip_len;
  // std::vector<size_t> index(use_size);
  // std::vector<int> vote_list(use_size);
  int max_vote_index = 0;
  int max_vote = 0;
  Eigen::Matrix3d best_rot;
  Eigen::Vector3d best_t;
  // for (size_t i = 0; i < index.size(); i++) {
  //   index[i] = i;
  // }
  // std::mutex mylock;
  // auto t0 = std::chrono::high_resolution_clock::now();
  // std::for_each(
  //     std::execution::par_unseq, index.begin(), index.end(),
  //     [&](const size_t &i) 
  // for(const size_t &i: index)
  for (size_t i = 0; i < use_size; i++)
  {
    auto &single_pair = candidate_matcher.match_list_[i * skip_len];
    int vote = 0;
    Eigen::Matrix3d test_rot;
    Eigen::Vector3d test_t;
    triangle_solver(single_pair, test_t, test_rot);
    for (size_t j = 0; j < candidate_matcher.match_list_.size(); j++) 
    {
      auto &verify_pair = candidate_matcher.match_list_[j];
      // Eigen::Vector3d A = verify_pair.first.binary_A_.location_;
      Eigen::Vector3d A_transform = test_rot * verify_pair.first.binary_A_.location_ + test_t;
      if((A_transform - verify_pair.second.binary_A_.location_).norm() >= dis_threshold)
        continue;
      // Eigen::Vector3d B = verify_pair.first.binary_B_.location_;
      Eigen::Vector3d B_transform = test_rot * verify_pair.first.binary_B_.location_ + test_t;
      if((B_transform - verify_pair.second.binary_B_.location_).norm() >= dis_threshold)
        continue;
      // Eigen::Vector3d C = verify_pair.first.binary_C_.location_;
      Eigen::Vector3d C_transform = test_rot * verify_pair.first.binary_C_.location_ + test_t;
      if((C_transform - verify_pair.second.binary_C_.location_).norm() < dis_threshold)
        vote++;
      // double dis_A =
      //     (A_transform - verify_pair.second.binary_A_.location_).norm();
      // double dis_B =
      //     (B_transform - verify_pair.second.binary_B_.location_).norm();
      // double dis_C =
      //     (C_transform - verify_pair.second.binary_C_.location_).norm();
      // if (dis_A < dis_threshold && dis_B < dis_threshold &&
      //     dis_C < dis_threshold) {
      //   vote++;
      // }
    }
    // mylock.lock();
    // vote_list[i] = vote;
    if (max_vote < vote) {
      max_vote_index = i;
      max_vote = vote;
      best_rot = test_rot;
      best_t = test_t;
    }
    // mylock.unlock();
  }
  // );

  // for (size_t i = 0; i < vote_list.size(); i++) {
  //   if (max_vote < vote_list[i]) {
  //     max_vote_index = i;
  //     max_vote = vote_list[i];
  //   }
  // }
  // printf("maxvote: %d %d\n", max_vote_index, max_vote);

  // old 4
  if (max_vote >= 4) 
  {
    // auto best_pair = candidate_matcher.match_list_[max_vote_index * skip_len];
    // int vote = 0;
    // Eigen::Matrix3d best_rot;
    // Eigen::Vector3d best_t;
    // triangle_solver(best_pair, best_t, best_rot);
    relative_pose.first = best_t;
    relative_pose.second = best_rot;
    // for (size_t j = 0; j < candidate_matcher.match_list_.size(); j++) {
    //   auto verify_pair = candidate_matcher.match_list_[j];
    //   Eigen::Vector3d A = verify_pair.first.binary_A_.location_;
    //   Eigen::Vector3d A_transform = best_rot * A + best_t;
    //   Eigen::Vector3d B = verify_pair.first.binary_B_.location_;
    //   Eigen::Vector3d B_transform = best_rot * B + best_t;
    //   Eigen::Vector3d C = verify_pair.first.binary_C_.location_;
    //   Eigen::Vector3d C_transform = best_rot * C + best_t;
    //   double dis_A =
    //       (A_transform - verify_pair.second.binary_A_.location_).norm();
    //   double dis_B =
    //       (B_transform - verify_pair.second.binary_B_.location_).norm();
    //   double dis_C =
    //       (C_transform - verify_pair.second.binary_C_.location_).norm();
    //   if (dis_A < dis_threshold && dis_B < dis_threshold &&
    //       dis_C < dis_threshold) {
    //     sucess_match_list.push_back(verify_pair);
    //   }
    // }
    verify_score = plane_geometric_verify(
        pl_cur,
        plane_cloud_vec_[candidate_matcher.match_id_.second], relative_pose);
  } else {
    verify_score = -1;
  }
  return;
}

void STDescManager::triangle_solver(std::pair<STD, STD> &std_pair,
                                    Eigen::Vector3d &t, Eigen::Matrix3d &rot) {
  Eigen::Matrix3d src = Eigen::Matrix3d::Zero();
  Eigen::Matrix3d ref = Eigen::Matrix3d::Zero();
  src.col(0) = std_pair.first.binary_A_.location_ - std_pair.first.center_;
  src.col(1) = std_pair.first.binary_B_.location_ - std_pair.first.center_;
  src.col(2) = std_pair.first.binary_C_.location_ - std_pair.first.center_;
  ref.col(0) = std_pair.second.binary_A_.location_ - std_pair.second.center_;
  ref.col(1) = std_pair.second.binary_B_.location_ - std_pair.second.center_;
  ref.col(2) = std_pair.second.binary_C_.location_ - std_pair.second.center_;
  Eigen::Matrix3d covariance = src * ref.transpose();
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(covariance, Eigen::ComputeThinU |
                                                        Eigen::ComputeThinV);
  Eigen::Matrix3d V = svd.matrixV();
  Eigen::Matrix3d U = svd.matrixU();
  rot = V * U.transpose();
  if (rot.determinant() < 0) {
    Eigen::Matrix3d K;
    K << 1, 0, 0, 0, 1, 0, 0, 0, -1;
    rot = V * K * U.transpose();
  }
  t = -rot * std_pair.first.center_ + std_pair.second.center_;
}

double STDescManager::plane_geometric_verify(
    const pcl::PointCloud<pcl::PointXYZINormal>::Ptr &source_cloud,
    const pcl::PointCloud<pcl::PointXYZINormal>::Ptr &target_cloud,
    const std::pair<Eigen::Vector3d, Eigen::Matrix3d> &transform) {
  Eigen::Vector3d t = transform.first;
  Eigen::Matrix3d rot = transform.second;
  pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr kd_tree(
      new pcl::KdTreeFLANN<pcl::PointXYZ>);
  pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud(
      new pcl::PointCloud<pcl::PointXYZ>);
  for (size_t i = 0; i < target_cloud->size(); i++) {
    pcl::PointXYZ pi;
    pi.x = target_cloud->points[i].x;
    pi.y = target_cloud->points[i].y;
    pi.z = target_cloud->points[i].z;
    input_cloud->push_back(pi);
  }

  kd_tree->setInputCloud(input_cloud);
  // 创建两个向量，分别存放近邻的索引值、近邻的中心距
  std::vector<int> pointIdxNKNSearch(1);
  std::vector<float> pointNKNSquaredDistance(1);
  double useful_match = 0;
  double normal_threshold = config_setting_.normal_threshold_;
  double dis_threshold = config_setting_.dis_threshold_;
  for (size_t i = 0; i < source_cloud->size(); i++) {
    pcl::PointXYZINormal searchPoint = source_cloud->points[i];
    pcl::PointXYZ use_search_point;
    use_search_point.x = searchPoint.x;
    use_search_point.y = searchPoint.y;
    use_search_point.z = searchPoint.z;
    Eigen::Vector3d pi(searchPoint.x, searchPoint.y, searchPoint.z);
    pi = rot * pi + t;
    use_search_point.x = pi[0];
    use_search_point.y = pi[1];
    use_search_point.z = pi[2];
    Eigen::Vector3d ni(searchPoint.normal_x, searchPoint.normal_y,
                       searchPoint.normal_z);
    ni = rot * ni;
    if (kd_tree->nearestKSearch(use_search_point, 1, pointIdxNKNSearch,
                                pointNKNSquaredDistance) > 0) {
      pcl::PointXYZINormal nearstPoint =
          target_cloud->points[pointIdxNKNSearch[0]];
      Eigen::Vector3d tpi(nearstPoint.x, nearstPoint.y, nearstPoint.z);
      Eigen::Vector3d tni(nearstPoint.normal_x, nearstPoint.normal_y,
                          nearstPoint.normal_z);
      Eigen::Vector3d normal_inc = ni - tni;
      Eigen::Vector3d normal_add = ni + tni;
      double point_to_plane = fabs(tni.transpose() * (pi - tpi));
      if ((normal_inc.norm() < normal_threshold ||
           normal_add.norm() < normal_threshold) &&
          point_to_plane < dis_threshold) {
        useful_match++;
      }
    }
  }
  return useful_match / source_cloud->size();
}
