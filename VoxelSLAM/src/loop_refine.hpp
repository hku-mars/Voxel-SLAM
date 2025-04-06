#ifndef LOOP_REFINE_HPP
#define LOOP_REFINE_HPP

#include "tools.hpp"
#include "voxel_map.hpp"

// #include "STDesc.h"
#include <gtsam/geometry/Pose3.h>
#include <gtsam/slam/PriorFactor.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/nonlinear/ISAM2.h>
#include <pcl/kdtree/kdtree_flann.h>

using namespace std;

struct ScanPose
{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  IMUST x;
  PVecPtr pvec;
  Eigen::Matrix<double, 6, 1> v6;

  ScanPose(IMUST &_x, PVecPtr _pvec): x(_x), pvec(_pvec)
  {
    v6.setZero();
  }

  void update(IMUST dx)
  {
    x.v = dx.R * x.v;
    x.p = dx.R * x.p + dx.p;
    x.R = dx.R * x.R;
  }

  void set_state(const gtsam::Pose3 &pose)
  {
    Eigen::Matrix3d rot = pose.rotation().matrix();
    rot = rot * x.R.transpose();
    x.R = pose.rotation().matrix();
    x.p = pose.translation();
    x.v = rot * x.v;
  }

};

bool icp_normal(pcl::PointCloud<PointType> &pl_src, pcl::PointCloud<PointType> &pl_tar, pair<Eigen::Vector3d, Eigen::Matrix3d> &pose, double icp_eigval)
{
  pcl::KdTreeFLANN<pcl::PointXYZ> kd_tree;
  pcl::PointCloud<pcl::PointXYZ> input_cloud;
  for(PointType &ap: pl_tar.points)
  {
    pcl::PointXYZ pi;
    pi.x = ap.x; pi.y = ap.y; pi.z = ap.z;
    input_cloud.push_back(pi);
  }
  kd_tree.setInputCloud(input_cloud.makeShared());

  vector<int> pointIdxNKNSearch(1);
  vector<float> pointNKNSquaredDistance(1);

  Eigen::Vector4d paras(0.2, 0.2, 0.5, 3);
  int is_converge = 0;

  int ssize = pl_src.size();
  int match_num = 0;
  Eigen::Matrix3d mat_norm;
  // for(int iterCount=0; iterCount<10; iterCount++)
  for(int iterCount=0; iterCount<20; iterCount++)
  {
    Eigen::Matrix<double, 6, 6> Hess; Hess.setZero();
    Eigen::Matrix<double, 6, 1> JacT; JacT.setZero();
    double resi = 0;
    match_num = 0;
    mat_norm.setZero();

    for(int i=0; i<ssize; i++)
    {
      PointType &searchPoint = pl_src[i];
      Eigen::Vector3d plocal(searchPoint.x, searchPoint.y, searchPoint.z);
      Eigen::Vector3d pi = pose.second * plocal + pose.first;
      pcl::PointXYZ use_search_point;
      use_search_point.x = pi[0];
      use_search_point.y = pi[1];
      use_search_point.z = pi[2];
      Eigen::Vector3d ni(searchPoint.normal_x, searchPoint.normal_y, searchPoint.normal_z);
      ni = pose.second * ni;
      if (kd_tree.nearestKSearch(use_search_point, 1, pointIdxNKNSearch,pointNKNSquaredDistance) > 0) 
      {
        pcl::PointXYZINormal nearstPoint = pl_tar[pointIdxNKNSearch[0]];
        Eigen::Vector3d tpi(nearstPoint.x, nearstPoint.y, nearstPoint.z);
        Eigen::Vector3d tni(nearstPoint.normal_x, nearstPoint.normal_y,nearstPoint.normal_z);
        Eigen::Vector3d normal_inc = ni - tni;
        Eigen::Vector3d normal_add = ni + tni;
        double point_to_point_dis = (pi - tpi).norm();
        double point_to_plane = fabs(tni.transpose() * (pi - tpi));
        if ((normal_inc.norm() < paras[0] || normal_add.norm() < paras[1]) && point_to_plane < paras[2] && point_to_point_dis < paras[3])
        // if ((normal_inc.norm() < 0.1 || normal_add.norm() < 0.1) &&   point_to_plane < 0.1 && point_to_point_dis < 1)
        {
          double rr = tni.dot(pi - tpi);
          Eigen::Matrix<double, 6, 1> jac;
          jac.head(3) = hat(plocal) * pose.second.transpose() * tni;
          jac.tail(3) = tni;

          Hess += jac * jac.transpose();
          JacT += jac * rr;
          resi += 0.5 * rr * rr;
          match_num++;
          mat_norm += tni * tni.transpose();
        }

      }
    }

    Eigen::Matrix<double, 6, 1> dxi = Hess.ldlt().solve(-JacT);
    pose.second = pose.second * Exp(dxi.head(3));
    pose.first = pose.first + dxi.tail(3);

    // printf("icp%d: %lf %d\n", iterCount, resi, match_num);

    if(dxi.head(3).norm()<1e-3 && dxi.tail(3).norm()<1e-3)
    {
      if(is_converge)
        break;
      else
      {
        paras << 0.1, 0.1, 0.1, 1;
        is_converge = 1;
      }
    }

    // if(is_converge == 0 && iterCount > 4)
    // {
    //   paras << 0.1, 0.1, 0.1, 1;
    //   is_converge = 1;
    // }
  }

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> saes(mat_norm);
  Eigen::Vector3d eig_vec = saes.eigenvalues();
  printf("eigvalue: %lf %lf %lf %d\n", eig_vec[0], eig_vec[1], eig_vec[2], is_converge);
  return eig_vec[0] > icp_eigval && is_converge == 1;
  // return eig_vec[0] > icp_eigval;

}

void add_edge(int pos1, int pos2, IMUST &x1, IMUST &x2, gtsam::NonlinearFactorGraph &graph, gtsam::noiseModel::Diagonal::shared_ptr odometryNoise)
{
  gtsam::Point3 tt(x1.R.transpose() * (x2.p - x1.p));
  gtsam::Rot3 RR(x1.R.transpose() * x2.R);
  gtsam::NonlinearFactor::shared_ptr factor(new gtsam::BetweenFactor<gtsam::Pose3>(pos1, pos2, gtsam::Pose3(RR, tt), odometryNoise));
  graph.push_back(factor);
}

void add_edge(int pos1, int pos2, Eigen::Matrix3d &rot, Eigen::Vector3d &tra, gtsam::NonlinearFactorGraph &graph, gtsam::noiseModel::Diagonal::shared_ptr odometryNoise)
{
  gtsam::Point3 tt(tra);
  gtsam::Rot3 RR(rot);
  gtsam::NonlinearFactor::shared_ptr factor(new gtsam::BetweenFactor<gtsam::Pose3>(pos1, pos2, gtsam::Pose3(RR, tt), odometryNoise));
  graph.push_back(factor);
}

struct PGO_Edge
{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  int m1, m2;
  vector<int> ids1, ids2;
  PLM(3) rots;
  PLV(3) tras;
  PLV(6) covs;

  PGO_Edge(int _m1, int _m2, int id1, int id2, Eigen::Matrix3d &rot, Eigen::Vector3d &tra, Eigen::Matrix<double, 6, 1> &v6): m1(_m1), m2(_m2) 
  {
    push(id1, id2, rot, tra, v6);
  }

  void push(int id1, int id2, Eigen::Matrix3d &rot, Eigen::Vector3d &tra, Eigen::Matrix<double, 6, 1> &v6)
  {
    ids1.push_back(id1); ids2.push_back(id2);
    rots.push_back(rot); tras.push_back(tra);
    covs.push_back(v6);
  }

  bool is_adapt(vector<int> &maps, vector<int> &step)
  {
    bool f1 = false, f2 = false;
    for(int i=0; i<maps.size(); i++)
    {
      if(m1 == maps[i])
      {
        f1 = true;
        step[0] = i;
      }
      if(m2 == maps[i])
      {
        f2 = true;
        step[1] = i;
      }
    }

    return f1 && f2;
  }

};

struct PGO_Edges
{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  vector<PGO_Edge> edges;
  vector<vector<int>> mates;

  void push(int _m1, int _m2, int _id1, int _id2, Eigen::Matrix3d &rot, Eigen::Vector3d &tra, Eigen::Matrix<double, 6, 1> &v6)
  {
    bool is_lack = true;
    for(PGO_Edge &e: edges)
    {
      if(e.m1 == _m1 && e.m2 == _m2)
      {
        is_lack = false;
        e.push(_id1, _id2, rot, tra, v6);
        break;
      }
    }

    if(is_lack)
    {
      edges.emplace_back(_m1, _m2, _id1, _id2, rot, tra, v6);
      int msize = mates.size();
      for(int i=msize; i<_m2+1; i++)
        mates.emplace_back();
      mates[_m1].push_back(_m2);
      mates[_m2].push_back(_m1);
    }
      
  }

  void connect(int root, vector<int> &ids)
  {
    ids.clear();
    ids.push_back(root);
    tras(root, ids);
    sort(ids.begin(), ids.end());
  }

  void tras(int ord, vector<int> &ids)
  {
    if(ord < mates.size())
    for(int id: mates[ord])
    {
      bool is_exist = false;
      for(int i: ids)
      if(id == i)
      {
        is_exist = true;
        break;
      }

      if(!is_exist)
      {
        ids.push_back(id);
        tras(id, ids);
      }
    }

  }

};

vector<double> gba_eigen_value_array;
double gba_min_eigen_value;
double gba_voxel_size;

class OctreeGBA
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  vector<PLV(3)> locals, worlds;
  PointCluster pcr_add;
  int layer, octo_state, wdsize;
  OctreeGBA* leaves[8];
  double voxel_center[3];
  float quater_length;
  bool is_plane;

  float ref;

  OctreeGBA(int _l, int _w): layer(_l), wdsize(_w)
  {
    locals.resize(wdsize);
    worlds.resize(wdsize);
    for(int i=0; i<8; i++) leaves[i] = nullptr;
    is_plane = false;
    octo_state = 0;

    ref = 255.0*rand()/(RAND_MAX + 1.0f);
  }

  void reset(int _l, int _w)
  {
    layer = _l; wdsize = _w;
    locals.clear(); worlds.clear();
    locals.resize(wdsize);
    worlds.resize(wdsize);
    pcr_add.clear();
    octo_state = 0;
    is_plane = false;
    for(int i=0; i<8; i++) leaves[i] = nullptr;
  }

  inline bool plane_judge(Eigen::Vector3d &eig_values)
  {
    // return (eig_values[0] < min_eigen_value);
    return (eig_values[0] < gba_min_eigen_value && (eig_values[0]/eig_values[2])<gba_eigen_value_array[layer]);
  }

  void push(int ord, Eigen::Vector3d &local, Eigen::Vector3d &world)
  {
    locals[ord].push_back(local);
    worlds[ord].push_back(world);
    pcr_add.push(world);
  }

  void subdivide(vector<OctreeGBA*> &oct_buf)
  {
    for(int i=0; i<wdsize; i++)
    for(int j=0; j<locals[i].size(); j++)
    {
      Eigen::Vector3d &pl = locals[i][j];
      Eigen::Vector3d &pw = worlds[i][j];
      int xyz[3] = {0, 0, 0};
      for(int k=0; k<3; k++)
        if(pw[k] > voxel_center[k])
          xyz[k] = 1;
      int leafnum = 4*xyz[0] + 2*xyz[1] + xyz[2];
      if(leaves[leafnum] == nullptr)
      {
        if(oct_buf.size() > 0)
        {
          leaves[leafnum] = oct_buf.back();
          leaves[leafnum]->reset(layer+1, wdsize);
          oct_buf.pop_back();
        }
        else
        {
          leaves[leafnum] = new OctreeGBA(layer+1, wdsize);
        }
        
        leaves[leafnum]->voxel_center[0] = voxel_center[0] + (2*xyz[0]-1)*quater_length;
        leaves[leafnum]->voxel_center[1] = voxel_center[1] + (2*xyz[1]-1)*quater_length;
        leaves[leafnum]->voxel_center[2] = voxel_center[2] + (2*xyz[2]-1)*quater_length;
        leaves[leafnum]->quater_length = quater_length / 2;
      }

      leaves[leafnum]->push(i, pl, pw);
    }
  }

  void recut(LidarFactor &vox_opt, vector<OctreeGBA*> &oct_buf)
  {
    if(pcr_add.N <= 10)
      return;
    
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> saes(pcr_add.cov());
    Eigen::Vector3d eig_value = saes.eigenvalues();
    Eigen::Matrix3d eig_vector = saes.eigenvectors();
    is_plane = plane_judge(eig_value);

    if(is_plane)
    {
      if(pcr_add.N < 10) return;

      int exi = 0;
      for(int i=0; i<wdsize; i++)
        if(locals[i].size() != 0)
          exi++;
      if(exi <= 1) return;

      if(eig_value[0]/eig_value[1] > 0.12) return;

      double coe = 1.0;
      // double coe = 1.0 / pcr_add.N;
      vector<PointCluster> pcrs(wdsize);
      for(int i=0; i<wdsize; i++)
        for(Eigen::Vector3d &v: locals[i])
          pcrs[i].push(v);
      PointCluster pcr_fix;
      vox_opt.push_voxel(pcrs, pcr_fix, coe, eig_value, eig_vector, pcr_add);

      return;
    }
    else if(layer >= max_layer)
    {
      return;
    }
    else
    {
      subdivide(oct_buf);
      octo_state = 1;
    }

    for(int i=0; i<8; i++)
      if(leaves[i] != nullptr)
        leaves[i]->recut(vox_opt, oct_buf);

  }

  void tras_ptr(vector<OctreeGBA*> &oct_buf)
  {
    if(octo_state == 1)
    for(int i=0; i<8; i++)
    if(leaves[i] != nullptr)
    {
      leaves[i]->tras_ptr(oct_buf);
      oct_buf.push_back(leaves[i]);
      leaves[i] = nullptr;
    }
  }

  void tras_display(pcl::PointCloud<PointType> &pl_send)
  {
    if(octo_state == 0)
    {
      PointType ap; ap.intensity = ref;
      for(int i=0; i<wdsize; i++)
      for(Eigen::Vector3d &v: worlds[i])
      {
        ap.x = v[0]; ap.y = v[1]; ap.z = v[2];
        pl_send.push_back(ap);
      }
    }
    else
    {
      for(int i=0; i<8; i++)
        if(leaves[i] != nullptr)
          leaves[i]->tras_display(pl_send);
    }
  }

  // ~OctreeGBA()
  // {
  //   for(int i=0; i<8; i++)
  //     if(leaves[i] != nullptr)
  //       delete leaves[i];
  // }

  static void cut_voxel(unordered_map<VOXEL_LOC, OctreeGBA*> &feat_map, IMUST &xc, pcl::PointCloud<PointType>::Ptr plptr, int win_count, int wdsize)
  {
    for(PointType &ap: plptr->points)
    {
      Eigen::Vector3d local(ap.x, ap.y, ap.z);
      Eigen::Vector3d world = xc.R * local + xc.p;
      float loc[3];
      for(int j=0; j<3; j++)
      {
        loc[j] = world[j] / gba_voxel_size;
        if(loc[j] < 0) loc[j] -= 1;
      }

      VOXEL_LOC position(loc[0], loc[1], loc[2]);
      auto iter = feat_map.find(position);
      if(iter != feat_map.end())
      {
        iter->second->push(win_count, local, world);
      }
      else
      {
        OctreeGBA *ot = new OctreeGBA(0, wdsize);
        ot->push(win_count, local, world);
        ot->voxel_center[0] = (0.5+position.x) * gba_voxel_size;
        ot->voxel_center[1] = (0.5+position.y) * gba_voxel_size;
        ot->voxel_center[2] = (0.5+position.z) * gba_voxel_size;
        ot->quater_length = gba_voxel_size / 4.0;
        feat_map[position] = ot;
      }

    }


  }

};

void OctreeGBA_multi_recut(unordered_map<VOXEL_LOC, OctreeGBA*> &feat_map, LidarFactor &voxhess, int thd_num)
{
  vector<vector<OctreeGBA*>> octss(thd_num);
  vector<LidarFactor> vec_voxhess(thd_num, voxhess);
  int g_size = feat_map.size();
  vector<thread*> mthreads(thd_num);
  double part = 1.0 * g_size / thd_num;
  int cnt = 0;
  for(auto iter=feat_map.begin(); iter!=feat_map.end(); iter++)
  {
    octss[cnt].push_back(iter->second);
    if(octss[cnt].size() >= part && cnt < thd_num-1)
      cnt++;
  }
  
  auto recut_func = [](vector<OctreeGBA*> &octs, LidarFactor &voxhess)
  {
    vector<OctreeGBA*> oct_buf;
    for(OctreeGBA *oc: octs)
    {
      oc->recut(voxhess, oct_buf);
      oc->tras_ptr(oct_buf);
      delete oc;
      for(OctreeGBA *oct: oct_buf)
        delete oct;
      oct_buf.clear();
    }

    // for(OctreeGBA *oc: oct_buf)
    //   delete oc;
  };
  
  for(int i=1; i<thd_num; i++)
    mthreads[i] = new thread(recut_func, ref(octss[i]), ref(vec_voxhess[i]));
  
  for(int i=0; i<thd_num; i++)
  {
    if(i == 0)
    {
      recut_func(octss[0], voxhess);
    }
    else
    {
      mthreads[i]->join();
      delete mthreads[i];

      voxhess.plvec_voxels.insert(voxhess.plvec_voxels.end(), vec_voxhess[i].plvec_voxels.begin(), vec_voxhess[i].plvec_voxels.end());
      voxhess.sig_vecs.insert(voxhess.sig_vecs.end(), vec_voxhess[i].sig_vecs.begin(), vec_voxhess[i].sig_vecs.end());
      voxhess.coeffs.insert(voxhess.coeffs.end(), vec_voxhess[i].coeffs.begin(), vec_voxhess[i].coeffs.end());
      voxhess.eig_values.insert(voxhess.eig_values.end(), vec_voxhess[i].eig_values.begin(), vec_voxhess[i].eig_values.end());
      voxhess.eig_vectors.insert(voxhess.eig_vectors.end(), vec_voxhess[i].eig_vectors.begin(), vec_voxhess[i].eig_vectors.end());
      voxhess.pcr_adds.insert(voxhess.pcr_adds.end(), vec_voxhess[i].pcr_adds.begin(), vec_voxhess[i].pcr_adds.end());
    }
  }
}

#endif
