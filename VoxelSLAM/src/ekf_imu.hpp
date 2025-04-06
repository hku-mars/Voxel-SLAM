#ifndef EKF_IMU_HPP
#define EKF_IMU_HPP

#include "tools.hpp"
#include <deque>
#include <sensor_msgs/Imu.h>

class IMUEKF
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  bool init_flag;
  double pcl_beg_time, pcl_end_time, last_pcl_end_time;
  int init_num;
  Eigen::Vector3d mean_acc, mean_gyr;
  sensor_msgs::Imu::Ptr last_imu;
  int min_init_num = 30;
  Eigen::Vector3d angvel_last, acc_s_last;

  Eigen::Vector3d cov_acc, cov_gyr;
  Eigen::Vector3d cov_bias_gyr, cov_bias_acc;
  // Eigen::Matrix<double, DIM, DIM> cov, cov_inv;

  Eigen::Matrix3d Lid_rot_to_IMU;
  Eigen::Vector3d Lid_offset_to_IMU;

  double scale_gravity = 1.0;
  vector<IMUST> imu_poses;
  string imu_topic = "/livox/imu";

  int point_notime = 0;

  IMUEKF()
  {
    init_flag = false;
    init_num = 0;
    mean_acc.setZero(); mean_gyr.setZero();
    angvel_last.setZero(); acc_s_last.setZero();
  }

  void motion_blur(IMUST &xc, pcl::PointCloud<PointType> &pcl_in, deque<sensor_msgs::Imu::Ptr> &imus)
  {
    imus.push_front(last_imu);

    if(last_pcl_end_time - pcl_beg_time > 0.01)
    {
      printf("%lf %lf\n", pcl_beg_time, last_pcl_end_time);
      printf("LiDAR time regress. Please check data\n"); exit(0);
    }

    imu_poses.clear();
    // imu_poses.emplace_back(0, xc.R, xc.p, xc.v, angvel_last, acc_s_last);

    Eigen::Vector3d acc_imu, angvel_avr, acc_avr, vel_imu(xc.v), pos_imu(xc.p);
    Eigen::Matrix3d R_imu(xc.R);
    Eigen::Matrix<double, DIM, DIM> F_x, cov_w;

    double dt = 0;
    for(auto it_imu=imus.begin(); it_imu!=imus.end()-1; it_imu++)
    {
      sensor_msgs::Imu &head = **(it_imu);
      sensor_msgs::Imu &tail = **(it_imu+1);

      if(tail.header.stamp.toSec() < last_pcl_end_time) continue;

      angvel_avr << 0.5*(head.angular_velocity.x + tail.angular_velocity.x), 
                    0.5*(head.angular_velocity.y + tail.angular_velocity.y), 
                    0.5*(head.angular_velocity.z + tail.angular_velocity.z);
      acc_avr << 0.5*(head.linear_acceleration.x + tail.linear_acceleration.x), 
                 0.5*(head.linear_acceleration.y + tail.linear_acceleration.y), 
                 0.5*(head.linear_acceleration.z + tail.linear_acceleration.z);

      angvel_avr -= xc.bg;
      acc_avr = acc_avr * scale_gravity - xc.ba;
      acc_imu = R_imu * acc_avr + xc.g;

      // if(head.header.stamp.toSec() < last_pcl_end_time)
      //   dt = tail.header.stamp.toSec() - last_pcl_end_time;
      // else
      //   dt = tail.header.stamp.toSec() - head.header.stamp.toSec();
      double cur_time = head.header.stamp.toSec();
      if(cur_time < last_pcl_end_time)
        cur_time = last_pcl_end_time;
      dt = tail.header.stamp.toSec() - cur_time;

      double offt = cur_time - pcl_beg_time;
      imu_poses.emplace_back(offt, R_imu, pos_imu, vel_imu, angvel_avr, acc_imu);
      
      Eigen::Matrix3d acc_avr_skew = hat(acc_avr);
      Eigen::Matrix3d Exp_f = Exp(angvel_avr, dt);

      F_x.setIdentity();
      cov_w.setZero();

      F_x.block<3,3>(0,0)  = Exp(angvel_avr, - dt);
      F_x.block<3,3>(0,9)  = -I33 * dt;
      F_x.block<3,3>(3,6)  = I33 * dt;
      F_x.block<3,3>(6,0)  = - R_imu * acc_avr_skew * dt;
      F_x.block<3,3>(6,12) = - R_imu * dt;
      cov_w.block<3,3>(0,0).diagonal() = cov_gyr * dt * dt;
      cov_w.block<3,3>(6,6) = R_imu * cov_acc.asDiagonal() * R_imu.transpose() * dt * dt;
      cov_w.block<3,3>(9,9).diagonal()   = cov_bias_gyr * dt * dt;
      cov_w.block<3,3>(12,12).diagonal() = cov_bias_acc * dt * dt;

      xc.cov = F_x * xc.cov * F_x.transpose() + cov_w;
      // R_imu = R_imu * Exp_f;
      // acc_imu = R_imu * acc_avr + xc.g;
      pos_imu = pos_imu + vel_imu * dt + 0.5 * acc_imu * dt * dt;
      vel_imu = vel_imu + acc_imu * dt;
      R_imu = R_imu * Exp_f;

      // double offt = tail.header.stamp.toSec() - pcl_beg_time;
      // double offt = max(head.header.stamp.toSec() - pcl_beg_time, 0.0);
      // imu_poses.emplace_back(offt, R_imu, pos_imu, vel_imu, angvel_avr, acc_imu);
    }

    double imu_end_time = imus.back()->header.stamp.toSec();
    double note = pcl_end_time > imu_end_time ? 1.0 : -1.0;
    dt = note * (pcl_end_time - imu_end_time);
    xc.v = vel_imu + note * acc_imu * dt;
    xc.R = R_imu * Exp(note*angvel_avr, dt);
    xc.p = pos_imu + note * vel_imu * dt + note * 0.5 * acc_imu * dt * dt;
    xc.t = pcl_end_time;

    sensor_msgs::ImuPtr imu1(new sensor_msgs::Imu(*imus.front()));
    sensor_msgs::ImuPtr imu2(new sensor_msgs::Imu(*imus.back()));
    imu1->header.stamp.fromSec(last_pcl_end_time);
    imu2->header.stamp.fromSec(pcl_end_time);
    // imus.pop_front();
    last_imu = imus.back();
    last_pcl_end_time = pcl_end_time;
    imus.front() = imu1;
    imus.back()  = imu2;

    if(point_notime)
      return;

    auto it_pcl = pcl_in.end() - 1;
    for(int i=imu_poses.size()-1; i>=0; i--)
    {
      IMUST &head = imu_poses[i];
      R_imu = head.R;
      acc_imu = head.ba;
      vel_imu = head.v;
      pos_imu = head.p;
      angvel_avr = head.bg;

      for(; it_pcl->curvature > head.t; it_pcl--)
      {
        dt = it_pcl->curvature - head.t;

        Eigen::Matrix3d R_i = R_imu * Exp(angvel_avr, dt);
        Eigen::Vector3d T_ei = pos_imu + vel_imu * dt + 0.5 * acc_imu * dt * dt - xc.p;

        Eigen::Vector3d P_i(it_pcl->x, it_pcl->y, it_pcl->z);
        Eigen::Vector3d P_compensate = Lid_rot_to_IMU.transpose() * (xc.R.transpose() * (R_i * (Lid_rot_to_IMU * P_i + Lid_offset_to_IMU) + T_ei) - Lid_offset_to_IMU);

        it_pcl->x = P_compensate(0);
        it_pcl->y = P_compensate(1);
        it_pcl->z = P_compensate(2);
        if(it_pcl == pcl_in.begin()) break;
      }
    }

  }

  void IMU_init(deque<sensor_msgs::Imu::Ptr> &imus)
  {
    Eigen::Vector3d cur_acc, cur_gyr;
    for(sensor_msgs::Imu::Ptr imu: imus)
    {
      cur_acc << imu->linear_acceleration.x,
                 imu->linear_acceleration.y,
                 imu->linear_acceleration.z;
      cur_gyr << imu->angular_velocity.x,
                 imu->angular_velocity.y,
                 imu->angular_velocity.z;

      if(init_num != 0)
      {
        mean_acc += (cur_acc - mean_acc) / init_num; 
        mean_gyr += (cur_gyr - mean_gyr) / init_num;
      }
      else
      {
        mean_acc = cur_acc; // modify
        mean_gyr = cur_gyr;
        init_num = 1;
      }
      
      init_num++;
    }

    last_imu = imus.back();
  }

  int process(IMUST &x_curr, pcl::PointCloud<PointType> &pcl_in, deque<sensor_msgs::Imu::Ptr> &imus)
  {
    if(!init_flag)
    {
      IMU_init(imus);
      if(mean_acc.norm() < 2 && imu_topic == "/livox/imu")
        scale_gravity = G_m_s2;
      printf("scale_gravity: %lf %lf %d\n", scale_gravity, mean_acc.norm(), init_num);
      x_curr.g = -mean_acc * scale_gravity;
      if(init_num > min_init_num) init_flag = true;
      last_pcl_end_time = pcl_end_time;

      return 0;
    }

    motion_blur(x_curr, pcl_in, imus);
    return 1;
  }

};

#endif

