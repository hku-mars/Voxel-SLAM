/*
 * Copyright (c) 2008, Willow Garage, Inc.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Willow Garage, Inc. nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef VOXELSLAM_PC2_H
#define VOXELSLAM_PC2_H

#include <sensor_msgs/PointCloud2.h>

#include <rviz/message_filter_display.h>

namespace rviz
{
class IntProperty;
class PointCloudCommon;
}

/**
 * \class PointCloud2Display
 * \brief Displays a point cloud of type sensor_msgs::PointCloud2
 *
 * By default it will assume channel 0 of the cloud is an intensity value, and will color them by
 * intensity.
 * If you set the channel's name to "rgb", it will interpret the channel as an integer rgb value, with r,
 * g and b
 * all being 8 bits.
 */
namespace voxelslam_pointcloud2
{

class PointCloud2Display : public rviz::MessageFilterDisplay<sensor_msgs::PointCloud2>
{
  Q_OBJECT
public:
  PointCloud2Display();
  ~PointCloud2Display() override;

  void reset() override;

  void update(float wall_dt, float ros_dt) override;

protected:
  /** @brief Do initialization. Overridden from MessageFilterDisplay. */
  void onInitialize() override;

  /** @brief Process a single message.  Overridden from MessageFilterDisplay. */
  void processMessage(const sensor_msgs::PointCloud2ConstPtr& cloud) override;

  rviz::PointCloudCommon* point_cloud_common_;
};

} // namespace rviz

#endif
