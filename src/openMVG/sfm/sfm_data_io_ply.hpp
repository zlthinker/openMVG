
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_IO_PLY_HPP
#define OPENMVG_SFM_DATA_IO_PLY_HPP

#include "openMVG/sfm/sfm_data_io.hpp"
#include <fstream>

namespace openMVG {
namespace sfm {

/// Save the structure and camera positions of a SfM_Data container as 3D points in a PLY ASCII file.
static bool Save_PLY(
  const SfM_Data & sfm_data,
  const std::string & filename,
  ESfM_Data flags_part)
{
  const bool b_structure = (flags_part & STRUCTURE) == STRUCTURE;
  const bool b_extrinsics = (flags_part & EXTRINSICS) == EXTRINSICS;

  if (!(b_structure || b_extrinsics))
    return false;

  //Create the stream and check it is ok
  std::ofstream stream(filename.c_str());
  if (!stream.is_open())
    return false;

  bool bOk = false;
  {
    // Count how many views having valid poses:
    IndexT view_with_pose_count = 0;
    if (b_extrinsics)
    {
      for (const auto & view : sfm_data.GetViews())
      {
        view_with_pose_count += sfm_data.IsPoseAndIntrinsicDefined(view.second.get());
      }
    }
    stream << "ply"
      << '\n' << "format ascii 1.0"
      << '\n' << "element vertex "
        // Vertex count: (#landmark + #view_with_valid_pose)
        << ((b_structure ? sfm_data.GetLandmarks().size() : 0) +
            view_with_pose_count * 5)
      << '\n' << "property float x"
      << '\n' << "property float y"
      << '\n' << "property float z"
      << '\n' << "property uchar red"
      << '\n' << "property uchar green"
      << '\n' << "property uchar blue"
      << '\n' << "element edge " << view_with_pose_count * 8
      << '\n' << "property int vertex1"
      << '\n' << "property int vertex2"
      << '\n' << "property uchar red"
      << '\n' << "property uchar green"
      << '\n' << "property uchar blue"
      << '\n' << "end_header" << std::endl;

      if (b_extrinsics)
      {
        for (const auto & view : sfm_data.GetViews())
        {
          if (sfm_data.IsPoseAndIntrinsicDefined(view.second.get()))
          {
            const geometry::Pose3 pose = sfm_data.GetPoseOrDie(view.second.get());
            Mat3 rotation = pose.rotation();
            Vec3 center = pose.center();
            float len = 0.5;
            Vec3 X1 = Vec3(len, len, len);
            Vec3 X2 = Vec3(len, -len, len);
            Vec3 X3 = Vec3(-len, -len, len);
            Vec3 X4 = Vec3(-len, len, len);
            Vec3 X1_trans = rotation.transpose() * X1 + center;
            Vec3 X2_trans = rotation.transpose() * X2 + center;
            Vec3 X3_trans = rotation.transpose() * X3 + center;
            Vec3 X4_trans = rotation.transpose() * X4 + center;
            stream << pose.center().transpose()
              << " 0 255 0" << "\n"
            << X1_trans.transpose() << " 0 255 0" << "\n"
            << X2_trans.transpose() << " 0 255 0" << "\n"
            << X3_trans.transpose() << " 0 255 0" << "\n"
            << X4_trans.transpose() << " 0 255 0" << "\n";
          }
        }
      }

      if (b_structure)
      {
        const Landmarks & landmarks = sfm_data.GetLandmarks();
        for (Landmarks::const_iterator iterLandmarks = landmarks.begin();
          iterLandmarks != landmarks.end();
          ++iterLandmarks)  {
          stream << iterLandmarks->second.X.transpose() << " " << iterLandmarks->second.RGB.transpose() << "\n";
        }
      }

      for(int i = 0; i < view_with_pose_count; i++) {
        stream << 5 * i << " " << 5 * i + 1 << " 255 255 255" << std::endl
        << 5 * i << " " << 5 * i + 2 << " 255 255 255" << std::endl
        << 5 * i << " " << 5 * i + 3 << " 255 255 255" << std::endl
        << 5 * i << " " << 5 * i + 4 << " 255 255 255" << std::endl
        << 5 * i + 1 << " " << 5 * i + 2 << " 255 255 255" << std::endl
        << 5 * i + 2 << " " << 5 * i + 3 << " 255 255 255" << std::endl
        << 5 * i + 3 << " " << 5 * i + 4 << " 255 255 255" << std::endl
        << 5 * i + 4 << " " << 5 * i + 1 << " 255 255 255" << std::endl;
      }
      stream.flush();
      bOk = stream.good();
      stream.close();
  }
  return bOk;
}

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_DATA_IO_PLY_HPP
