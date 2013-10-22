// This file is part of TrajPlanning.
//
// TrajPlanning is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// TrajPlanning is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with TrajPlanning.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

// includes
// Boost
#include <boost/multi_array.hpp>

// Eigen
#include <Eigen/Core>

namespace tpg
{

class ObsPen
{
public:
  ObsPen& operator=(const ObsPen& op);

  void setPen(const Eigen::Vector3d& start, const Eigen::Vector3d& scale,
      int sizeX, int sizeY, int sizeZ,
      const std::vector<double>& penality,
      const std::vector<double>& penalityGradX,
      const std::vector<double>& penalityGradY,
      const std::vector<double>& penalityGradZ);

  double penality(const Eigen::Vector3d& pos);
  Eigen::Vector3d penalityGrad(const Eigen::Vector3d& pos);

private:
  boost::multi_array<double, 3> pen_;
  boost::multi_array<Eigen::Vector3d, 3> penGrad_;
  Eigen::Vector3d start_;
  Eigen::Array3d scale_;
};

} // tpg
