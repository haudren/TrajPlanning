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

// include
// Eigen
#include <Eigen/Core>

// sva
#include <SpaceVecAlg/SpaceVecAlg>

// forward declaration
namespace rbd
{
class MultiBody;
}

namespace tpg
{
class ObsPen;
class OptimizerData;


class CollisionCost
{
public:
  CollisionCost(const rbd::MultiBody& mb, int nrWp, const ObsPen& obsPen);

  void update(const OptimizerData& data);

  double cost();
  const Eigen::VectorXd& costGrad();

private:
  struct SphereData
  {
    int bodyIndex;
    double radius;
    sva::PTransformd pos;
  };

private:
  const ObsPen& obsPen_;
  std::vector<SphereData> spheres_;
  double cost_;
  Eigen::VectorXd grad_;
};


} // tpg
