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

// associated header
#include "Cost.h"

// includes
// RBDyn
#include <RBDyn/MultiBody.h>

namespace tpg
{

CollisionCost::CollisionCost(const rbd::MultiBody& mb, int nrWp, const ObsPen& op)
  : obsPen_(op)
  , spheres_(mb.nrBodies())
  , cost_(0.)
  , grad_(mb.nrParams()*nrWp)
{
}


void CollisionCost::update(const OptimizerData& data)
{
  cost_ = 0.;
  grad_.setZero();
}


double CollisionCost::cost()
{
  return cost_;
}


const Eigen::VectorXd& CollisionCost::costGrad()
{
  return grad_;
}

} // tpg
