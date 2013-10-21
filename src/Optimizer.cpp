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
#include "Optimizer.h"

// includes
// RBDyn
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>

// TrajPlanning
#include "SmoothnessMatrix.h"

namespace tpg
{

void Optimizer::init(const OptimizerConfig& oc)
{
  collisionSphere_.clear();
  wpData_.clear();

  mb_ = oc.mb;
  mbc_ = rbd::MultiBodyConfig(mb_);
  mbc_.zero(mb_);
  rbd::forwardKinematics(mb_, mbc_);
  rbd::forwardVelocity(mb_, mbc_);

  // collision data
  pen_ = oc.pen;
  for(const auto& s: oc.collisionSpheres)
  {
    collisionSphere_.push_back({mb_.bodyIndexById(s.bodyId), s.radius,
                                s.position});
  }

  // WP data
  WPData wpdata;
  wpdata.mbc = mbc_;
  wpdata.bodyVel.resize(mb_.nrBodies());
  wpdata.bodyVelNorm.resize(mb_.nrBodies());
  for(const rbd::Body& b: mb_.bodies())
  {
    wpdata.jacVel.emplace_back(mb_, b.id());
  }
  for(const Sphere& s: oc.collisionSpheres)
  {
    wpdata.jacSphere.emplace_back(mb_, s.bodyId, s.position);
  }
  wpData_.resize(oc.nrWp, wpdata);

  // start-end poses
  // start
  rbd::MultiBodyConfig mbcStartEnd(oc.start);
  rbd::forwardKinematics(mb_, mbcStartEnd);
  startPoses_ = mbcStartEnd.bodyPosW;
  // end
  mbcStartEnd = oc.end;
  rbd::forwardKinematics(mb_, mbcStartEnd);
  endPoses_ = mbcStartEnd.bodyPosW;

  // smoothness
  auto jSM = jerkSmoothness(oc.nrWp,
                            rbd::paramToVector(mb_, oc.start.q),
                            rbd::paramToVector(mb_, oc.end.q),
                            oc.velWeight, oc.accWeight, oc.jerkWeight);
  A_ = std::get<0>(jSM);
  b_ = std::get<1>(jSM);
  c_ = std::get<2>(jSM);
  AInv_.compute(A_);

  // optim
  path_.resize(mb_.nrParams()*oc.nrWp);
  obsGrad_.resize(path_.cols());
  speedGrad_.resize(path_.cols());
}


void Optimizer::optimize(int nrIter, double learningRate, const Eigen::VectorXd& initPath)
{
  for(int iter = 0; iter < nrIter; ++iter)
  {
    computeFK();
    computeVel();
    computeJac();
  }
}


void Optimizer::computeFK()
{
  for(int wp = 0; wp < int(wpData_.size()); ++wp)
  {
    WPData& data = wpData_[wp];
    rbd::vectorToParam(path_.segment(wp*mb_.nrParams(), (wp+1)*mb_.nrParams()),
                       data.mbc.q);
    rbd::forwardKinematics(mb_, data.mbc);
  }
}


void Optimizer::computeVel()
{
  for(int wp = 0; wp < int(wpData_.size() - 1); ++wp)
  {
    WPData& data = wpData_[wp];
    WPData& dataNext = wpData_[wp + 1];
    computeVel(data, dataNext.mbc.bodyPosW);
  }
  computeVel(wpData_.back(), endPoses_);
}


void Optimizer::computeVel(WPData& data, const std::vector<sva::PTransformd>& nextPoses)
{
  for(std::size_t bi = 0; bi < data.mbc.bodyPosW.size(); ++bi)
  {
    data.bodyVel[bi] = nextPoses[bi].translation() - data.mbc.bodyPosW[bi].translation();
    data.bodyVelNorm[bi] = data.bodyVel[bi].norm();
  }
}


void Optimizer::computeJac()
{
  for(int wp = 0; wp < int(wpData_.size()); ++wp)
  {
    WPData& data = wpData_[wp];
    for(std::size_t ji = 0; ji < data.jacVel.size(); ++ji)
    {
      data.jacVel[ji].jacobian(mb_, data.mbc);
    }
    for(std::size_t js = 0; js < data.jacSphere.size(); ++js)
    {
      data.jacSphere[js].jacobian(mb_, data.mbc);
    }
  }
}


void Optimizer::computeSmCost()
{
  smCost_ = (path_.transpose()*A_*path_ + path_.transpose()*b_)[0] + c_;
}


void Optimizer::computeSpeedCost()
{

}


void Optimizer::computeObsCost()
{

}


void Optimizer::computeSpeedGrad()
{

}


void Optimizer::computeObsGrad()
{

}

} // tpg
