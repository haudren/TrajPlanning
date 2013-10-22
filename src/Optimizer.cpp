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
    wpdata.jacVelMat.emplace_back(3, mb_.nrDof());
  }
  for(const Sphere& s: oc.collisionSpheres)
  {
    wpdata.jacSphere.emplace_back(mb_, s.bodyId, s.position);
    wpdata.jacSphereMat.emplace_back(3, mb_.nrDof());
  }
  wpData_.resize(oc.nrWp, wpdata);

  // start-end poses
  // start
  rbd::MultiBodyConfig mbcStartEnd(oc.start);
  rbd::forwardKinematics(mb_, mbcStartEnd);
  startPoses_ = mbcStartEnd.bodyPosW;
  startBodyVel_.resize(mb_.nrBodies());
  startBodyVelNorm_.resize(mb_.nrBodies());
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
  obsGrad_.resize(path_.rows());
  speedGrad_.resize(path_.rows());
}


Eigen::VectorXd Optimizer::path() const
{
  return path_;
}


std::vector<rbd::MultiBodyConfig> Optimizer::pathMbc() const
{
  std::vector<rbd::MultiBodyConfig> ret(wpData_.size());
  for(int wp = 0; wp < int(wpData_.size()); ++wp)
  {
    const WPData& data = wpData_[wp];
    ret[wp] = data.mbc;
    rbd::vectorToParam(path_.segment(wp*mb_.nrParams(), mb_.nrParams()),
                       ret[wp].q);
    rbd::forwardKinematics(mb_, ret[wp]);
  }
  return std::move(ret);
}


std::vector<IterResult> Optimizer::iters() const
{
  return iters_;
}


void Optimizer::optimize(int nrIter, double learningRate, const Eigen::VectorXd& initPath)
{
  iters_.clear();
  iters_.reserve(nrIter);

  path_ = initPath;

  for(int iter = 0; iter < nrIter; ++iter)
  {
    computeFK();
    computeVel();
    computeJac();

    computeObsCost();
    computeSpeedCost();
    computeSmCost();

    computeObsGrad();
    computeSpeedGrad();

    path_ += learningRate*AInv_.solve(obsGrad_ + speedGrad_);

    iters_.push_back({obsCost_, speedCost_, smCost_});
  }
}


void Optimizer::computeFK()
{
  for(int wp = 0; wp < int(wpData_.size()); ++wp)
  {
    WPData& data = wpData_[wp];
    rbd::vectorToParam(path_.segment(wp*mb_.nrParams(), mb_.nrParams()),
                       data.mbc.q);
    rbd::forwardKinematics(mb_, data.mbc);
  }
}


void Optimizer::computeVel()
{
  // start
  computeVel(startBodyVel_, startBodyVelNorm_, startPoses_, wpData_.front().mbc.bodyPosW);
  for(std::size_t wp = 0; wp < wpData_.size() - 1; ++wp)
  {
    WPData& data = wpData_[wp];
    WPData& dataNext = wpData_[wp + 1];
    computeVel(data.bodyVel, data.bodyVelNorm, data.mbc.bodyPosW, dataNext.mbc.bodyPosW);
  }
  // last wp
  computeVel(wpData_.back().bodyVel, wpData_.back().bodyVelNorm,
             wpData_.back().mbc.bodyPosW, endPoses_);
}


void Optimizer::computeVel(std::vector<Eigen::Vector3d>& bodyVel,
                           Eigen::ArrayXd& bodyVelNorm,
                           const std::vector<sva::PTransformd>& poses,
                           const std::vector<sva::PTransformd>& nextPoses)
{
  for(std::size_t bi = 0; bi < poses.size(); ++bi)
  {
    bodyVel[bi] = nextPoses[bi].translation() - poses[bi].translation();
    bodyVelNorm[bi] = std::sqrt(bodyVel[bi].squaredNorm() + 1);
  }
}


void Optimizer::computeJac()
{
  for(int wp = 0; wp < int(wpData_.size()); ++wp)
  {
    WPData& data = wpData_[wp];
    for(std::size_t ji = 0; ji < data.jacVel.size(); ++ji)
    {
      const Eigen::MatrixXd& jac = data.jacVel[ji].jacobian(mb_, data.mbc);
      data.jacVel[ji].fullJacobian(mb_, jac.block(3, 0, 3, jac.cols()), data.jacVelMat[ji]);
    }
    for(std::size_t js = 0; js < data.jacSphere.size(); ++js)
    {
      const Eigen::MatrixXd& jac = data.jacSphere[js].jacobian(mb_, data.mbc);
      data.jacSphere[js].fullJacobian(mb_, jac.block(3, 0, 3, jac.cols()), data.jacSphereMat[js]);
    }
  }
}


void Optimizer::computeSmCost()
{
  smCost_ = (path_.transpose()*A_*path_ + path_.transpose()*b_)[0] + c_;
}


void Optimizer::computeSpeedCost()
{
  speedCost_ = startBodyVelNorm_.sum();
  for(const WPData& wp: wpData_)
  {
    speedCost_ += wp.bodyVelNorm.sum();
  }
}


void Optimizer::computeObsCost()
{
  obsCost_ = 0.;
  for(const WPData& wp: wpData_)
  {
    for(const SphereData& sd: collisionSphere_)
    {
      obsCost_ += pen_.penality((sd.position*wp.mbc.bodyPosW[sd.bodyIndex]).translation()) - sd.radius;
    }
  }
}


void Optimizer::computeSpeedGrad()
{
  speedGrad_.setZero();
  // start
  {
    int startq0 = 0;
    const WPData& wp0 = wpData_.front();
    for(std::size_t j = 0; j < wp0.jacVelMat.size(); ++j)
    {
      Eigen::Vector3d velN = startBodyVel_[j]/startBodyVelNorm_[j];
      speedGrad_.segment(startq0, mb_.nrParams()) += wp0.jacVelMat[j].transpose()*velN;
    }
  }

  {
    int startqi = 0;
    int startqip1 =  mb_.nrParams();
    for(std::size_t i = 0; i < wpData_.size() - 1; ++i)
    {
      const WPData& wpi = wpData_[i];
      const WPData& wpip1 = wpData_[i + 1];
      for(std::size_t j = 0; j < wpi.jacVelMat.size(); ++j)
      {
        Eigen::Vector3d velN = wpi.bodyVel[j]/wpi.bodyVelNorm[j];
        speedGrad_.segment(startqi, mb_.nrParams()) -= wpi.jacVelMat[j].transpose()*velN;
        speedGrad_.segment(startqip1, mb_.nrParams()) += wpip1.jacVelMat[j].transpose()*velN;
      }
      startqi += mb_.nrParams();
      startqip1 += mb_.nrParams();
    }
  }

  // last wp
  {
    int startql = int(mb_.nrParams()*(wpData_.size() - 1));
    const WPData& wpl = wpData_.back();
    for(std::size_t j = 0; j < wpl.jacVelMat.size(); ++j)
    {
      Eigen::Vector3d velN = wpl.bodyVel[j]/wpl.bodyVelNorm[j];
      speedGrad_.segment(startql, mb_.nrParams()) -= wpl.jacVelMat[j].transpose()*velN;
    }
  }
}


void Optimizer::computeObsGrad()
{
  obsGrad_.setZero();
  int start = 0;
  for(const WPData& wp: wpData_)
  {
    for(std::size_t i = 0; i < wp.jacSphere.size(); ++i)
    {
      const SphereData& sd = collisionSphere_[i];
      Eigen::Vector3d position = (sd.position*wp.mbc.bodyPosW[sd.bodyIndex]).translation();
      obsGrad_.segment(start, mb_.nrParams()) +=
          wp.jacSphereMat[i].transpose()*pen_.penalityGrad(position);
    }
    start += mb_.nrParams();
  }
}

} // tpg
