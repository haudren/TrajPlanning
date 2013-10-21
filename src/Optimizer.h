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
// stl
#include <vector>

// Eigen
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

// RBDyn
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>
#include <RBDyn/Jacobian.h>

// TrajPlanning
#include "ObsPen.h"

namespace tpg
{

struct Sphere
{
  int bodyId;
  double radius;
  Eigen::Vector3d position;
};


struct OptimizerConfig
{
  int nrWp;
  ObsPen pen;
  rbd::MultiBody mb;
  rbd::MultiBodyConfig start;
  rbd::MultiBodyConfig end;
  std::vector<Sphere> collisionSpheres;
  double velWeight;
  double accWeight;
  double jerkWeight;
};


class Optimizer
{
public:
  void init(const OptimizerConfig& oc);
  void optimize(int nrIter, double learningRate, const Eigen::VectorXd& initPath);

private:
  struct WPData
  {
    rbd::MultiBodyConfig mbc;
    std::vector<Eigen::Vector3d> bodyVel;
    std::vector<double> bodyVelNorm;
    std::vector<rbd::Jacobian> jacVel;
    std::vector<rbd::Jacobian> jacSphere;
  };

  struct SphereData
  {
    int bodyIndex;
    double radius;
    sva::PTransformd position;
  };

private:
  void computeFK();
  void computeVel();
  void computeVel(WPData& data, const std::vector<sva::PTransformd>& nextPoses);
  void computeJac();

  void computeSmCost();
  void computeSpeedCost();
  void computeObsCost();

  void computeSpeedGrad();
  void computeObsGrad();

private:
  rbd::MultiBody mb_;
  rbd::MultiBodyConfig mbc_;
  ObsPen pen_;
  std::vector<SphereData> collisionSphere_;
  std::vector<WPData> wpData_;
  std::vector<sva::PTransformd> startPoses_, endPoses_;

  // optim
  Eigen::VectorXd path_;
  double obsCost_, speedCost_, smCost_;
  Eigen::VectorXd obsGrad_, speedGrad_;

  // smoothness matrix
  Eigen::SparseMatrix<double> A_;
  Eigen::SparseVector<double> b_;
  double c_;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > AInv_;
};

} // tpg
