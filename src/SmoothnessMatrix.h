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
// std
#include <tuple>

// Eigen
#include <Eigen/Sparse>

namespace tpg
{

Eigen::SparseMatrix<double> smoothnessRec(int nrWp, int der);
Eigen::SparseMatrix<double> smoothness(int nrWp, int der);
Eigen::SparseMatrix<double> smoothness(int nrParams, int nrWp, int der);

typedef std::tuple<Eigen::SparseMatrix<double>,
                   Eigen::SparseVector<double>,
                   double> QuadTuple;

QuadTuple velocitySmoothness(int nrWp,
                             const Eigen::VectorXd& iQ, const Eigen::VectorXd& fQ,
                             double velW);
QuadTuple accelerationSmoothness(int nrWp,
                                 const Eigen::VectorXd& iQ, const Eigen::VectorXd& fQ,
                                 double velW, double accW);
QuadTuple jerkSmoothness(int nrWp,
                             const Eigen::VectorXd& iQ, const Eigen::VectorXd& fQ,
                             double velW, double accW, double jerkW);

} // tpg
