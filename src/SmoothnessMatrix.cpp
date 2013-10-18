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
#include "SmoothnessMatrix.h"

#include <iostream>
// Eigen
#include <unsupported/Eigen/KroneckerProduct>

namespace tpg
{

Eigen::SparseMatrix<double> smoothnessRec(int nrWp, int der)
{
  typedef Eigen::Triplet<double> Trip;
  std::vector<Trip> tripletList;
  tripletList.reserve(2*(nrWp - 1));
  for(int i = 0; i < (nrWp - 1); ++i)
  {
    tripletList.push_back(Trip(i, i, -1.));
    tripletList.push_back(Trip(i, i + 1,  1.));
  }

  Eigen::SparseMatrix<double> sm(nrWp - 1, nrWp);
  sm.setFromTriplets(tripletList.begin(), tripletList.end());

  if((der - 1) > 0)
  {
    return std::move(smoothnessRec(nrWp - 1, der - 1)*sm);
  }
  else
  {
    return std::move(sm);
  }
}


Eigen::SparseMatrix<double> smoothness(int nrWp, int der)
{
  return std::move(smoothnessRec(nrWp + 2, der).block(0, 1, nrWp + 2 - der, nrWp));
}


Eigen::SparseMatrix<double> smoothness(int nrParams, int nrWp, int der)
{
  Eigen::SparseMatrix<double> I(nrParams, nrParams);
  I.setIdentity();
  return std::move(Eigen::kroneckerProduct(smoothness(nrWp, der), I));
}


Eigen::SparseVector<double> smoothnessVector(int nrWp, int der,
    const Eigen::VectorXd& iQ, const Eigen::VectorXd& fQ, double iC, double fC)
{
  int nrParams = int(iQ.rows());
  typedef Eigen::Triplet<double> Trip;
  std::vector<Trip> tripletList;
  tripletList.reserve(nrParams*2);

  for(int i = 0; i < nrParams; ++i)
  {
    tripletList.push_back(Trip(i, 0, iC*iQ(i)));
  }
  int start = nrParams*(nrWp - der + 1);
  for(int i = 0; i < nrParams; ++i)
  {
    tripletList.push_back(Trip(start + i, 0, fC*fQ(i)));
  }

  Eigen::SparseMatrix<double> d1(nrParams*(nrWp - der + 2), 1);
  d1.setFromTriplets(tripletList.begin(), tripletList.end());

  return std::move(Eigen::SparseMatrix<double>(d1));
}


QuadTuple velocitySmoothness(int nrWp,
                             const Eigen::VectorXd& iQ, const Eigen::VectorXd& fQ,
                             double velW)
{
  int nrParams = int(iQ.rows());
  auto K1 = smoothness(nrParams, nrWp, 1);
  auto d1 = smoothnessVector(nrWp, 1, iQ, fQ, -1., 1.);

  auto A = velW*K1.transpose()*K1;
  auto b = 2.*velW*K1.transpose()*d1;
  auto c = velW*d1.dot(d1);

  return std::make_tuple(A, b, c);
}


QuadTuple accelerationSmoothness(int nrWp,
                                 const Eigen::VectorXd& iQ, const Eigen::VectorXd& fQ,
                                 double velW, double accW)
{
  int nrParams = int(iQ.rows());
  auto K2 = smoothness(nrParams, nrWp, 2);
  auto d2 = smoothnessVector(nrWp, 2, iQ, fQ, 1., 1.);
  auto velSmooth = velocitySmoothness(nrWp, iQ, fQ, velW);

  auto A = accW*K2.transpose()*K2 + std::get<0>(velSmooth);
  auto b = accW*2.*K2.transpose()*d2 + std::get<1>(velSmooth);
  auto c = accW*d2.dot(d2) + std::get<2>(velSmooth);

  return std::make_tuple(A, b, c);
}


QuadTuple jerkSmoothness(int nrWp,
                         const Eigen::VectorXd& iQ, const Eigen::VectorXd& fQ,
                         double velW, double accW, double jerkW)
{
  int nrParams = int(iQ.rows());
  auto K3 = smoothness(nrParams, nrWp, 3);
  auto d3 = smoothnessVector(nrWp, 3, iQ, fQ, -1., 1.);
  auto accSmooth = accelerationSmoothness(nrWp, iQ, fQ, velW, accW);

  auto A = jerkW*K3.transpose()*K3 + std::get<0>(accSmooth);
  auto b = jerkW*2.*K3.transpose()*d3 + std::get<1>(accSmooth);
  auto c = jerkW*d3.dot(d3) + std::get<2>(accSmooth);

  return std::make_tuple(A, b, c);
}

} // tpg
