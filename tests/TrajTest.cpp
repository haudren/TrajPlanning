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

// include
// std
#include <fstream>
#include <iostream>
#include <tuple>

// boost
#define BOOST_TEST_MODULE Algo test
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>

// RBDyn
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>
#include <RBDyn/MultiBodyGraph.h>

// TrajPlanning
#include "ObsPen.h"
#include "Optimizer.h"
#include "SmoothnessMatrix.h"

const Eigen::Vector3d gravity(0., 9.81, 0.);


/// @return An simple Z*6 arm with Y as up axis.
std::tuple<rbd::MultiBody, rbd::MultiBodyConfig> makeZ6Arm(bool isFixed=true)
{
  using namespace Eigen;
  using namespace sva;
  using namespace rbd;

  MultiBodyGraph mbg;

  double mass = 1.;
  Matrix3d I = Matrix3d::Identity();
  Vector3d h = Vector3d::Zero();

  RBInertiad rbi(mass, h, I);

  for(int i = 0; i < 7; ++i)
  {
    std::stringstream ss;
    ss << "b" << i;
    mbg.addBody({rbi, i, ss.str()});
  }

  for(int i = 0; i < 6; ++i)
  {
    std::stringstream ss;
    ss << "j" << i;
    mbg.addJoint({Joint::RevZ, true, i, ss.str()});
  }

  PTransformd to(Vector3d(0., 1., 0.));
  PTransformd from(Vector3d(0., 0., 0.));

  mbg.linkBodies(0, from, 1, from, 0);
  for(int i = 1; i < 6; ++i)
  {
    mbg.linkBodies(i, to, i + 1, from, i);
  }

  MultiBody mb = mbg.makeMultiBody(0, isFixed);

  MultiBodyConfig mbc(mb);
  mbc.zero(mb);

  return std::make_tuple(mb, mbc);
}


BOOST_AUTO_TEST_CASE(ObsPenTest)
{
  using namespace Eigen;

  std::vector<double> pen =
  {
    0.19338496,  0.08683781,  0.76232272,  0.10032556,  0.0311701 ,
    0.74030221,  0.4986186 ,  0.58865215,  0.63947176,  0.37107554,
    0.77703448,  0.94472095,  0.16495522,  0.2537881 ,  0.12636114,
    0.89338157,  0.18361576,  0.7980018 ,  0.58179607,  0.19751129,
    0.2025195 ,  0.70539315,  0.06764872,  0.90060331,  0.23950046,
    0.30006224,  0.33118872
  };
  std::vector<double> penGradX =
  {
    0.17769058,  0.69019667,  0.18239823,  0.06462966,  0.222618  ,
   -0.61394108,  0.39476297, -0.40503639,  0.15853003,  0.19420555,
    0.05533674, -0.27990161,  0.30253379,  0.01823931,  0.08015055,
   -0.12955907, -0.14429495, -0.15414152,  0.21072053, -0.57952318,
   -0.74220145,  0.54043792, -0.18613938,  0.77424217, -0.65388111,
    0.11644648, -0.46681308
  };
  std::vector<double> penGradY =
  {
    -0.0930594 , -0.0556677 , -0.02202051,  0.15261682,  0.25090717,
    -0.06142548,  0.39829304,  0.55748205, -0.10083045, -0.20612032,
    -0.52324638, -0.81835981,  0.26115301, -0.29670936, -0.07335958,
     0.72842634, -0.07017234,  0.67164066,  0.12359707, -0.12986258,
     0.6980838 , -0.17114781,  0.05127547,  0.06433461, -0.46589269,
     0.23241353, -0.56941459
  };
  std::vector<double> penGradZ =
  {
    -0.10654716,  0.28446888,  0.67548491, -0.06915546,  0.31998833,
     0.70913211,  0.09003355,  0.07042658,  0.05081961,  0.40595893,
     0.2868227 ,  0.16768647,  0.08883288, -0.01929704, -0.12742697,
    -0.70976581, -0.04768988,  0.61438604, -0.38428478, -0.18963829,
     0.00500821, -0.63774443,  0.09760508,  0.83295459,  0.06056178,
     0.04584413,  0.03112647
  };

  std::vector<Vector3d> points =
  {
    {0.23352768013842229, 1.7767243859557733, 0.26601814042029426},
    {2.2485555785893334, 0.55614005239213393, 1.269660653465275},
    {2.1185875346180856, 2.071751514751571, 1.509121136794852},
    {1.8537692322878447, 1.2411053144925532, 1.4575584431698618},
    {1.3542112341981949, 0.85807431414450286, 2.6060901323821746}
  };

  std::vector<double> penRes =
  {
    0.46275887,  0.        ,  0.        ,  0.39354775,  0.
  };

  std::vector<Vector3d> penGradRes =
  {
    {0.11422513837340077, 0.39239411239478911, -0.037287388077032885},
    {0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0},
    {0.13306870948137786, -0.0026793299369999649, 0.29066137572319539},
    {0.0, 0.0, 0.0}
  };


  tpg::ObsPen obsp;
  obsp.setPen(Vector3d::Zero(), Vector3d(1., 1., 1.), 3, 3, 3,
              pen, penGradX, penGradY, penGradZ);

  for(std::size_t i = 0; i < points.size(); ++i)
  {
    BOOST_CHECK_SMALL(obsp.penality(points[i]) - penRes[i], 1e-4);
    BOOST_CHECK_SMALL((obsp.penalityGrad(points[i]) - penGradRes[i]).norm(), 1e-4);
  }
}


BOOST_AUTO_TEST_CASE(VelSmoothnessTest)
{
  using namespace Eigen;

  MatrixXd resA(9,9);
  VectorXd resB(9);
  double resC;

  resA << 10.0, 0.0, 0.0, -5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0,
          -5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, -5.0, 0.0, 0.0,
          0.0, -5.0, 0.0, 0.0, 10.0, 0.0, 0.0, -5.0, 0.0, 0.0, 0.0, -5.0, 0.0,
          0.0, 10.0, 0.0, 0.0, -5.0, 0.0, 0.0, 0.0, -5.0, 0.0, 0.0, 10.0, 0.0,
          0.0, -5.0, 0.0, 0.0, 0.0, -5.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, -5.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -5.0, 0.0, 0.0, 10.0;
  resB << 0.,   0.,   0.,   0.,   0.,   0., -10., -10., -10.;
  resC = 15.;

  VectorXd initQ(3);
  initQ << 0., 0., 0.;
  VectorXd finalQ(3);
  finalQ << 1., 1., 1.;
  auto abc = tpg::velocitySmoothness(3, initQ, finalQ, 5.);

  BOOST_CHECK_SMALL((resA - MatrixXd(std::get<0>(abc))).norm(), 1e-4);
  BOOST_CHECK_SMALL((resB - VectorXd(std::get<1>(abc))).norm(), 1e-4);
  BOOST_CHECK_SMALL(std::abs(resC - std::get<2>(abc)), 1e-4);
}


BOOST_AUTO_TEST_CASE(AccSmoothnessTest)
{
  using namespace Eigen;

  MatrixXd resA(12,12);
  VectorXd resB(12);
  double resC;

  resA << 22.5, 0.0, 0.0, 0.0, -16.5, 0.0, 0.0, 0.0, 3.5, 0.0, 0.0, 0.0, 0.0,
      22.5, 0.0, 0.0, 0.0, -16.5, 0.0, 0.0, 0.0, 3.5, 0.0, 0.0, 0.0, 0.0, 22.5,
      0.0, 0.0, 0.0, -16.5, 0.0, 0.0, 0.0, 3.5, 0.0, 0.0, 0.0, 0.0, 22.5, 0.0,
      0.0, 0.0, -16.5, 0.0, 0.0, 0.0, 3.5, -16.5, 0.0, 0.0, 0.0, 26.0, 0.0,
      0.0, 0.0, -16.5, 0.0, 0.0, 0.0, 0.0, -16.5, 0.0, 0.0, 0.0, 26.0, 0.0, 0.0,
      0.0, -16.5, 0.0, 0.0, 0.0, 0.0, -16.5, 0.0, 0.0, 0.0, 26.0, 0.0, 0.0, 0.0,
      -16.5, 0.0, 0.0, 0.0, 0.0, -16.5, 0.0, 0.0, 0.0, 26.0, 0.0, 0.0, 0.0, -16.5,
      3.5, 0.0, 0.0, 0.0, -16.5, 0.0, 0.0, 0.0, 22.5, 0.0, 0.0, 0.0, 0.0, 3.5, 0.0,
      0.0, 0.0, -16.5, 0.0, 0.0, 0.0, 22.5, 0.0, 0.0, 0.0, 0.0, 3.5, 0.0, 0.0, 0.0,
      -16.5, 0.0, 0.0, 0.0, 22.5, 0.0, 0.0, 0.0, 0.0, 3.5, 0.0, 0.0, 0.0, -16.5, 0.0,
      0.0, 0.0, 22.5;
  resB <<-57.0, -38.0, -19.0, 0.0, 35.0, 28.0, 14.0, 14.0, -38.0, -38.0, -19.0, -38.0;
  resC = 162.;

  VectorXd initQ(4);
  initQ << 3.,2.,1.,0.;
  VectorXd finalQ(4);
  finalQ << 2.,2.,1.,2.;
  auto abc = tpg::accelerationSmoothness(3, initQ, finalQ, 2.5, 3.5);

  BOOST_CHECK_SMALL((resA - MatrixXd(std::get<0>(abc))).norm(), 1e-4);
  BOOST_CHECK_SMALL((resB - VectorXd(std::get<1>(abc))).norm(), 1e-4);
  BOOST_CHECK_SMALL(std::abs(resC - std::get<2>(abc)), 1e-4);
}


BOOST_AUTO_TEST_CASE(JerkSmoothnessTest)
{
  using namespace Eigen;

  MatrixXd resA(12,12);
  VectorXd resB(12);
  double resC;

  resA << 62.5, 0.0, 0.0, -64.5, 0.0, 0.0, 27.5, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0,
      62.5, 0.0, 0.0, -64.5, 0.0, 0.0, 27.5, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 62.5,
      0.0, 0.0, -64.5, 0.0, 0.0, 27.5, 0.0, 0.0, -4.0, -64.5, 0.0, 0.0, 102.0,
      0.0, 0.0, -76.5, 0.0, 0.0, 27.5, 0.0, 0.0, 0.0, -64.5, 0.0, 0.0, 102.0, 0.0,
      0.0, -76.5, 0.0, 0.0, 27.5, 0.0, 0.0, 0.0, -64.5, 0.0, 0.0, 102.0, 0.0, 0.0,
      -76.5, 0.0, 0.0, 27.5, 27.5, 0.0, 0.0, -76.5, 0.0, 0.0, 102.0, 0.0, 0.0, -64.5,
      0.0, 0.0, 0.0, 27.5, 0.0, 0.0, -76.5, 0.0, 0.0, 102.0, 0.0, 0.0, -64.5, 0.0, 0.0,
      0.0, 27.5, 0.0, 0.0, -76.5, 0.0, 0.0, 102.0, 0.0, 0.0, -64.5, -4.0, 0.0, 0.0,
      27.5, 0.0, 0.0, -64.5, 0.0, 0.0, 62.5, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 27.5, 0.0,
      0.0, -64.5, 0.0, 0.0, 62.5, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 27.5, 0.0, 0.0,
      -64.5, 0.0, 0.0, 62.5;
  resB << -129.0, -86.0, -43.0, 53.0, 38.0, 23.0, 131.0, 77.0, 23.0, -215.0, -129.0, -43.0;
  resC = 490;

  VectorXd initQ(3);
  initQ << 3.,2.,1.;
  VectorXd finalQ(3);
  finalQ << 5.,3.,1.;
  auto abc = tpg::jerkSmoothness(4, initQ, finalQ, 2.5, 3.5, 4.);

  BOOST_CHECK_SMALL((resA - MatrixXd(std::get<0>(abc))).norm(), 1e-4);
  BOOST_CHECK_SMALL((resB - VectorXd(std::get<1>(abc))).norm(), 1e-4);
  BOOST_CHECK_SMALL(std::abs(resC - std::get<2>(abc)), 1e-4);
}


struct PenMap
{
  int sizeX, sizeY, sizeZ;
  Eigen::Vector3d scale;
  std::vector<double> penality;
  std::vector<double> penalityGradX;
  std::vector<double> penalityGradY;
  std::vector<double> penalityGradZ;
};


PenMap loadPenalityMap(const std::string& filename)
{
  std::ifstream mapfile(filename);
  BOOST_REQUIRE(mapfile.is_open());

  PenMap map;
  mapfile >> map.sizeX >> map.sizeY >> map.sizeZ;
  mapfile >> map.scale.x() >> map.scale.y() >> map.scale.z();
  map.penality.resize(map.sizeX*map.sizeY*map.sizeZ);
  map.penalityGradX.resize(map.penality.size());
  map.penalityGradY.resize(map.penality.size());
  map.penalityGradZ.resize(map.penality.size());

  for(std::size_t i = 0; i < map.penality.size(); ++i)
  {
    mapfile >> map.penality[i];
  }
  for(std::size_t i = 0; i < map.penality.size(); ++i)
  {
    mapfile >> map.penalityGradX[i];
  }
  for(std::size_t i = 0; i < map.penality.size(); ++i)
  {
    mapfile >> map.penalityGradY[i];
  }
  for(std::size_t i = 0; i < map.penality.size(); ++i)
  {
    mapfile >> map.penalityGradZ[i];
  }

  return std::move(map);
}


void writeFK(std::ofstream& outfile, const rbd::MultiBodyConfig& mbc)
{
  outfile << "[";
  for(const sva::PTransformd& pt: mbc.bodyPosW)
  {
    outfile << "[";
    outfile << pt.translation()[0] << "," << pt.translation()[1] << "," <<
               pt.translation()[2];
    outfile << "],";
  }
  outfile << "],";
  outfile << std::endl;
}


void writeResult(const std::string& filename, const tpg::Optimizer& opt,
                 const rbd::MultiBodyConfig& start, const rbd::MultiBodyConfig& end)
{
  std::ofstream outfile(filename);
  BOOST_REQUIRE(outfile.is_open());

  outfile << "obsCost=[";
  for(const tpg::IterResult& it: opt.iters())
  {
    outfile << it.obsCost << ",";
  }
  outfile << "]" << std::endl;
  outfile << "smCost=[";
  for(const tpg::IterResult& it: opt.iters())
  {
    outfile << it.smCost << ",";
  }
  outfile << "]" << std::endl;
  outfile << "speedCost=[";
  for(const tpg::IterResult& it: opt.iters())
  {
    outfile << it.speedCost << ",";
  }
  outfile << "]" << std::endl;

  outfile << "bodyPos=[";
  writeFK(outfile, start);
  for(const rbd::MultiBodyConfig& mbc: opt.pathMbc())
  {
    writeFK(outfile, mbc);
  }
  writeFK(outfile, end);
  outfile << "]" << std::endl;

  outfile << "penMap=[";
  const tpg::ObsPen& op = opt.obsPen();
  for(int y = 0; y < 10; ++y)
  {
    outfile << "[";
    for(int x = 0; x < 10; ++x)
    {
      outfile << op.penality(Eigen::Vector3d(x, y, 0.)) << ",";
    }
    outfile << "]," << std::endl;
  }
  outfile << "]" << std::endl;
}


Eigen::VectorXd basicInterp(const rbd::MultiBody& mb,
                            const rbd::MultiBodyConfig& mbcStart,
                            const rbd::MultiBodyConfig& mbcEnd,
                            int nrWp)
{
  Eigen::VectorXd ret(mb.nrParams()*nrWp);
  Eigen::VectorXd start = rbd::paramToVector(mb, mbcStart.q);
  Eigen::VectorXd end = rbd::paramToVector(mb, mbcEnd.q);

  for(int i = 0; i < nrWp; ++i)
  {
    ret.segment(i*mb.nrParams(), mb.nrParams()) =
        double(i + 1)/(nrWp + 1)*(end - start);
  }

  return ret;
}


void testBasicTrajectory(const std::string& map, const std::string& out)
{
  using namespace Eigen;
  using namespace sva;
  using namespace rbd;
  namespace cst = boost::math::constants;

  MultiBody mb;
  MultiBodyConfig mbcInit, mbcWork;

  std::tie(mb, mbcInit) = makeZ6Arm();
  mbcWork = mbcInit;
  PenMap penMap = loadPenalityMap(map);

  tpg::ObsPen obsPen;
  obsPen.setPen(Vector3d(0., 0., -5.), penMap.scale,
                penMap.sizeX, penMap.sizeY, penMap.sizeZ,
                penMap.penality, penMap.penalityGradX, penMap.penalityGradY,
                penMap.penalityGradZ);

  tpg::Optimizer opt;
  tpg::OptimizerConfig conf;
  conf.nrWp = 30;
  conf.pen = obsPen;
  conf.mb = mb;
  conf.start = mbcWork;
  rbd::forwardKinematics(mb, conf.start);
  mbcWork.q[1][0] = -cst::pi<double>()/2.;
  conf.end = mbcWork;
  rbd::forwardKinematics(mb, conf.end);
  for(int i = 0; i < 6; ++i)
  {
    conf.collisionSpheres.push_back({i + 1, 0., Vector3d(0., 0.5, 0.)});
  }
  conf.velWeight = 1.;
  conf.accWeight = 0.;
  conf.jerkWeight = 0.;

  opt.init(conf);

  opt.optimize(100, 0.01, basicInterp(mb, conf.start, conf.end, conf.nrWp));
  writeResult(out, opt, conf.start, conf.end);
}


BOOST_AUTO_TEST_CASE(FreeTraj)
{
  testBasicTrajectory("empty1.map", "freetraj.py");
}


BOOST_AUTO_TEST_CASE(ObsTraj)
{
  testBasicTrajectory("obs1.map", "obs1traj.py");
  testBasicTrajectory("obs2.map", "obs2traj.py");
  testBasicTrajectory("obs3.map", "obs3traj.py");
}

