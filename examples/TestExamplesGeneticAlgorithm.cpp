// Unfit: Data fitting and optimization software
//
// Copyright (C) 2012- Dr Martin Buist & Dr Alberto Corrias
// Contacts: martin.buist _at_ nus.edu.sg; alberto _at_ nus.edu.sg
//
// See the 'Contributors' file for a list of those who have contributed
// to this work.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
#include <chrono>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "CardiacAlphaN.hpp"
#include "DataFileReader.hpp"
#include "Exponential.hpp"
#include "GaussianEquation.hpp"
#include "GeneticAlgorithm.hpp"
#include "GeneticSwitch.hpp"
#include "HelicalValleyCostFunction.hpp"
#include "HodgkinHuxleyBetaN.hpp"
#include "LevenbergMarquardt.hpp"
#include "Matrix.hpp"
#include "Modroslam.hpp"
#include "NelderMead.hpp"
#include "NonStationaryMarkov.hpp"
#include "ODE2D.hpp"
#include "ODE3DVariant.hpp"
#include "Osborne.hpp"
#include "Parabolic.hpp"
#include "Powell.hpp"
#include "Rosenbrock.hpp"
#include "ThreeNaMarkov.hpp"
#include "UnitTest++.h"
#include "Woods.hpp"

namespace  // file scope
{
typedef std::chrono::high_resolution_clock hrclock_t;

// Needed because at present mingw has not implemented this (C++11).
// Replace to_string with std::to_string when possible
template <typename T> std::string to_string(T duration)
{
  std::stringstream converter;
  converter << duration;
  return converter.str();
}

std::string TestTime(hrclock_t::time_point t1, hrclock_t::time_point t2)
{
  using namespace std::chrono;
  auto elapsed = duration_cast<nanoseconds>(t2 - t1).count();
  if (elapsed < 10000) return {to_string(elapsed) + "\t nsec for : "};
  elapsed = duration_cast<microseconds>(t2 - t1).count();
  if (elapsed < 10000) return {to_string(elapsed) + "\t usec for : "};
  elapsed = duration_cast<milliseconds>(t2 - t1).count();
  if (elapsed < 10000) return {to_string(elapsed) + "\t msec for : "};
  elapsed = duration_cast<seconds>(t2 - t1).count();
  return {to_string(elapsed) + "\t sec  for : "};
}
}  // namespace

SUITE(ExamplesGeneticAlgorithm)
{
TEST(GeneticAlgorithm_CardiacAlphaN)  // From Unfit 1
{
  Unfit::GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, -100.0, 100.0);
  ga.bounds.SetBounds(1, -100.0, 100.0);
  ga.bounds.SetBounds(2, -10.0, 10.0);
  ga.options.SetGamma(0.4);
  ga.options.SetCostTolerance(1e-3);
  // Read in the experimental data
  std::vector<double> vm;
  std::vector<double> alpha_n;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/alphadata.csv"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, vm));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, alpha_n));
  CHECK_EQUAL(vm.size(), alpha_n.size());
  // Create the cost function
  Unfit::Examples::CardiacAlphaN cost_func(vm, alpha_n);

  // Initial guess
  std::vector<double> min_point {0.08, -0.04, 3.0};
  // Check the function calculates the correct cost at the initial guess
  auto residual = cost_func(min_point);
  auto sum_sq_residual = Unfit::SumOfSquares(residual);
  CHECK_CLOSE(0.000853649, sum_sq_residual, 1e-8);

  // Minimise
  auto t1 = hrclock_t::now();  // Start time
  auto rc = ga.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "CardiacAlphaN (GeneticAlgorithm)";
  std::cout << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.028540881, min_point[0], 1e-4);
  CHECK_CLOSE(-0.19883019, min_point[1], 1e-4);
  CHECK_CLOSE(9.7910215, min_point[2], 1e-4);

  // Check we can find the minimum from the GA initial guess with NM
  auto nm_min_point = min_point;
  Unfit::NelderMead nm;
  rc = nm.FindMin(cost_func, nm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0731685, nm_min_point[0], 1e-4);
  CHECK_CLOSE(-0.033983, nm_min_point[1], 1e-4);
  CHECK_CLOSE(3.10956, nm_min_point[2], 1e-4);

  // Check we can find the minimum from the GA initial guess with LM
  auto lm_min_point = min_point;
  Unfit::LevenbergMarquardt lm;
  rc = lm.FindMin(cost_func, lm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0731685, lm_min_point[0], 1e-4);
  CHECK_CLOSE(-0.033983, lm_min_point[1], 1e-4);
  CHECK_CLOSE(3.10956, lm_min_point[2], 1e-4);
}

TEST(GeneticAlgorithm_Exponential)
{
  Unfit::GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, 0.0, 100.0);
  ga.bounds.SetBounds(1, -10.0, 10.0);
  ga.bounds.SetBounds(2, -100.0, 100.0);
  ga.options.SetGamma(0.25);
  // Read in the experimental data
  std::vector<double> x;
  std::vector<double> y;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/exponential_exp_data.dat"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, x));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, y));
  CHECK_EQUAL(x.size(), y.size());
  // Create the cost function
  Unfit::Examples::Exponential cost_func(x, y);

  // Initial guess
  std::vector<double> min_point {1.0, 1.0, 1.0};
  // Check the function calculates the correct cost at the initial guess
  auto residual = cost_func(min_point);
  auto sum_sq_residual = Unfit::SumOfSquares(residual);
  CHECK_CLOSE(1793.39, sum_sq_residual, 1e-2);

  // Minimise
  auto t1 = hrclock_t::now();  // Start time
  auto rc = ga.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Exponential (GeneticAlgorithm)";
  std::cout << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(2.5060727, min_point[0], 0.03);
  CHECK_CLOSE(0.071810527, min_point[1], 0.03);
  CHECK_CLOSE(2.9201527, min_point[2], 0.2);

  // Check we can find the minimum from the GA initial guess with NM
  auto nm_min_point = min_point;
  Unfit::NelderMead nm;
  rc = nm.FindMin(cost_func, nm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(3.66696, nm_min_point[0], 1e-2);
  CHECK_CLOSE(0.0885905, nm_min_point[1], 1e-2);
  CHECK_CLOSE(3.06235, nm_min_point[2], 1e-2);

  // Check we can find the minimum from the GA initial guess with LM
  auto lm_min_point = min_point;
  Unfit::LevenbergMarquardt lm;
  rc = lm.FindMin(cost_func, lm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(3.66696, lm_min_point[0], 1e-2);
  CHECK_CLOSE(0.0885905, lm_min_point[1], 1e-2);
  CHECK_CLOSE(3.06235, lm_min_point[2], 1e-2);
}

TEST(GeneticAlgorithm_Gaussian)
{
  Unfit::GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, 0.0, 100.0);
  ga.bounds.SetBounds(1, -100.0, 100.0);
  ga.bounds.SetBounds(2, 0.0, 100.0);
  ga.bounds.SetBounds(3, -100.0, 100.0);
  ga.options.SetGeometricTolerance(1e-8);
  // Read in the experimental data
  std::vector<double> x;
  std::vector<double> y;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/gaussian_exp_data.dat"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, x));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, y));
  CHECK_EQUAL(x.size(), y.size());
  // Create the cost function
  Unfit::Examples::GaussianEquation cost_func(x, y);

  // Initial guess
  std::vector<double> min_point {1.0, 22.0, 3.0, 12.0};
  // Check the function calculates the correct cost at the initial guess
  auto residual = cost_func(min_point);
  auto sum_sq_residual = Unfit::SumOfSquares(residual);
  CHECK_CLOSE(2636.63, sum_sq_residual, 1e-2);

  // Minimise
  auto t1 = hrclock_t::now();  // Start time
  auto rc = ga.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Gaussian (GeneticAlgorithm)";
  std::cout << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(12.397833, min_point[0], 0.2);
  CHECK_CLOSE(22.745196, min_point[1], 0.2);
  CHECK_CLOSE(3.9811465, min_point[2], 0.2);
  CHECK_CLOSE(12.628567, min_point[3], 0.25);

  // Check we can find the minimum from the GA initial guess with NM
  auto nm_min_point = min_point;
  Unfit::NelderMead nm;
  rc = nm.FindMin(cost_func, nm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(12.4144, nm_min_point[0], 1e-4);
  CHECK_CLOSE(22.7355, nm_min_point[1], 1e-4);
  CHECK_CLOSE(4.0337, nm_min_point[2], 1e-4);
  CHECK_CLOSE(12.5486, nm_min_point[3], 1e-4);

  // Check we can find the minimum from the GA initial guess with LM
  auto lm_min_point = min_point;
  Unfit::LevenbergMarquardt lm;
  rc = lm.FindMin(cost_func, lm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(12.4144, lm_min_point[0], 1e-4);
  CHECK_CLOSE(22.7355, lm_min_point[1], 1e-4);
  CHECK_CLOSE(4.0337, lm_min_point[2], 1e-4);
  CHECK_CLOSE(12.5486, lm_min_point[3], 1e-4);
}

// MB: This test example does run, but it takes a long time and so has
//     been omitted from the regular tests for now (3x longer than all of
//     the other tests combined).
//
// TEST(GeneticAlgorithm_GeneticSwitch)  // From Unfit 1
// {
//  Unfit::GeneticAlgorithm ga;
//  ga.bounds.SetBounds(0, 0.0, 1000.0);
//  ga.bounds.SetBounds(1, 0.0, 10.0);
//  ga.bounds.SetBounds(2, -10.0, 10.0);
//  ga.bounds.SetBounds(3, -1.0, 1.0);
//  ga.options.SetCostTolerance(2500.0);
//  ga.options.SetPopulationSize(20u);
//  ga.options.SetGamma(0.25);
//  // Read in the experimental data
//  std::vector<double> t;
//  std::vector<double> z_a;
//  Unfit::DataFileReader<double> dfr;
//  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/training_data_zA.txt"));
//  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, t));
//  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, z_a));
//  auto dt = t[1] - t[0];
//  // Create the cost function
//  Unfit::Examples::GeneticSwitch cost_func(z_a, dt);
//
//  // Initial guess
//  std::vector<double> min_point = {90.0, 60.0, -1.0, -0.1};
//  // Check the function calculates the correct cost at the initial guess
//  auto residual = cost_func(min_point);
//  auto sum_sq_residual = Unfit::SumOfSquares(residual);
//  CHECK_CLOSE(1.12863e7, sum_sq_residual, 1e2);
//
//  // Minimise
//  auto t1 = hrclock_t::now();  // Start time
//  auto rc = ga.FindMin(cost_func, min_point);
//  auto t2 = hrclock_t::now();  // End time
//  std::cout << TestTime(t1, t2) << "GeneticSwitch (GeneticAlgorithm)";
//  std::cout << std::endl;
//
//  // Check the result matches what we expect
//  CHECK_EQUAL(0, rc);
//  CHECK_CLOSE(253.59291, min_point[0], 1e-3);
//  CHECK_CLOSE(2.3086556, min_point[1], 1e-3);
//  CHECK_CLOSE(-1.8484436, min_point[2], 1e-3);
//  CHECK_CLOSE(-0.52417421, min_point[3], 1e-3);
//
//  // Check we can find the minimum from the GA initial guess with NM
//  auto nm_min_point = min_point;
//  Unfit::NelderMead nm;
//  rc = nm.FindMin(cost_func, nm_min_point);
//  CHECK_EQUAL(0, rc);
//  CHECK_CLOSE(216.433109901728, nm_min_point[0], 1e-3);
//  CHECK_CLOSE(2.00476590272927, nm_min_point[1], 1e-3);
//  CHECK_CLOSE(-2.17034564257086, nm_min_point[2], 1e-3);
//  CHECK_CLOSE(-0.226712906915202, nm_min_point[3], 1e-3);
//
//  // Check we can find the minimum from the GA initial guess with LM
//  auto lm_min_point = min_point;
//  Unfit::LevenbergMarquardt lm;
//  rc = lm.FindMin(cost_func, lm_min_point);
//  CHECK_EQUAL(0, rc);
//  CHECK_CLOSE(216.433109901728, lm_min_point[0], 1e-3);
//  CHECK_CLOSE(2.00476590272927, lm_min_point[1], 1e-3);
//  CHECK_CLOSE(-2.17034564257086, lm_min_point[2], 1e-3);
//  CHECK_CLOSE(-0.226712906915202, lm_min_point[3], 1e-3);
// }

TEST(GeneticAlgorithm_HelicalValley)  // From LEVMAR
{
  Unfit::GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, -10.0, 10.0);
  ga.bounds.SetBounds(1, -10.0, 10.0);
  ga.bounds.SetBounds(2, -10.0, 10.0);
  ga.options.SetGamma(0.5);
  ga.options.SetMaxIterations(20000);
  ga.options.SetMaxFunctionEvaluations(1000000);
  ga.options.SetCostTolerance(1e-4);
  // Create the cost function
  Unfit::Examples::HelicalValleyCostFunction cost_func;

  // Initial guess
  std::vector<double> min_point {-1.0, 0.0, 0.0};
  // Check the function calculates the correct cost at the initial guess
  auto residual = cost_func(min_point);
  CHECK_CLOSE(2500.0, residual[0], 1e-8);

  // Minimise
  auto t1 = hrclock_t::now();  // Start time
  auto rc = ga.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "HelicalValley (GeneticAlgorithm)";
  std::cout << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1.0, min_point[0], 0.01);
  CHECK_CLOSE(0.0, min_point[1], 0.04);
  CHECK_CLOSE(0.0, min_point[2], 0.06);
}

TEST(GeneticAlgorithm_HodgkinHuxleyBetaN)  // From Unfit 1
{
  Unfit::GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, -100.0, 100.0);
  ga.bounds.SetBounds(1, -100.0, 100.0);
  ga.options.SetCostTolerance(0.0165);
  ga.options.SetGamma(0.5);
  // Read in the experimental data
  std::vector<double> vm;
  std::vector<double> beta_n;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/hh_beta_n_data.txt"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, vm));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, beta_n));
  CHECK_EQUAL(vm.size(), beta_n.size());
  // Create the cost function
  Unfit::Examples::HodgkinHuxleyBetaN cost_func(vm, beta_n);

  // Initial guess
  std::vector<double> min_point = {0.120, -81.0};
  // Check the function calculates the correct cost at the initial guess
  auto residual = cost_func(min_point);
  auto sum_sq_residual = Unfit::SumOfSquares(residual);
  CHECK_CLOSE(0.0179026, sum_sq_residual, 1e-6);

  // Minimise
  auto t1 = hrclock_t::now();  // Start time
  auto rc = ga.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "HodgkinHuxleyBetaN ";
  std::cout << "(GeneticAlgorithm)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.11853534, min_point[0], 1e-3);
  CHECK_CLOSE(-91.2798, min_point[1], 1e-3);

  // Check we can find the minimum from the GA initial guess with NM
  auto nm_min_point = min_point;
  Unfit::NelderMead nm;
  rc = nm.FindMin(cost_func, nm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.118753, nm_min_point[0], 1e-3);
  CHECK_CLOSE(-91.5825, nm_min_point[1], 1e-2);

  // Check we can find the minimum from the GA initial guess with LM
  auto lm_min_point = min_point;
  Unfit::LevenbergMarquardt lm;
  rc = lm.FindMin(cost_func, lm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.118753, lm_min_point[0], 1e-3);
  CHECK_CLOSE(-91.5825, lm_min_point[1], 1e-2);
}

TEST(GeneticAlgorithm_ModifiedRosenbrock)  // From LEVMAR
{
  Unfit::GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, -10.0, 10.0);
  ga.bounds.SetBounds(1, -10.0, 10.0);
  ga.options.SetCostTolerance(10000.0001);
  ga.options.SetGamma(0.5);
  // Create the cost function
  Unfit::Examples::Modroslam cost_func(3);

  // Initial guess
  std::vector<double> min_point {-1.2, 1.0};
  // Check the function calculates the correct cost at the initial guess
  auto residual = cost_func(min_point);
  auto sum_sq_residual = Unfit::SumOfSquares(residual);
  CHECK_CLOSE(10024.2, sum_sq_residual, 1e-1);

  // Minimise
  auto t1 = hrclock_t::now();  // Start time
  auto rc = ga.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "ModifiedRosenbrock ";
  std::cout << "(GeneticAlgorithm)" <<std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1.0, min_point[0], 2e-2);
  CHECK_CLOSE(1.0, min_point[1], 2e-2);

  // Check we can find the minimum from the GA initial guess with NM
  auto nm_min_point = min_point;
  Unfit::NelderMead nm;
  rc = nm.FindMin(cost_func, nm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1.0, nm_min_point[0], 1e-3);
  CHECK_CLOSE(1.0, nm_min_point[1], 1e-3);

  // Check we can find the minimum from the GA initial guess with LM
  auto lm_min_point = min_point;
  Unfit::LevenbergMarquardt lm;
  rc = lm.FindMin(cost_func, lm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1.0, lm_min_point[0], 1e-3);
  CHECK_CLOSE(1.0, lm_min_point[1], 1e-3);
}

TEST(GeneticAlgorithm_NonStationaryMarkov)  // From Unfit 1
{
  Unfit::GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, -100.0, 100.0);
  ga.bounds.SetBounds(1, -100.0, 100.0);
  ga.bounds.SetBounds(2, -10.0, 10.0);
  ga.bounds.SetBounds(3, 0.0, 10.0);
  ga.options.SetCostTolerance(141.0);
  ga.options.SetPopulationSize(30u);
  ga.options.SetGamma(0.25);
  ga.options.SetMaxFunctionEvaluations(200000);
  // Read in the experimental data
  std::vector<std::vector<double>> x;
  x.assign(8, std::vector<double> (20, 0));
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/nsMarkovData.txt"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, x[0]));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, x[1]));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(2, x[2]));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(3, x[3]));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(4, x[4]));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(5, x[5]));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(6, x[6]));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(7, x[7]));
  // Create the cost function
  Unfit::Examples::NonStationaryMarkov cost_func(x);

  // Initial guess
  std::vector<double> min_point {1.0, 10.0, 1.0, 1.0};
  // Check the function calculates the correct cost at the initial guess
  auto residual = cost_func(min_point);
  auto sum_sq_residual = Unfit::SumOfSquares(residual);
  CHECK_CLOSE(2576.29, sum_sq_residual, 1e-2);

  // Minimise
  auto t1 = hrclock_t::now();  // Start time
  auto rc = ga.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "NonStationaryMarkov ";
  std::cout << "(GeneticAlgorithm)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(20.868722, min_point[0], 0.3);
  CHECK_CLOSE(24.065139, min_point[1], 0.05);
  CHECK_CLOSE(0.29473423, min_point[2], 0.005);
  CHECK_CLOSE(0.19128591, fabs(min_point[3]), 0.01);

  // Check we can find the minimum from the GA initial guess with NM
  auto nm_min_point = min_point;
  Unfit::NelderMead nm;
  rc = nm.FindMin(cost_func, nm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(21.243, nm_min_point[0], 1e-3);
  CHECK_CLOSE(23.826, nm_min_point[1], 1e-3);
  CHECK_CLOSE(0.29907, nm_min_point[2], 1e-3);
  CHECK_CLOSE(0.18659, fabs(nm_min_point[3]), 1e-3);

  // Check we can find the minimum from the GA initial guess with LM
  auto lm_min_point = min_point;
  Unfit::LevenbergMarquardt lm;
  rc = lm.FindMin(cost_func, lm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(21.243, lm_min_point[0], 1e-3);
  CHECK_CLOSE(23.826, lm_min_point[1], 1e-3);
  CHECK_CLOSE(0.29907, lm_min_point[2], 1e-3);
  CHECK_CLOSE(0.18659, fabs(lm_min_point[3]), 1e-3);
}

TEST(GeneticAlgorithm_ODE2D)  // From Unfit 1
{
  Unfit::GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, -100.0, 100.0);
  ga.bounds.SetBounds(1, -100.0, 100.0);
  ga.options.SetMaxIterations(100000);
  ga.options.SetMaxFunctionEvaluations(1000000);
  ga.options.SetGamma(0.5);
  // Read in the experimental data
  std::vector<double> t;
  std::vector<double> x;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/ode_data.txt"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, t));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, x));
  auto dt = t[1] - t[0];
  // Create the cost function
  Unfit::Examples::ODE2D cost_func(x, dt);

  // Initial guess
  std::vector<double> min_point = {1.9, 3.1};
  // Check the function calculates the correct cost at the initial guess
  auto residual = cost_func(min_point);
  auto sum_sq_residual = Unfit::SumOfSquares(residual);
  CHECK_CLOSE(5754.83, sum_sq_residual, 1e-2);

  // Minimise
  auto t1 = hrclock_t::now();  // Start time
  auto rc = ga.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "ODE2D ";
  std::cout << "(GeneticAlgorithm)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(2.1925218, min_point[0], 1e-4);
  CHECK_CLOSE(3.4114404, min_point[1], 1e-4);

  // Check we can find the minimum from the GA initial guess with NM
  auto nm_min_point = min_point;
  Unfit::NelderMead nm;
  rc = nm.FindMin(cost_func, nm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(2.21403, nm_min_point[0], 1e-2);
  CHECK_CLOSE(3.32104, nm_min_point[1], 1e-2);

  // Check we can find the minimum from the GA initial guess with LM
  auto lm_min_point = min_point;
  Unfit::LevenbergMarquardt lm;
  rc = lm.FindMin(cost_func, lm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(2.21403, lm_min_point[0], 1e-2);
  CHECK_CLOSE(3.32104, lm_min_point[1], 1e-2);
}

TEST(GeneticAlgorithm_ODE3DVariant)  // From Unfit 1
{
  Unfit::GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, -100.0, 100.0);
  ga.bounds.SetBounds(1, -100.0, 100.0);
  ga.bounds.SetBounds(2, -10.0, 10.0);
  ga.options.SetCostTolerance(2.0);
  ga.options.SetMaxIterations(100000);
  ga.options.SetMaxFunctionEvaluations(1000000);
  ga.options.SetElitism(3u);
  ga.options.SetGamma(0.5);
  // Read in the experimental data
  std::vector<double> t;
  std::vector<double> x;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/ode_data.txt"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, t));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, x));
  auto dt = t[1] - t[0];
  // Create the cost function
  Unfit::Examples::ODE3DVariant cost_func(x, dt);

  // Initial guess
  std::vector<double> min_point = {1.9, 3.1, 6.0};
  // Check the function calculates the correct cost at the initial guess
  auto residual = cost_func(min_point);
  auto sum_sq_residual = Unfit::SumOfSquares(residual);
  CHECK_CLOSE(188319, sum_sq_residual, 1.0);

  // Minimise
  auto t1 = hrclock_t::now();  // Start time
  auto rc = ga.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "ODE3DVariant (GeneticAlgorithm)"
      << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(2.1782233, min_point[0], 1e-4);
  CHECK_CLOSE(3.7736588, min_point[1], 1e-4);
  CHECK_CLOSE(-0.14156534, min_point[2], 1e-4);

  // Check we can find the minimum from the GA initial guess with NM
  auto nm_min_point = min_point;
  Unfit::NelderMead nm;
  rc = nm.FindMin(cost_func, nm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(2.21403, nm_min_point[0], 1e-2);
  CHECK_CLOSE(3.32104, nm_min_point[1], 1e-2);
  CHECK_CLOSE(0.0, nm_min_point[2], 1e-2);

  // Check we can find the minimum from the GA initial guess with LM
  auto lm_min_point = min_point;
  Unfit::LevenbergMarquardt lm;
  rc = lm.FindMin(cost_func, lm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(2.21403, lm_min_point[0], 1e-2);
  CHECK_CLOSE(3.32104, lm_min_point[1], 1e-2);
  CHECK_CLOSE(0.0, lm_min_point[2], 1e-2);
}

TEST(GeneticAlgorithm_Osborne)  // From LEVMAR
{
  Unfit::GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, -10.0, 10.0);
  ga.bounds.SetBounds(1, -10.0, 10.0);
  ga.bounds.SetBounds(2, -10.0, 10.0);
  ga.bounds.SetBounds(3, -10.0, 10.0);
  ga.bounds.SetBounds(4, -10.0, 10.0);
  ga.options.SetPopulationSize(40u);
  ga.options.SetCostTolerance(0.05);
  ga.options.SetMaxFunctionEvaluations(1000000);
  // Set the experimental data
  std::vector<double> x {8.44e-1, 9.08e-1, 9.32e-1, 9.36e-1, 9.25e-1, 9.08e-1,
      8.81e-1, 8.5e-1, 8.18e-1, 7.84e-1, 7.51e-1, 7.18e-1, 6.85e-1, 6.58e-1,
      6.28e-1, 6.03e-1, 5.8e-1, 5.58e-1, 5.38e-1, 5.22e-1, 5.06e-1, 4.9e-1,
      4.78e-1, 4.67e-1, 4.57e-1, 4.48e-1, 4.38e-1, 4.31e-1, 4.24e-1, 4.2e-1,
      4.14e-1, 4.11e-1, 4.06e-1};
  // Create the cost function
  Unfit::Examples::Osborne cost_func(x);

  // Initial guess
  std::vector<double> min_point = {0.5, 1.5, -1.0, 0.01, 0.02};
  // Check the function calculates the correct cost at the initial guess
  auto residual = cost_func(min_point);
  auto sum_sq_residual = Unfit::SumOfSquares(residual);
  CHECK_CLOSE(0.879026, sum_sq_residual, 1e-6);

  // Minimise
  auto t1 = hrclock_t::now();  // Start time
  auto rc = ga.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Osborne (GeneticAlgorithm)";
  std::cout << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.45247469, min_point[0], 1e-3);
  CHECK_CLOSE(0.95389702, min_point[1], 1e-3);
  CHECK_CLOSE(-0.57352404, min_point[2], 1e-3);
  CHECK_CLOSE(0.013624981, min_point[3], 1e-3);
  CHECK_CLOSE(0.053401499, min_point[4], 1e-3);

  // Check we can find the minimum from the GA initial guess with NM
  auto nm_min_point = min_point;
  Unfit::NelderMead nm;
  rc = nm.FindMin(cost_func, nm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.3754, nm_min_point[0], 1e-3);
  CHECK_CLOSE(1.9358, nm_min_point[1], 1e-3);
  CHECK_CLOSE(-1.4647, nm_min_point[2], 1e-3);
  CHECK_CLOSE(0.0129, nm_min_point[3], 1e-3);
  CHECK_CLOSE(0.0221, nm_min_point[4], 1e-3);

  // Check we can find the minimum from the GA initial guess with LM
  auto lm_min_point = min_point;
  Unfit::LevenbergMarquardt lm;
  rc = lm.FindMin(cost_func, lm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.3754, lm_min_point[0], 1e-3);
  CHECK_CLOSE(1.9358, lm_min_point[1], 1e-3);
  CHECK_CLOSE(-1.4647, lm_min_point[2], 1e-3);
  CHECK_CLOSE(0.0129, lm_min_point[3], 1e-3);
  CHECK_CLOSE(0.0221, lm_min_point[4], 1e-3);
}

TEST(GeneticAlgorithm_Parabolic)
{
  Unfit::GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, -100.0, 100.0);
  ga.bounds.SetBounds(1, -100.0, 100.0);
  ga.bounds.SetBounds(2, -100.0, 100.0);
  ga.options.SetCostTolerance(25.0);
  ga.options.SetGamma(0.5);
  // Read in the experimental data
  std::vector<double> x;
  std::vector<double> y;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/parabolic_data.txt"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, x));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, y));
  CHECK_EQUAL(x.size(), y.size());
  // Create the cost function
  Unfit::Examples::Parabolic cost_func(x, y);

  // Initial guess
  std::vector<double> min_point = {1.0, 1.0, 1.0};
  // Check the function calculates the correct cost at the initial guess
  auto residual = cost_func(min_point);
  auto sum_sq_residual = Unfit::SumOfSquares(residual);
  CHECK_CLOSE(6850.23, sum_sq_residual, 1e-2);

  // Minimise
  auto t1 = hrclock_t::now();  // Start time
  auto rc = ga.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Parabolic (GeneticAlgorithm)";
  std::cout << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(-21.047296, min_point[0], 1e-2);
  CHECK_CLOSE(60.745038, min_point[1], 1e-2);
  CHECK_CLOSE(-9.3505343, min_point[2], 1e-2);

  // Check we can find the minimum from the GA initial guess with NM
  auto nm_min_point = min_point;
  Unfit::NelderMead nm;
  rc = nm.FindMin(cost_func, nm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(-21.1144, nm_min_point[0], 1e-2);
  CHECK_CLOSE(61.0351, nm_min_point[1], 1e-2);
  CHECK_CLOSE(-9.61397, nm_min_point[2], 1e-2);

  // Check we can find the minimum from the GA initial guess with LM
  auto lm_min_point = min_point;
  Unfit::LevenbergMarquardt lm;
  rc = lm.FindMin(cost_func, lm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(-21.1144, lm_min_point[0], 1e-2);
  CHECK_CLOSE(61.0351, lm_min_point[1], 1e-2);
  CHECK_CLOSE(-9.61397, lm_min_point[2], 1e-2);
}

TEST(GeneticAlgorithm_Powell)  // From LEVMAR
{
  Unfit::GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, -100.0, 100.0);
  ga.bounds.SetBounds(1, -100.0, 100.0);
  ga.options.SetElitism(2u);
  ga.options.SetGamma(0.5);
  // Create the cost function
  Unfit::Examples::Powell cost_func(4);

  // Initial guess
  std::vector<double> min_point {3.0, 1.0};
  // Check the function calculates the correct cost at the initial guess
  auto residual = cost_func(min_point);
  auto sum_sq_residual = Unfit::SumOfSquares(residual);
  CHECK_CLOSE(290.724, sum_sq_residual, 1e-3);

  // Minimise
  auto t1 = hrclock_t::now();  // Start time
  auto rc = ga.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Powell (GeneticAlgorithm)";
  std::cout << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 0.1);
  CHECK_CLOSE(0.0, min_point[1], 0.1);

  // Check we can find the minimum from the GA initial guess with NM
  auto nm_min_point = min_point;
  Unfit::NelderMead nm;
  rc = nm.FindMin(cost_func, nm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, nm_min_point[0], 1e-3);
  CHECK_CLOSE(0.0, nm_min_point[1], 1e-2);

  // Check we can find the minimum from the GA initial guess with LM
  auto lm_min_point = min_point;
  Unfit::LevenbergMarquardt lm;
  rc = lm.FindMin(cost_func, lm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, lm_min_point[0], 1e-3);
  CHECK_CLOSE(0.0, lm_min_point[1], 1e-3);
}

TEST(GeneticAlgorithm_Rosenbrock)  // From LEVMAR
{
  Unfit::GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, -10.0, 10.0);
  ga.bounds.SetBounds(1, -10.0, 10.0);
  ga.options.SetGamma(0.5);
  // Create the cost function
  Unfit::Examples::Rosenbrock cost_func(2);

  // Initial guess
  std::vector<double> min_point {-1.2, 1.0};
  // Check the function calculates the correct cost at the initial guess
  auto residual = cost_func(min_point);
  auto sum_sq_residual = Unfit::SumOfSquares(residual);
  CHECK_CLOSE(1266.86, sum_sq_residual, 1e-2);

  // Minimise
  auto t1 = hrclock_t::now();  // Start time
  auto rc = ga.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Rosenbrock (GeneticAlgorithm)";
  std::cout << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1.0, min_point[0], 0.01);
  CHECK_CLOSE(1.0, min_point[1], 0.01);
}

TEST(GeneticAlgorithm_ThreeNaMarkov)  // From Unfit 1
{
  Unfit::GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, 0.0, 10.0);
  ga.bounds.SetBounds(1, 0.0, 10.0);
  ga.bounds.SetBounds(2, 0.0, 10.0);
  ga.bounds.SetBounds(3, 0.0, 10.0);
  ga.options.SetCostTolerance(0.006);
  ga.options.SetGamma(0.25);
  // Read in the experimental data
  std::vector<double> t;
  std::vector<double> po;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/ThreeNaMarkov.csv"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, t));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, po));
  auto dt = t[1] - t[0];
  // Create the cost function
  Unfit::Examples::ThreeNaMarkov cost_func(po, dt);

  // Initial guess
  std::vector<double> min_point = {1.0, 0.01, 1.0, 0.01};
  // Check the function calculates the correct cost at the initial guess
  auto residual = cost_func(min_point);
  auto sum_sq_residual = Unfit::SumOfSquares(residual);
  CHECK_CLOSE(1.09713, sum_sq_residual, 1e-4);

  // Minimise
  auto t1 = hrclock_t::now();  // Start time
  auto rc = ga.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "ThreeNaMarkov (GeneticAlgorithm)";
  std::cout << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.060244368, min_point[0], 1e-4);
  CHECK_CLOSE(0.10951119, min_point[1], 1e-4);
  CHECK_CLOSE(0.31749085, min_point[2], 1e-4);
  CHECK_CLOSE(0.0022059417, min_point[3], 1e-4);

  // Check we can find the minimum from the GA initial guess with NM
  auto nm_min_point = min_point;
  Unfit::NelderMead nm;
  rc = nm.FindMin(cost_func, nm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0494253, nm_min_point[0], 1e-6);
  CHECK_CLOSE(0.0283517, nm_min_point[1], 1e-5);
  CHECK_CLOSE(0.321102, nm_min_point[2], 1e-6);
  CHECK_CLOSE(0.00174363, nm_min_point[3], 1e-6);

  // Check we can find the minimum from the GA initial guess with LM
  auto lm_min_point = min_point;
  Unfit::LevenbergMarquardt lm;
  rc = lm.FindMin(cost_func, lm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0494253, lm_min_point[0], 1e-6);
  CHECK_CLOSE(0.0283517, lm_min_point[1], 1e-6);
  CHECK_CLOSE(0.321102, lm_min_point[2], 1e-6);
  CHECK_CLOSE(0.00174363, lm_min_point[3], 1e-6);
}

TEST(GeneticAlgorithm_Woods)  // From LEVMAR
{
  Unfit::GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, -10.0, 10.0);
  ga.bounds.SetBounds(1, -10.0, 10.0);
  ga.bounds.SetBounds(2, -10.0, 10.0);
  ga.bounds.SetBounds(3, -10.0, 10.0);
  ga.options.SetElitism(4u);
  ga.options.SetCostTolerance(0.01);
  ga.options.SetGamma(0.25);
  // Create the cost function
  Unfit::Examples::Woods cost_func;

  // Initial guess
  std::vector<double> min_point {-3.0, -1.0, -3.0, -1.0};
  // Check the function calculates the correct cost at the initial guess
  auto residual = cost_func(min_point);
  auto sum_sq_residual = Unfit::SumOfSquares(residual);
  CHECK_CLOSE(19192, sum_sq_residual, 1.0);

  // Minimise
  auto t1 = hrclock_t::now();  // Start time
  auto rc = ga.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Woods (GeneticAlgorithm)";
  std::cout << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1, min_point[0], 0.05);
  CHECK_CLOSE(1, min_point[1], 0.05);
  CHECK_CLOSE(1, min_point[2], 0.05);
  CHECK_CLOSE(1, min_point[3], 0.05);

  // Check we can find the minimum from the GA initial guess with NM
  auto nm_min_point = min_point;
  Unfit::NelderMead nm;
  rc = nm.FindMin(cost_func, nm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1, nm_min_point[0], 1e-3);
  CHECK_CLOSE(1, nm_min_point[1], 1e-3);
  CHECK_CLOSE(1, nm_min_point[2], 1e-3);
  CHECK_CLOSE(1, nm_min_point[3], 1e-3);

  // Check we can find the minimum from the GA initial guess with LM
  auto lm_min_point = min_point;
  Unfit::LevenbergMarquardt lm;
  rc = lm.FindMin(cost_func, lm_min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1, lm_min_point[0], 1e-3);
  CHECK_CLOSE(1, lm_min_point[1], 1e-3);
  CHECK_CLOSE(1, lm_min_point[2], 1e-3);
  CHECK_CLOSE(1, lm_min_point[3], 1e-3);
}
}  // suite ExamplesGeneticAlgorithm
