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
#include "GeneticSwitch.hpp"
#include "HelicalValleyCostFunction.hpp"
#include "HodgkinHuxleyBetaN.hpp"
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

SUITE(ExamplesNelderMead)
{
TEST(NelderMead_CardiacAlphaN)  // From Unfit 1
{
  Unfit::NelderMead nm;
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
  auto rc = nm.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "CardiacAlphaN (NelderMead)";
  std::cout << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0731685, min_point[0], 1e-4);
  CHECK_CLOSE(-0.033983, min_point[1], 1e-4);
  CHECK_CLOSE(3.10956, min_point[2], 1e-4);
}

TEST(NelderMead_Exponential)
{
  Unfit::NelderMead nm;
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
  auto rc = nm.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Exponential (NelderMead)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(3.66696, min_point[0], 1e-2);
  CHECK_CLOSE(0.0885905, min_point[1], 1e-2);
  CHECK_CLOSE(3.06235, min_point[2], 1e-2);
}

TEST(NelderMead_Gaussian)
{
  Unfit::NelderMead nm;
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
  auto rc = nm.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Gaussian (NelderMead)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(12.4144, min_point[0], 1e-4);
  CHECK_CLOSE(22.7355, min_point[1], 1e-4);
  CHECK_CLOSE(4.0337, min_point[2], 1e-4);
  CHECK_CLOSE(12.5486, min_point[3], 1e-4);
}

TEST(NelderMead_GeneticSwitch)  // From Unfit 1
{
  Unfit::NelderMead nm;
  // Read in the experimental data
  std::vector<double> t;
  std::vector<double> z_a;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/training_data_zA.txt"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, t));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, z_a));
  auto dt = t[1] - t[0];
  // Create the cost function
  Unfit::Examples::GeneticSwitch cost_func(z_a, dt);

  // Initial guess
  std::vector<double> min_point = {90.0, 60.0, -1.0, -0.1};
  // Check the function calculates the correct cost at the initial guess
  auto residual = cost_func(min_point);
  auto sum_sq_residual = Unfit::SumOfSquares(residual);
  CHECK_CLOSE(1.12863e7, sum_sq_residual, 1e2);

  // Minimise
  auto t1 = hrclock_t::now();  // Start time
  auto rc = nm.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "GeneticSwitch (NelderMead)";
  std::cout << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(216.433109901728, min_point[0], 1e-3);
  CHECK_CLOSE(2.00476590272927, min_point[1], 1e-3);
  CHECK_CLOSE(-2.17034564257086, min_point[2], 1e-3);
  CHECK_CLOSE(-0.226712906915202, min_point[3], 1e-3);
}

TEST(NelderMead_HelicalValley)  // From LEVMAR
{
  Unfit::NelderMead nm;
  // Create the cost function
  Unfit::Examples::HelicalValleyCostFunction cost_func;

  // Initial guess
  std::vector<double> min_point {-1.0, 0.0, 0.0};
  // Check the function calculates the correct cost at the initial guess
  auto residual = cost_func(min_point);
  CHECK_CLOSE(2500.0, residual[0], 1e-8);

  // Minimise
  auto t1 = hrclock_t::now();  // Start time
  auto rc = nm.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "HelicalValley (NelderMead)";
  std::cout << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1.0, min_point[0], 1e-3);
  CHECK_CLOSE(0.0, min_point[1], 1e-3);
  CHECK_CLOSE(0.0, min_point[2], 1e-3);
}

TEST(NelderMead_HodgkinHuxleyBetaN)  // From Unfit 1
{
  Unfit::NelderMead nm;
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
  auto rc = nm.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "HodgkinHuxleyBetaN (NelderMead)";
  std::cout << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.118753, min_point[0], 1e-3);
  CHECK_CLOSE(-91.5825, min_point[1], 1e-3);

  // Now, try choosing an alternative initial guess
  min_point = {-100.0, -100.0};

  // Minimise
  t1 = hrclock_t::now();  // Start time
  rc = nm.FindMin(cost_func, min_point);
  t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "HodgkinHuxleyBetaN_2 (NelderMead)";
  std::cout << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.118753, min_point[0], 1e-3);
  CHECK_CLOSE(-91.5825, min_point[1], 1e-3);
}

TEST(NelderMead_ModifiedRosenbrock)  // From LEVMAR
{
  Unfit::NelderMead nm;
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
  auto rc = nm.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "ModifiedRosenbrock (NelderMead)";
  std::cout << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1.0, min_point[0], 1e-3);
  CHECK_CLOSE(1.0, min_point[1], 1e-3);
}

TEST(NelderMead_NonStationaryMarkov)  // From Unfit 1
{
  Unfit::NelderMead nm;
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
  auto rc = nm.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "NonStationaryMarkov (NelderMead)";
  std::cout << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(21.243, min_point[0], 1e-3);
  CHECK_CLOSE(23.826, min_point[1], 1e-3);
  CHECK_CLOSE(0.29907, min_point[2], 1e-3);
  CHECK_CLOSE(0.18659, fabs(min_point[3]), 1e-3);
}

TEST(NelderMead_ODE2D)  // From Unfit 1
{
  Unfit::NelderMead nm;
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
  auto rc = nm.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "ODE2D (NelderMead)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(2.21403, min_point[0], 1e-2);
  CHECK_CLOSE(3.32104, min_point[1], 1e-2);
}

TEST(NelderMead_ODE3DVariant)  // From Unfit 1
{
  Unfit::NelderMead nm;
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
  auto rc = nm.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "ODE3DVariant (NelderMead)"
      << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(2.21403, min_point[0], 1e-2);
  CHECK_CLOSE(3.32104, min_point[1], 1e-2);
  CHECK_CLOSE(0.0, min_point[2], 1e-2);
}

TEST(NelderMead_Osborne)  // From LEVMAR
{
  Unfit::NelderMead nm;
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
  auto rc = nm.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Osborne (NelderMead)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.3754, min_point[0], 1e-3);
  CHECK_CLOSE(1.9358, min_point[1], 1e-3);
  CHECK_CLOSE(-1.4647, min_point[2], 1e-3);
  CHECK_CLOSE(0.0129, min_point[3], 1e-3);
  CHECK_CLOSE(0.0221, min_point[4], 1e-3);
}

TEST(NelderMead_Parabolic)
{
  Unfit::NelderMead nm;
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
  auto rc = nm.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Parabolic (NelderMead)";
  std::cout << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(-21.1144, min_point[0], 1e-2);
  CHECK_CLOSE(61.0351, min_point[1], 1e-2);
  CHECK_CLOSE(-9.61397, min_point[2], 1e-2);
}

TEST(NelderMead_Powell)  // From LEVMAR
{
  Unfit::NelderMead nm;
  nm.options.SetCostTolerance(1.0e-16);
  nm.options.SetGeometricTolerance(1.0e-6);
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
  auto rc = nm.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Powell (NelderMead)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-3);
  CHECK_CLOSE(0.0, min_point[1], 1e-3);
}

TEST(NelderMead_Rosenbrock)  // From LEVMAR
{
  Unfit::NelderMead nm;
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
  auto rc = nm.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Rosenbrock (NelderMead)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1.0, min_point[0], 1e-3);
  CHECK_CLOSE(1.0, min_point[1], 1e-3);
}

TEST(NelderMead_ThreeNaMarkov)  // From Unfit 1
{
  Unfit::NelderMead nm;
  nm.options.SetGeometricTolerance(1.0e-8);
  nm.options.SetDegenerateTolerance(1.0e-16);
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
  auto rc = nm.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "ThreeNaMarkov (NelderMead)";
  std::cout << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0494253, min_point[0], 1e-6);
  CHECK_CLOSE(0.0283517, min_point[1], 1e-6);
  CHECK_CLOSE(0.321102, min_point[2], 1e-6);
  CHECK_CLOSE(0.00174363, min_point[3], 1e-6);
}

TEST(NelderMead_Woods)  // From LEVMAR
{
  Unfit::NelderMead nm;
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
  auto rc = nm.FindMin(cost_func, min_point);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Woods (NelderMead)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1, min_point[0], 1e-3);
  CHECK_CLOSE(1, min_point[1], 1e-3);
  CHECK_CLOSE(1, min_point[2], 1e-3);
  CHECK_CLOSE(1, min_point[3], 1e-3);
}
}  // suite ExamplesNelderMead
