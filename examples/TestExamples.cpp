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
#include <vector>
#include "DataFileReader.hpp"
#include "GenericNDCostFunction.hpp"
#include "HelicalValleyCostFunction.hpp"
#include "HodgkinHuxleyBetaNModel.hpp"
#include "ParabolicModel.hpp"
#include "SimpleParabolicCostFunction.hpp"
#include "Unfit.hpp"
#include "UnitTest++.h"

// Introduction:
//   This file contains a number of examples showing different ways of doing
// things with Unfit, whether you are doing an optimization, fitting data to
// a function, or something else. In parallel with this, you should look through
// the models/cost functions that are called so you can see how they are
// implemented.
//   We have chosen to provide all of these examples in terms of tests, but when
// you are writing your own models (or cost functions) it is up to you whether
// you use UnitTest++ or not. To that end, the keywords that show something
// UnitTest++ is doing within a test are those that say CHECK, CHECK_EQUAL, or
// CHECK_CLOSE. Independent of whether or not you chose to use UnitTest++, we
// strongly encourage you to check the return code of anything that returns one,
// be it the DataFileReader or the optimizers. Only by checking the optimizer
// return code will you know if it converged.

SUITE(UnfitExamples)
{
TEST(MinimiseASimpleFunctionNoData)
{
  // Here we have no data, no model, just a function to minimise, so instead
  // of creating a model, we have to create a cost function directly.
  Unfit::Examples::SimpleParabolicCostFunction sp_cost;
  // Choose an initial guess for our parameter (c0)
  std::vector<double> c {1.0};

  // Try the Nelder-Mead simplex optimisation algorithm
  Unfit::NelderMead nm_opt;
  // Pass in the cost function and our initial guess for c
  auto rc = nm_opt.FindMin(sp_cost, c);
  // Check the optimiser came back with a converged solution
  CHECK_EQUAL(0, rc);
  // Check the solution (we know the exact solution in this case)
  CHECK_CLOSE(0.0, c[0], 1e-3);
}

TEST(MinimiseATrickyFunctionNoData)
{
  // Here we have no data, no model, just a function to minimise, so instead
  // of creating a model, we have to create a cost function directly.
  Unfit::Examples::HelicalValleyCostFunction hv_cost;
  // Choose an initial guess for our parameters (c0, c1, c2)
  std::vector<double> c {0.5, 0.5, 0.5};

  // Try the Nelder-Mead simplex optimisation algorithm
  Unfit::NelderMead nm_opt;
  // Pass in the cost function and our initial guess for c
  auto rc = nm_opt.FindMin(hv_cost, c);
  // Check the optimiser came back with a converged solution
  CHECK_EQUAL(0, rc);
  // Check the solution (we know the exact solution in this case)
  CHECK_CLOSE(1.0, c[0], 1e-3);
  CHECK_CLOSE(0.0, c[1], 1e-3);
  CHECK_CLOSE(0.0, c[2], 1e-3);
}

TEST(ReadDataFromFileAndOptimise)
{
  // Read in the experimental data; first create the DataFileReader
  Unfit::DataFileReader<double> dfr;
  // Get the reader to open the file and import the data. Check the return code
  auto rc = dfr.ReadFile("examples/data/parabolic_data.txt");
  CHECK_EQUAL(0u, rc);
  // Put the data into x and y. We must set the number of vectors in x.
  std::vector<std::vector<double>> x(1);
  rc = dfr.RetrieveColumn(0, x[0]);
  CHECK_EQUAL(0u, rc);
  std::vector<double> y;
  rc = dfr.RetrieveColumn(1, y);
  CHECK_EQUAL(0u, rc);
  // A quick check x and y have the same number of data points
  CHECK_EQUAL(x[0].size(), y.size());

  // Make ourselves an object of the desired model type, in this case a parabola
  Unfit::Examples::ParabolicModel parabola;
  // We want to fit the data by optimizing our parameters, so we need a function
  // to calculate the cost of our model given this data set (x, y)
  Unfit::GenericNDCostFunction parabola_cost(parabola, x, y);
  // Choose an initial guess for our parameters (c0, c1, c2)
  std::vector<double> c {1.0, 1.0, 1.0};

  // Try the Nelder-Mead simplex optimisation algorithm
  Unfit::NelderMead nm_opt;
  // Pass in the cost function and our initial guess for c
  auto rc2 = nm_opt.FindMin(parabola_cost, c);
  // Check the optimiser came back with a converged solution
  CHECK_EQUAL(0, rc2);
}

TEST(TwoDimensionTwoParameterModel)
{
  // Read in the experimental data
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/parabolic_data.txt"));
  std::vector<std::vector<double>> x(1);
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, x[0]));
  std::vector<double> y;
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, y));
  CHECK_EQUAL(x[0].size(), y.size());

  // Make ourselves an object of the desired model type
  Unfit::Examples::HodgkinHuxleyBetaNModel hh_beta_n;
  // We want to fit the data by optimizing our parameters, so we need a function
  // to calculate the cost of our model given this data set (x, y)
  Unfit::GenericNDCostFunction hh_beta_n_cost(hh_beta_n, x, y);
  // Choose an initial guess for our parameters (c0, c1)
  std::vector<double> c {1.0, 1.0};

  // Try the Nelder-Mead simplex optimisation algorithm
  Unfit::NelderMead nm_opt;
  // Pass in the cost function and our initial guess for c
  auto rc = nm_opt.FindMin(hh_beta_n_cost, c);
  // Check the optimiser came back with a converged solution
  CHECK_EQUAL(0, rc);
}

TEST(TwoDimensionThreeParameterModel)
{
  // Set the experimental data
  std::vector<std::vector<double>> x {{0.498531, 0.622145, 0.746551, 0.899687,
      0.995019, 1.24803, 1.49695, 1.7464, 1.86737, 1.92478, 2.07206, 2.12789,
      2.23212}};
  std::vector<double> y {17.0676, 20.0356, 24.0914, 27.5963, 28.9598, 31.9736,
      34.6866, 33.7931, 31.9415, 30.6897, 28.3853, 23.9687, 18.5146};
  // Make ourselves an object of the desired model type, in this case a parabola
  Unfit::Examples::ParabolicModel parabola;
  // We want to fit the data by optimizing our parameters, so we need a function
  // to calculate the cost of our model given this data set (x, y)
  Unfit::GenericNDCostFunction parabola_cost(parabola, x, y);
  // Choose an initial guess for our parameters (c0, c1, c2)
  std::vector<double> c {1.0, 1.0, 1.0};

  // Try the Nelder-Mead simplex optimisation algorithm
  Unfit::NelderMead nm_opt;
  // Pass in the cost function and our initial guess for c
  auto rc = nm_opt.FindMin(parabola_cost, c);
  // Check the optimiser came back with a converged solution
  CHECK_EQUAL(0, rc);
}


}  // suite UnfitExamples

//TEST(LevenbergMarquardt_FindMinParabolicSetBounds)
//{
//  Unfit::LevenbergMarquardt test;
//  test.bounds.SetBounds({0.0, 0.0, 0.0}, {0.9, 1.9, 2.9});
//  // test.SetStepSize(0.001);
//  test.options.SetMaxIterations(30000);
//  test.options.SetCostTolerance(1e-8);
//  std::vector<double> x {-5.0, -4.0, -3.0, -2.0, -1.0, 0.0,
//                         1.0, 2.0, 3.0, 4.0, 5.0};
//  std::vector<double> y {18.0, 11.0, 6.0, 3.0, 2.0, 3.0,
//                         6.0, 11.0, 18.0, 27.0, 38.0};
//  Unfit::Examples::Parabolic para_eq(x, y);
//  std::vector<double> init_guess {0, 0, 0};
//  int rc = test.FindMin(para_eq, init_guess);
//  CHECK_EQUAL(0, rc);
//  CHECK_CLOSE(0.9, init_guess[0], 1e-8);
//  CHECK_CLOSE(1.9, init_guess[1], 1e-4);
//  CHECK_CLOSE(2.9, init_guess[2], 1e-4);
//
//  // Check the function calculates the correct cost at the initial guess
//  std::vector<double> residual = para_eq({0.9, 1.9, 2.9});
//  double sum {0.0};
//  for (unsigned i = 0; i < residual.size(); ++i) sum += residual[i]*residual[i];
//  CHECK_CLOSE(22.99, sum, 1e-4);
//
//  residual = para_eq(init_guess);
//  sum = 0.0;
//  for (unsigned i = 0; i < residual.size(); ++i) sum += residual[i]*residual[i];
//  CHECK_CLOSE(22.99, sum, 1e-4);
//}
//
//TEST(LevenbergMarquardt_FindMinGaussianSetBounds)
//{
//  Unfit::LevenbergMarquardt test;
//  test.bounds.SetBounds({1.1, 2.1, 3.1, 0.1}, {10.0, 10.0, 10.0, 10.0});
//  test.options.SetMaxIterations(30000);
//  test.options.SetCostTolerance(1e-8);
//  std::vector<double> x {-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0,
//                         5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
//  std::vector<double> y {0.0657285286165, 0.135335283236, 0.249352208777,
//                         0.411112290507,  0.606530659712, 0.800737402916,
//                         0.945959468906,  1.0, 0.945959468906, 0.800737402916,
//                         0.606530659712, 0.411112290507, 0.249352208777,
//                         0.135335283236, 0.0657285286165, 0.0285655007845};
//  Unfit::Examples::GaussianEquation gauss_eq(x, y);
//  std::vector<double> init_guess {1.5, 2.7, 3.3, 7.0};
//  int rc = test.FindMin(gauss_eq, init_guess);
//  CHECK_EQUAL(0, rc);
//  CHECK_CLOSE(1.1, init_guess[0], 1e-4);
//  CHECK_CLOSE(2.1, init_guess[1], 1e-4);
//  CHECK_CLOSE(3.1, init_guess[2], 1e-4);
//  CHECK_CLOSE(0.1, init_guess[3], 1e-4);
//  std::vector<double> residual = gauss_eq({1.1, 2.1, 3.1, 0.1});
//  double sum {0.0};
//  for (unsigned i = 0; i < residual.size(); ++i) sum += residual[i]*residual[i];
//  CHECK_CLOSE(0.441077, sum, 1e-4);
//
//  residual = gauss_eq(init_guess);
//  sum = 0.0;
//  for (unsigned i = 0; i < residual.size(); ++i) sum += residual[i]*residual[i];
//  CHECK_CLOSE(0.441077, sum, 1e-4);
//}
//
//TEST(LevenbergMarquardt_SampleHimmelblauSetBounds)
//{
//  Unfit::LevenbergMarquardt object;
//  object.options.SetMaxIterations(300000);
//  object.bounds.SetBounds({-10.0, -10.0}, {2.9, 1.9});
//  object.options.SetCostTolerance(1e-24);
//  Unfit::Examples::Himmelblau cost_func;
//  // Initial guess
//  // std::vector<double> min_point = {2.0, 1.0};
//  std::vector<double> min_point = {-0.3, -2.0};
//  // Check the function calculates the correct cost at the initial guess
//// std::vector<double> residual = cost_func(min_point);
//// double sum {0.0};
//// for(unsigned i = 0; i < residual.size(); ++i) sum += residual[i]*residual[i];
////  CHECK_CLOSE(512, sum, 1e-4);
//  // Minimise
//  int rc = object.FindMin(cost_func, min_point);
//  // Check the result matches what we expect
//  CHECK_EQUAL(0, rc);
//  CHECK_CLOSE(-3.77931, min_point[0], 1e-3);
//  CHECK_CLOSE(-3.283186, min_point[1], 1e-3);
//}
//
//TEST(LevenbergMarquardt_FindMinGaussian)
//{
//  Unfit::LevenbergMarquardt test;
//  test.options.SetMaxIterations(3000);
//  std::vector<double> x {-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0,
//                         5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
//  std::vector<double> y {0.0657285286165, 0.135335283236, 0.249352208777,
//                         0.411112290507,  0.606530659712, 0.800737402916,
//                         0.945959468906,  1.0, 0.945959468906, 0.800737402916,
//                         0.606530659712, 0.411112290507, 0.249352208777,
//                         0.135335283236, 0.0657285286165, 0.0285655007845};
//  Unfit::Examples::GaussianEquation gauss_eq(x, y);
//  std::vector<double> init_guess {1.5, 2.7, 3.3, 7.0};
//  auto rc =   test.FindMin(gauss_eq, init_guess);
//  CHECK_EQUAL(0, rc);
//  CHECK_CLOSE(1.0, init_guess[0], 1e-4);
//  CHECK_CLOSE(2.0, init_guess[1], 1e-4);
//  CHECK_CLOSE(3.0, init_guess[2], 1e-4);
//  CHECK_CLOSE(0.0, init_guess[3], 1e-4);
//}
//
//TEST(LevenbergMarquardt_FindMinExponential)
//{
//  Unfit::LevenbergMarquardt test;
//  test.options.SetMaxIterations(3000);
//  std::vector<double> x {-10.0, -9.5, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0,
//                         -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5,
//                         -1.0, -0.5, 0.0, 0.5};
//  std::vector<double> y {-0.0111089965382, -0.0132335500965, -0.0157644164848,
//                         -0.0187793014946, -0.0223707718561, -0.0266490973363,
//                         -0.0317456363780, -0.0378168692293, -0.0450492023935,
//                         -0.0536646919127, -0.0639278612067, -0.0761538227986,
//                         -0.0907179532894, -0.108067418634,  -0.128734903587,
//                         -0.153354966844,  -0.182683524052,  -0.217621056865,
//                         -0.259240260645,  -0.30881897968,   -0.367879441171,
//                         -0.438234992464};
//  Unfit::Examples::Exponential exponential_func(x, y);
////   initial guess is very sensitive to 2nd and 3rd element
//  std::vector<double> init_guess {1.0, 0.2, -0.2};
//  auto rc = test.FindMin(exponential_func, init_guess);
//  CHECK_EQUAL(0, rc);
//  CHECK_CLOSE(0.0, init_guess[0], 1e-5);
//  CHECK_CLOSE(0.35, init_guess[1], 1e-5);
//  CHECK_CLOSE(-1.0, init_guess[2], 1e-5);
//}
//
//TEST(LevenbergMarquardt_SampleHimmelblau)
//{
//  Unfit::LevenbergMarquardt object;
//  Unfit::Examples::Himmelblau cost_func;
//  object.options.SetCostTolerance(1e-24);
//  // Initial guess
//  std::vector<double> min_point = {2.0, 3.0};
//  // Check the function calculates the correct cost at the initial guess
//  std::vector<double> residual = cost_func(min_point);
//  double sum {0.0};
//  for (unsigned i = 0; i < residual.size(); ++i) sum += residual[i]*residual[i];
//  CHECK_CLOSE(512, sum, 1e-4);
//  // Minimise
//  auto rc = object.FindMin(cost_func, min_point);
//  // Check the result matches what we expect
//  CHECK_EQUAL(0, rc);
//  CHECK_CLOSE(3.0, min_point[0], 1e-3);
//  CHECK_CLOSE(2.0, min_point[1], 1e-3);  // cant find yet
//}
//
//TEST(LevenbergMarquardt_SampleHimmelblau1)
//{
//  Unfit::LevenbergMarquardt object;
//  Unfit::Examples::Himmelblau cost_func;
//  object.options.SetCostTolerance(1e-16);
//  // Initial guess
//  std::vector<double> min_point = {-2.0, -3.0};
//  // Check the function calculates the correct cost at the initial guess
//  std::vector<double> residual = cost_func(min_point);
//  double sum {0.0};
//  for (unsigned i = 0; i < residual.size(); ++i) sum += residual[i]*residual[i];
//  CHECK_CLOSE(10000, sum, 1e-4);
//  // Minimise
//  auto rc = object.FindMin(cost_func, min_point);
//  // Check the result matches what we expect
//  CHECK_EQUAL(0, rc);
//  CHECK_CLOSE(-3.77931, min_point[0], 1e-4);
//  CHECK_CLOSE(-3.28319, min_point[1], 1e-4);  // cant find yet
//}
