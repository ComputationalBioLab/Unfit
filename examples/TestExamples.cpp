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
#include "Exponential.hpp"
#include "GaussianEquation.hpp"
#include "Himmelblau.hpp"
#include "LevenbergMarquardt.hpp"
#include "Parabolic.hpp"
#include "UnitTest++.h"

TEST(LevenbergMarquardt_FindMinParabolicSetBounds)
{
  Unfit::LevenbergMarquardt test;
  test.bounds.SetBounds({0.0, 0.0, 0.0}, {0.9, 1.9, 2.9});
  // test.SetStepSize(0.001);
  test.options.SetMaxIterations(30000);
  test.options.SetCostTolerance(1e-8);
  std::vector<double> x {-5.0, -4.0, -3.0, -2.0, -1.0, 0.0,
                         1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> y {18.0, 11.0, 6.0, 3.0, 2.0, 3.0,
                         6.0, 11.0, 18.0, 27.0, 38.0};
  Unfit::Examples::Parabolic para_eq(x, y);
  std::vector<double> init_guess {0, 0, 0};
  int rc = test.FindMin(para_eq, init_guess);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.9, init_guess[0], 1e-8);
  CHECK_CLOSE(1.9, init_guess[1], 1e-4);
  CHECK_CLOSE(2.9, init_guess[2], 1e-4);

  // Check the function calculates the correct cost at the initial guess
  std::vector<double> residual = para_eq({0.9, 1.9, 2.9});
  double sum {0.0};
  for (unsigned i = 0; i < residual.size(); ++i) sum += residual[i]*residual[i];
  CHECK_CLOSE(22.99, sum, 1e-4);

  residual = para_eq(init_guess);
  sum = 0.0;
  for (unsigned i = 0; i < residual.size(); ++i) sum += residual[i]*residual[i];
  CHECK_CLOSE(22.99, sum, 1e-4);
}

TEST(LevenbergMarquardt_FindMinGaussianSetBounds)
{
  Unfit::LevenbergMarquardt test;
  test.bounds.SetBounds({1.1, 2.1, 3.1, 0.1}, {10.0, 10.0, 10.0, 10.0});
  test.options.SetMaxIterations(30000);
  test.options.SetCostTolerance(1e-8);
  std::vector<double> x {-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0,
                         5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  std::vector<double> y {0.0657285286165, 0.135335283236, 0.249352208777,
                         0.411112290507,  0.606530659712, 0.800737402916,
                         0.945959468906,  1.0, 0.945959468906, 0.800737402916,
                         0.606530659712, 0.411112290507, 0.249352208777,
                         0.135335283236, 0.0657285286165, 0.0285655007845};
  Unfit::Examples::GaussianEquation gauss_eq(x, y);
  std::vector<double> init_guess {1.5, 2.7, 3.3, 7.0};
  int rc = test.FindMin(gauss_eq, init_guess);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1.1, init_guess[0], 1e-4);
  CHECK_CLOSE(2.1, init_guess[1], 1e-4);
  CHECK_CLOSE(3.1, init_guess[2], 1e-4);
  CHECK_CLOSE(0.1, init_guess[3], 1e-4);
  std::vector<double> residual = gauss_eq({1.1, 2.1, 3.1, 0.1});
  double sum {0.0};
  for (unsigned i = 0; i < residual.size(); ++i) sum += residual[i]*residual[i];
  CHECK_CLOSE(0.441077, sum, 1e-4);

  residual = gauss_eq(init_guess);
  sum = 0.0;
  for (unsigned i = 0; i < residual.size(); ++i) sum += residual[i]*residual[i];
  CHECK_CLOSE(0.441077, sum, 1e-4);
}

TEST(LevenbergMarquardt_SampleHimmelblauSetBounds)
{
  Unfit::LevenbergMarquardt object;
  object.options.SetMaxIterations(300000);
  object.bounds.SetBounds({-10.0, -10.0}, {2.9, 1.9});
  object.options.SetCostTolerance(1e-24);
  Unfit::Examples::Himmelblau cost_func;
  // Initial guess
  // std::vector<double> min_point = {2.0, 1.0};
  std::vector<double> min_point = {-0.3, -2.0};
  // Check the function calculates the correct cost at the initial guess
// std::vector<double> residual = cost_func(min_point);
// double sum {0.0};
// for(unsigned i = 0; i < residual.size(); ++i) sum += residual[i]*residual[i];
//  CHECK_CLOSE(512, sum, 1e-4);
  // Minimise
  int rc = object.FindMin(cost_func, min_point);
  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(-3.77931, min_point[0], 1e-3);
  CHECK_CLOSE(-3.283186, min_point[1], 1e-3);
}

TEST(LevenbergMarquardt_FindMinGaussian)
{
  Unfit::LevenbergMarquardt test;
  test.options.SetMaxIterations(3000);
  std::vector<double> x {-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0,
                         5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  std::vector<double> y {0.0657285286165, 0.135335283236, 0.249352208777,
                         0.411112290507,  0.606530659712, 0.800737402916,
                         0.945959468906,  1.0, 0.945959468906, 0.800737402916,
                         0.606530659712, 0.411112290507, 0.249352208777,
                         0.135335283236, 0.0657285286165, 0.0285655007845};
  Unfit::Examples::GaussianEquation gauss_eq(x, y);
  std::vector<double> init_guess {1.5, 2.7, 3.3, 7.0};
  auto rc =   test.FindMin(gauss_eq, init_guess);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1.0, init_guess[0], 1e-4);
  CHECK_CLOSE(2.0, init_guess[1], 1e-4);
  CHECK_CLOSE(3.0, init_guess[2], 1e-4);
  CHECK_CLOSE(0.0, init_guess[3], 1e-4);
}

TEST(LevenbergMarquardt_FindMinExponential)
{
  Unfit::LevenbergMarquardt test;
  test.options.SetMaxIterations(3000);
  std::vector<double> x {-10.0, -9.5, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0,
                         -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5,
                         -1.0, -0.5, 0.0, 0.5};
  std::vector<double> y {-0.0111089965382, -0.0132335500965, -0.0157644164848,
                         -0.0187793014946, -0.0223707718561, -0.0266490973363,
                         -0.0317456363780, -0.0378168692293, -0.0450492023935,
                         -0.0536646919127, -0.0639278612067, -0.0761538227986,
                         -0.0907179532894, -0.108067418634,  -0.128734903587,
                         -0.153354966844,  -0.182683524052,  -0.217621056865,
                         -0.259240260645,  -0.30881897968,   -0.367879441171,
                         -0.438234992464};
  Unfit::Examples::Exponential exponential_func(x, y);
//   initial guess is very sensitive to 2nd and 3rd element
  std::vector<double> init_guess {1.0, 0.2, -0.2};
  auto rc = test.FindMin(exponential_func, init_guess);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, init_guess[0], 1e-5);
  CHECK_CLOSE(0.35, init_guess[1], 1e-5);
  CHECK_CLOSE(-1.0, init_guess[2], 1e-5);
}

TEST(LevenbergMarquardt_SampleHimmelblau)
{
  Unfit::LevenbergMarquardt object;
  Unfit::Examples::Himmelblau cost_func;
  object.options.SetCostTolerance(1e-24);
  // Initial guess
  std::vector<double> min_point = {2.0, 3.0};
  // Check the function calculates the correct cost at the initial guess
  std::vector<double> residual = cost_func(min_point);
  double sum {0.0};
  for (unsigned i = 0; i < residual.size(); ++i) sum += residual[i]*residual[i];
  CHECK_CLOSE(512, sum, 1e-4);
  // Minimise
  auto rc = object.FindMin(cost_func, min_point);
  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(3.0, min_point[0], 1e-3);
  CHECK_CLOSE(2.0, min_point[1], 1e-3);  // cant find yet
}

TEST(LevenbergMarquardt_SampleHimmelblau1)
{
  Unfit::LevenbergMarquardt object;
  Unfit::Examples::Himmelblau cost_func;
  object.options.SetCostTolerance(1e-16);
  // Initial guess
  std::vector<double> min_point = {-2.0, -3.0};
  // Check the function calculates the correct cost at the initial guess
  std::vector<double> residual = cost_func(min_point);
  double sum {0.0};
  for (unsigned i = 0; i < residual.size(); ++i) sum += residual[i]*residual[i];
  CHECK_CLOSE(10000, sum, 1e-4);
  // Minimise
  auto rc = object.FindMin(cost_func, min_point);
  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(-3.77931, min_point[0], 1e-4);
  CHECK_CLOSE(-3.28319, min_point[1], 1e-4);  // cant find yet
}
