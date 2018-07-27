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
#include "LinearCostFunction.hpp"
#include "LinearModel.hpp"
#include "ParabolicModel.hpp"
#include "SimpleParabolicCostFunction.hpp"
#include "Unfit.hpp"
#include "UnitTest++.h"

#include "Parabolic.hpp"
#include "GaussianEquation.hpp"
#include "Himmelblau.hpp"
#include "Exponential.hpp"

// This file contains a number of examples showing different ways of doing
// things with Unfit, whether you are doing an optimization, fitting data to
// a function, or something else. In parallel with this, you should look through
// the models/cost functions that are called so you can see how they are
// implemented.
//
// We have chosen to provide all of these examples in terms of tests, but when
// you are writing your own models (or cost functions) it is up to you whether
// you use UnitTest++ or not. To that end, the keywords that show something
// UnitTest++ is doing within a test are those that say CHECK, CHECK_EQUAL, or
// CHECK_CLOSE. Independent of whether or not you chose to use UnitTest++, we
// strongly encourage you to check the return code of anything that returns one,
// be it the DataFileReader or the optimizers. Only by checking the optimizer
// return code will you know if it converged.

SUITE(UnfitExamples)
{
// Here we have no data, no model, just a simple function to minimise, so
// instead of creating a model, we have to create a cost function directly.
// In this case our function is just f(c) = |c*c|, a simple parabola.
TEST(MinimiseASimpleFunctionNoData)
{
  // Create an object that contains our function
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

// Here again we have no data, no model, just a simple function to minimise, so
// instead of creating a model, we have to create a cost function directly.
// In this case our function is quite tricky to minimise. Look at the definition
// in the header to see what the function looks likes.
TEST(MinimiseATrickyFunctionNoData)
{
  // Create an object that contains our function
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

// When you have a model to fit to some data, we recommend you create the model
// by deriving from our GenericNDModel class, and then using our
// GenericNDCostFunction class to turn that into a cost function. We are in the
// process of developing more useful tools for models...
TEST(StraightLine_TwoParameterModelWithData)
{
  // Set the experimental data
  std::vector<std::vector<double>> x {{-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0,
      2.0, 3.0, 4.0, 5.0}};
  std::vector<double> y {-6.88, -4.71, -2.02, -0.98, 1.53, 3.80, 5.49, 7.75,
      9.90, 11.97, 13.98};

  // Make ourselves an object of the model type, in this case a straight line
  Unfit::Examples::LinearModel linear;
  // We want to fit the data by optimizing our parameters, so we need a function
  // to calculate the cost of our model given this data set (x, y)
  Unfit::GenericNDCostFunction linear_cost(linear, x, y);
  // Choose an initial guess for our parameters (c0, c1)
  std::vector<double> c {1.0, 1.0};

  // Try the Nelder-Mead simplex optimisation algorithm
  Unfit::NelderMead nm_opt;
  // Pass in the cost function and our initial guess for c
  auto rc = nm_opt.FindMin(linear_cost, c);
  // Check the optimiser came back with a converged solution
  CHECK_EQUAL(0, rc);
}

// For completeness, here is how you could do the same as the above example,
// but embedding the model within the cost function. As mentioned above, we
// prefer creating models, but if you cannot fit your problem into our model
// framework this gives you a very flexible alternative.
TEST(StraightLine_TwoParameterModelCostFunctionWithData)
{
  // Set the experimental data
  std::vector<double> x {-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0,
      5.0};
  std::vector<double> y {-6.88, -4.71, -2.02, -0.98, 1.53, 3.80, 5.49, 7.75,
      9.90, 11.97, 13.98};

  // Make ourselves a cost function, in this case a straight line
  Unfit::Examples::LinearCostFunction linear_cost(x, y);
  // Choose an initial guess for our parameters (c0, c1)
  std::vector<double> c {1.0, 1.0};

  // Try the Nelder-Mead simplex optimisation algorithm
  Unfit::NelderMead nm_opt;
  // Pass in the cost function and our initial guess for c
  auto rc = nm_opt.FindMin(linear_cost, c);
  // Check the optimiser came back with a converged solution
  CHECK_EQUAL(0, rc);
}

// Back to using our preferred model framework, here is how to do the same
// thing again, but this time reading the data in from a file using the
// DataFileReader class that comes with Unfit.
TEST(StraightLine_ReadDataFromFile)
{
  // Read in the experimental data; first create the DataFileReader
  Unfit::DataFileReader<double> dfr;
  // Get the reader to open the file and import the data. Check the return code
  auto rc = dfr.ReadFile("examples/data/straight_line_data.txt");
  CHECK_EQUAL(0u, rc);
  // Put the data into x and y. We must set the number of dimensions for x.
  std::vector<std::vector<double>> x(1);
  rc = dfr.RetrieveColumn(0, x[0]);
  CHECK_EQUAL(0u, rc);
  std::vector<double> y;
  rc = dfr.RetrieveColumn(1, y);
  CHECK_EQUAL(0u, rc);
  // A quick check x and y have the same number of data points
  CHECK_EQUAL(x[0].size(), y.size());

  // Make ourselves an object of the model type, in this case a straight line
  Unfit::Examples::LinearModel linear;
  // We want to fit the data by optimizing our parameters, so we need a function
  // to calculate the cost of our model given this data set (x, y)
  Unfit::GenericNDCostFunction linear_cost(linear, x, y);
  // Choose an initial guess for our parameters (c0, c1)
  std::vector<double> c {1.0, 1.0};

  // Try the Nelder-Mead simplex optimisation algorithm
  Unfit::NelderMead nm_opt;
  // Pass in the cost function and our initial guess for c
  auto rc2 = nm_opt.FindMin(linear_cost, c);
  // Check the optimiser came back with a converged solution
  CHECK_EQUAL(0, rc2);
}

// ---------------------------


TEST(TwoDimensionTwoParameterModel)
{
  // Read in the experimental data
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/parabolic_data.txt")); //!!!!!
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
