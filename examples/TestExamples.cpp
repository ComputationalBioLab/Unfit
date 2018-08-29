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
#include <iostream>
#include <vector>
#include "DataFileReader.hpp"
#include "GenerateInitialGuess.hpp"
#include "GenericModelCostFunction.hpp"
#include "HelicalValleyCostFunction.hpp"
#include "HodgkinHuxleyBetaNModel.hpp"
#include "LinearCostFunction.hpp"
#include "LinearModel.hpp"
#include "MyFirstModel.hpp"
#include "NonstationaryMarkovModel.hpp"
#include "Ode2DModel.hpp"
#include "Ode3DModel.hpp"
#include "ParabolicModel.hpp"
#include "SimpleParabolicCostFunction.hpp"
#include "Unfit.hpp"
#include "UnitTest++.h"

// This file contains a number of examples showing different ways of doing
// things with Unfit, whether you are doing an optimization, fitting data to
// a function, or something else. In parallel with this, you should look through
// the models/cost functions that are called so you can see how they are
// implemented. These tests are not exhaustive. To see what else you can do you
// should look at the code and documentation.
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
// First, here is a test that matches the example given in the Unfit tutorial on
// how to create your own model and fit data. The corresponding model code can
// be found in MyFirstModel.hpp in the examples directory.
TEST(TutorialExample)
{
  std::vector<std::vector<double>> t_data {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}};
  std::vector<double> y_data {1.0, 0.6065, 0.3679, 0.2231, 0.1353, 0.0821,
                              0.0498, 0.0302, 0.0183, 0.0111, 0.0067};
  MyFirstModel model;
  Unfit::GenericModelCostFunction model_cost(model, t_data, y_data);
  std::vector<double> c_guess {1.0, 1.0};
  Unfit::NelderMead nm_opt;
  auto rc = nm_opt.FindMin(model_cost, c_guess);
  if (rc == 0) {
    std::cout << "Answer: " << c_guess[0] << " " << c_guess[1] << std::endl;
  } else {
    std::cout << "Oh dear. RC = " << rc << std::endl;
  }
  // An extra check that will cause the test to fail if the optimisation fails
  CHECK_EQUAL(0, rc);
}

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
// by deriving from our GenericModel class, and then using our
// GenericModelCostFunction class to turn that into a cost function. We are in
// the process of developing more useful tools for models...
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
  Unfit::GenericModelCostFunction linear_cost(linear, x, y);
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
  Unfit::GenericModelCostFunction linear_cost(linear, x, y);
  // Choose an initial guess for our parameters (c0, c1)
  std::vector<double> c {1.0, 1.0};

  // Try the Nelder-Mead simplex optimisation algorithm
  Unfit::NelderMead nm_opt;
  // Pass in the cost function and our initial guess for c
  auto rc2 = nm_opt.FindMin(linear_cost, c);
  // Check the optimiser came back with a converged solution
  CHECK_EQUAL(0, rc2);
}

// This example comes from data digitised from Hodgin & Huxley's 1952 model of
// a squid axon. Here we read in the data and fit a two parameter exponential
// model to it.
TEST(TwoDimensionTwoParameterModelMoreData)
{
  // Read in the experimental data
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/hh_beta_n_data.txt"));
  std::vector<std::vector<double>> x(1);
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, x[0]));
  std::vector<double> y;
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, y));
  CHECK_EQUAL(x[0].size(), y.size());

  // Make ourselves an object of the desired model type
  Unfit::Examples::HodgkinHuxleyBetaNModel hh_beta_n;
  // We want to fit the data by optimizing our parameters, so we need a function
  // to calculate the cost of our model given this data set (x, y)
  Unfit::GenericModelCostFunction hh_beta_n_cost(hh_beta_n, x, y);
  // Choose an initial guess for our parameters (c0, c1)
  std::vector<double> c {1.0, 1.0};

  // Try the Nelder-Mead simplex optimisation algorithm
  Unfit::NelderMead nm_opt;
  // Pass in the cost function and our initial guess for c
  auto rc = nm_opt.FindMin(hh_beta_n_cost, c);
  // Check the optimiser came back with a converged solution
  CHECK_EQUAL(0, rc);
}

// This example comes from data digitised from a cardiac ion channel. Here we
// read in the data and fit a three parameter logistic sigmoid function to it.
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
  Unfit::GenericModelCostFunction parabola_cost(parabola, x, y);
  // Choose an initial guess for our parameters (c0, c1, c2)
  std::vector<double> c {1.0, 1.0, 1.0};

  // Try the Nelder-Mead simplex optimisation algorithm
  Unfit::NelderMead nm_opt;
  // Pass in the cost function and our initial guess for c
  auto rc = nm_opt.FindMin(parabola_cost, c);
  // Check the optimiser came back with a converged solution
  CHECK_EQUAL(0, rc);
}

// This model is a bit more complex. Here we have a 3D model, eta = f(t, vm),
// or if you would like it in a more generic way, y = f(x0, x1). Have a look at
// the header file to see the function - it is a type of probability
// distribution. You can see that apart from a few more lines of code organising
// the data, nothing else changes as the model becomes more difficult.
TEST(ThreeDimensionFourParameterModel)
{
  // Read in the experimental data
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/nonstationary_markov_data.txt"));
  // We need to initialise our x with the number of x dimensions we have
  std::vector<std::vector<double>> x(2);
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, x[0]));  // time = x[0]
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, x[1]));  // vm = x[1]
  std::vector<double> eta;
  CHECK_EQUAL(0u, dfr.RetrieveColumn(2, eta));
  CHECK_EQUAL(x[0].size(), eta.size());
  CHECK_EQUAL(x[1].size(), eta.size());

  // Create the model then the cost function with our data
  Unfit::Examples::NonstationaryMarkovModel nsm;
  Unfit::GenericModelCostFunction nsm_cost(nsm, x, eta);
  // The initial guess for our model parameters
  std::vector<double> c {10.0, 10.0, 1.0, 1.0};

  // Find the best fit parameters
  Unfit::NelderMead nm;
  auto rc = nm.FindMin(nsm_cost, c);
  CHECK_EQUAL(0, rc);
}

// This model shows an example of how you could approach a problem that is
// in the form of an ordinary differential equation (ODE). Here we have a
// problem of the form dy/dt = f(y), along with some data we want to try to
// match. Here we choose the forward Euler method to integrate the ODE and
// thus obtain model values we can compare to the experimental data. Note that
// as this is an initial value problem, we pass in the initial value for y
// along with the time step when we construct the model.
TEST(OrdinaryDifferentialEquationModel)
{
  // Read in the experimental data
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/ode_data.txt"));
  // We need to initialise our t with the number of t dimensions we have
  std::vector<std::vector<double>> t(1);
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, t[0]));
  std::vector<double> y;
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, y));

  // Create the model then the cost function with our data
  auto ic = 0.0;  // initial condition
  // Here we will get the time step from the time points in the data file
  // We could also do this within the model if we wanted to
  auto dt = t[0][1] - t[0][0];
  Unfit::Examples::Ode2DModel ode2d(dt, ic);
  Unfit::GenericModelCostFunction ode2d_cost(ode2d, t, y);
  // The initial guess for our model parameters
  std::vector<double> c = {1.0, 1.0};

  // Find the best fit parameters
  Unfit::NelderMead nm;
  auto rc = nm.FindMin(ode2d_cost, c);
  CHECK_EQUAL(0, rc);
}

// This model shows is a variation of the above ODE example whereby we change
// two things. First, we assume that we do not know what the initial value
// should be for our initial value problem. We therefore add it as an additional
// variable to our optimisation problem. The other change is that we no longer
// pass in a time step. Instead, this version looks at the times we read in from
// the data file and uses those to calculated the time step on the fly.
TEST(OrdinaryDifferentialEquationModelUnknownInitialCondition)
{
  // Read in the experimental data
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/ode_data.txt"));
  // We need to initialise our t with the number of t dimensions we have
  std::vector<std::vector<double>> t(1);
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, t[0]));
  std::vector<double> y;
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, y));

  // Create the model then the cost function with our data
  Unfit::Examples::Ode3DModel ode3d;
  Unfit::GenericModelCostFunction ode3d_cost(ode3d, t, y);
  // The initial guess for our model parameters
  std::vector<double> c = {1.0, 1.0, 1.0};

  // Find the best fit parameters
  Unfit::NelderMead nm;
  auto rc = nm.FindMin(ode3d_cost, c);
  CHECK_EQUAL(0, rc);
}

// This example demonstrates how to get the information you need out of the
// Unfit optimisers.
TEST(GettinOutputFromUnfit)
{
  // See the StraightLine_ReadDataFromFile example for details about this part
  Unfit::DataFileReader<double> dfr;
  dfr.ReadFile("examples/data/straight_line_data.txt");
  std::vector<std::vector<double>> x(1);
  dfr.RetrieveColumn(0, x[0]);
  std::vector<double> y;
  dfr.RetrieveColumn(1, y);
  Unfit::Examples::LinearModel linear;
  Unfit::GenericModelCostFunction linear_cost(linear, x, y);
  std::vector<double> c {1.0, 1.0};
  Unfit::NelderMead nm_opt;

  // Q. Can I get output as I go to see if the optimiser is working?
  // A. Yes, all optimisers have a SetOutputLevel method as an option. You can
  //    change the output level before you run the optimisation. Unfit does not
  //    store intermediate states, so changing if afterwards will not do
  //    anything.
  nm_opt.options.SetOutputLevel(1);  // Iteration counter
  nm_opt.options.SetOutputLevel(2);  // Iteration by iteration information
  nm_opt.options.SetOutputLevel(3);  // Level 2 plus a summary at the end
  nm_opt.options.SetOutputLevel(0);  // No output (default)

  // Q. How to I know that the optimisation worked?
  // A. Check the return code from FindMin. If the return code is zero the
  //    optimisation process did not encounter any problems. This does not rule
  //    out you hitting a bound (if you have them) or them finding a local
  //    minimum, so you may need to used your common sense.
  auto rc = nm_opt.FindMin(linear_cost, c);
  CHECK_EQUAL(0, rc);

  // Q. How to you know what the answer is?
  // A. They are returned in the same vector you used for your initial guess.
  //    In this case they are returned in c[0] and c[1], which you can access
  //    directly.
  CHECK_CLOSE(2.07455, c[0], 1e-2);
  CHECK_CLOSE(3.62091, c[1], 1e-2);

  // Q. How to you know what the final cost was?
  // A. All optimisers have a GetCost method.
  CHECK_CLOSE(0.64202, nm_opt.GetCost(), 1e-2);
}

// There are several common things we might like to do when running an
// optimisation problem, such as controlling the tolerance we need to obtain to
// say we have an optimal solution. In addition, Unfit counts the number of
// iterations and function evaluations as the optimisation proceeds, and
// terminates if there are too many. Increasing these limits can also be useful.
// In addition, the optimisers  have several control parameters that you can
// play with to control their behaviour. You can see them all in Options.hpp but
// here we given an example of how to the use this functionality.
TEST(WorkingWithOptimiserOptions)
{
  // See the StraightLine_ReadDataFromFile example for details about this part
  Unfit::DataFileReader<double> dfr;
  dfr.ReadFile("examples/data/straight_line_data.txt");
  std::vector<std::vector<double>> x(1);
  dfr.RetrieveColumn(0, x[0]);
  std::vector<double> y;
  dfr.RetrieveColumn(1, y);
  Unfit::Examples::LinearModel linear;
  Unfit::GenericModelCostFunction linear_cost(linear, x, y);
  std::vector<double> c {1.0, 1.0};

  // Here we will choose to use the Differential Evolution method
  Unfit::DifferentialEvolution de_opt;
  // First we set up the options that we want to use.
  de_opt.options.SetMaxIterations(10000);
  de_opt.options.SetMaxFunctionEvaluations(10000);
  de_opt.options.SetCostTolerance(1.0e-6);
  de_opt.options.SetGeometricTolerance(1.0e-6);
  de_opt.options.SetPopulationSize(15);

  // Now that all of our options are in place, we can run the optimisation
  auto rc = de_opt.FindMin(linear_cost, c);
  CHECK_EQUAL(0, rc);
}

// Sometimes we have a model and we know, for example, that the parameters must
// be positive. It could be a mechanics problem when we know the material
// constants must be positive or something else. Unfit provides what are known
// as Box constraints - you can set an upper and/or lower bound for each
// parameter. Here we choose a simple parabola with one parameter that has a
// minimum at x=0, but we impose a lower bound that x must be >2. Note that
// Unfit will not tell you if you hit a bound, you should check that yourself
// upon exit.
TEST(SimpleExampleWithBounds)
{
  // Create an object that contains our function
  Unfit::Examples::SimpleParabolicCostFunction sp_cost;
  // Choose an initial guess for our parameter (c0)
  std::vector<double> c {3.0};

  // Try the Levenberg-Marquardt optimisation algorithm
  Unfit::LevenbergMarquardt lm_opt;
  // Add a lower bound to our c0 parameter so it must be > 2
  lm_opt.bounds.SetLowerBound(0, 2.0);
  // Pass in the cost function and our initial guess for c
  auto rc = lm_opt.FindMin(sp_cost, c);
  // Check the optimiser came back with a converged solution
  CHECK_EQUAL(0, rc);
  // Check the solution (we know the exact solution in this case)
  CHECK_CLOSE(2.0, c[0], 1e-3);
  // (Note that Nelder-Mead will complain about a degenerate simplex for this
  // example (rc=4), as the simplex collapses on itself at the bound, but it
  // still gets the correct answer)
}

// In some cases we may want to set bounds on all of our parameters. We can do
// this by creating two vectors, one containing the lower bounds and one
// containing the upper bounds. (Note: the terms lower_bound and upper_bound are
// reserved in C++, so call your vectors something different). Here we have an
// optimisation with two parameters, and say we want to ensure both parameters
// are positive. We can make a lower bound vector containing two zeros, and an
// upper bound vector containing suitably large values. We use bound.SetBounds
// to impose these as bounds on the optimisation.
TEST(ExampleWithUpperAndLowerBounds)
{
  std::vector<std::vector<double>> x {{-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0,
      2.0, 3.0, 4.0, 5.0}};
  std::vector<double> y {-6.88, -4.71, -2.02, -0.98, 1.53, 3.80, 5.49, 7.75,
      9.90, 11.97, 13.98};
  Unfit::Examples::LinearModel linear;
  Unfit::GenericModelCostFunction linear_cost(linear, x, y);
  std::vector<double> c {1.0, 1.0};

  // Use the Nelder-Mead simplex optimisation algorithm
  Unfit::NelderMead nm_opt;
  // The bounds on all of the parameters can be set via a pair of vectors.
  // The length of these vectors must match the length of the c vector.
  std::vector<double> lower_bnd {0.0, 0.0};
  std::vector<double> upper_bnd {1.0e6, 1.0e6};
  nm_opt.bounds.SetBounds(lower_bnd, upper_bnd);
  // Pass in the cost function and our initial guess for c
  auto rc = nm_opt.FindMin(linear_cost, c);
  // Check the optimiser came back with a converged solution
  CHECK_EQUAL(0, rc);
}

// For the global optimisers, Differential Evolution, Genetic Algorithm,
// Particle Swarm, and Simulated Annealing, you should set bounds on all of your
// parameters. DE, GA & PS are population-based, so use these bounds to generate
// an initial population. SA uses the bounds to generate new trial points. For
// all of these, the default bounds are very very large (+ve and -ve) so unless
// you set bounds yourself, based on the scale of your problem, the algorithms
// will likely take an age to converge.
//
// Note also that due to their algorithms, DE & PS can generate points outside
// the initial bounds that are given during their iterations. By default we
// allow them to do this. If you want to strictly keep these methods within the
// bounds you set, then set the options.SetUseHardBounds to true.
//
// (The optimisers that do not require bounds (but can still handle them) are
// Nelder Mead and Levenberg Marquardt).
TEST(GlobalOptimizersRequireBounds)
{
  std::vector<std::vector<double>> x {{-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0,
      2.0, 3.0, 4.0, 5.0}};
  std::vector<double> y {-6.88, -4.71, -2.02, -0.98, 1.53, 3.80, 5.49, 7.75,
      9.90, 11.97, 13.98};
  Unfit::Examples::LinearModel linear;
  Unfit::GenericModelCostFunction linear_cost(linear, x, y);
  std::vector<double> c {1.0, 1.0};

  // Use the Differential Evolution optimisation algorithm
  Unfit::DifferentialEvolution de_opt;
  // Set upper and lower bounds on each parameter
  std::vector<double> lower_bnd {-100.0, -100.0};
  std::vector<double> upper_bnd {100.0, 100.0};
  de_opt.bounds.SetBounds(lower_bnd, upper_bnd);
  // Restrict DE such that it enforces the bounds for all time
  de_opt.options.SetUseHardBounds(true);
  // Pass in the cost function and our initial guess for c
  auto rc = de_opt.FindMin(linear_cost, c);
  // Check the optimiser came back with a converged solution
  CHECK_EQUAL(0, rc);
}

// By default, the optimisation methods in Unfit (with the exception of Nelder
// Mead and Levenberg Marquardt) will generate a random initial guess or
// population of guesses for the unknown parameters. These will be generated
// within the bounds that have been set. Sometimes you may already have an
// idea of where the minimum is (e.g. from a previous time step or solution) so
// you may want to take advantage of this. There is an option that allows you
// to add the guess you pass in to the initial population, and this example
// shows it in action. Note that the population based methods also allow you
// to set the whole population, and Nelder Mead allows you to set the initial
// simplex. See the GenericOptimizer SetPopulation documentation for details.
TEST(AddInitialGuessToAPopulation)
{
  std::vector<std::vector<double>> x {{-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0,
      2.0, 3.0, 4.0, 5.0}};
  std::vector<double> y {-6.88, -4.71, -2.02, -0.98, 1.53, 3.80, 5.49, 7.75,
      9.90, 11.97, 13.98};
  Unfit::Examples::LinearModel linear;
  Unfit::GenericModelCostFunction linear_cost(linear, x, y);
  std::vector<double> c {1.0, 1.0};

  // Use the Differential Evolution optimisation algorithm
  Unfit::DifferentialEvolution de_opt;
  // Set upper and lower bounds on each parameter
  std::vector<double> lower_bnd {-100.0, -100.0};
  std::vector<double> upper_bnd {100.0, 100.0};
  de_opt.bounds.SetBounds(lower_bnd, upper_bnd);
  // Tell the optimiser to add our initial guess (1.0, 1.0) to the initial
  // population. The remainder of the initial population will be randomly
  // generated.
  de_opt.options.SetAddInitialToPopulation(true);
  // Pass in the cost function and our initial guess for c
  auto rc = de_opt.FindMin(linear_cost, c);
  // Check the optimiser came back with a converged solution
  CHECK_EQUAL(0, rc);
}

// If you have a very large data set that is taking a long time to optimise, it
// may be useful to downsample the data set and run an initial optimisation on
// a smaller data set to get an estimate of the parameter values (which should
// be much faster, depending on how aggressive the downsampling). You can then
// use those parameters as an initial guess with the larger data set to get the
// proper solution. This example shows how to do just this, and could also be
// applicable if you have multiple data sets but expect the parameters to be
// similar for each set.
TEST(MultipleDataSetsAndOptimisations)
{
  // A sample low resolution data set
  std::vector<std::vector<double>> x_low {{-5.0, -1.0, 1.0, 5.0}};
  std::vector<double> y_low {-6.88, 1.53, 5.49, 13.98};
  Unfit::Examples::LinearModel linear;
  Unfit::GenericModelCostFunction linear_cost(linear, x_low, y_low);
  std::vector<double> c {1.0, 1.0};

  // Use the Differential Evolution optimisation algorithm
  Unfit::DifferentialEvolution de_opt;
  // Set upper and lower bounds on each parameter
  std::vector<double> lower_bnd {-100.0, -100.0};
  std::vector<double> upper_bnd {100.0, 100.0};
  de_opt.bounds.SetBounds(lower_bnd, upper_bnd);
  // Pass in the cost function and our initial guess for c
  auto rc = de_opt.FindMin(linear_cost, c);
  // Check the optimiser came back with a converged solution
  CHECK_EQUAL(0, rc);

  // Now we have an approximate solution based on our low resolution data set,
  // (stored in c), we can fit our high resolution data set
  std::vector<std::vector<double>> x_high {{-5.0, -4.0, -3.0, -2.0, -1.0, 0.0,
      1.0, 2.0, 3.0, 4.0, 5.0}};
  std::vector<double> y_high {-6.88, -4.71, -2.02, -0.98, 1.53, 3.80, 5.49,
      7.75, 9.90, 11.97, 13.98};
  // Change the data in our cost function to the high resolution version
  linear_cost.SetData(x_high, y_high);
  // Reset the optimiser (we don't need to create another)
  de_opt.Reset();
  // Re-add the bounds
  de_opt.bounds.SetBounds(lower_bnd, upper_bnd);
  // The c vector already contains the best fit from our first optimisation, so
  // we must add it to the population to make use of this information.
  de_opt.options.SetAddInitialToPopulation(true);
  // Pass in the cost function and our current estimate for c
  rc = de_opt.FindMin(linear_cost, c);
  // Check the optimiser came back with a converged solution
  CHECK_EQUAL(0, rc);
}

// For certain problems it can be advantageous to use a global optimiser to get
// close to the optimal parameter values before switching to a local optimiser
// to refine this solution and return the optimised parameter set. This is one
// form of what is known as multi-level optimisation and is demonstrated in the
// example below.
TEST(MultiLevelOptimisationGlobalThenLocal)
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
  Unfit::GenericModelCostFunction parabola_cost(parabola, x, y);
  // Choose an initial guess for our parameters (c0, c1, c2)
  std::vector<double> c {1.0, 1.0, 1.0};

  // First use a global optimiser to get close to the solution. In this case
  // we have chosen to use Particle Swarm.
  Unfit::ParticleSwarm ps_opt;
  // Set upper and lower bounds on each parameter
  std::vector<double> lower_bnd {-100.0, -100.0, -100.0};
  std::vector<double> upper_bnd {100.0, 100.0, 100.0};
  ps_opt.bounds.SetBounds(lower_bnd, upper_bnd);
  // As we are only looking to get in the vicinity of the correct answer, we
  // can choose to use loose tolerances (optional)
  ps_opt.options.SetCostTolerance(1.0e-2);
  ps_opt.options.SetGeometricTolerance(1.0e-2);
  // Perform the optimisation
  auto rc = ps_opt.FindMin(parabola_cost, c);
  // Check the optimiser came back with a converged solution
  CHECK_EQUAL(0, rc);

  // Now the vector c contains the resulting parameters from the Particle Swarm
  // optimisation, we can use a local optimiser (Nelder Mead, Levenberg
  // Marquardt) to refine our solution.
  Unfit::LevenbergMarquardt lm_opt;
  // Pass in the cost function and our current guess for c
  rc = lm_opt.FindMin(parabola_cost, c);
  // Check the optimiser came back with a converged solution
  CHECK_EQUAL(0, rc);
}

// Certain problems, particularly those which have multiple local minima, can
// benefit from a more complete coverage of the parameter space as can be
// achieved with a multi-start optimisation approach. Here we run multiple
// instances of the Differential Evolution method from different starting points
// and then choose the best solution across all of the instances as an initial
// guess for a local optimiser (Nelder Mead). The multiple instances of
// Differential Evolution run in parallel, if parallel hardware is available,
// which can significantly speed up the optimisation process.
TEST(MultiLevelMultiStartInitialGuess)
{
  // Set the experimental data
  std::vector<std::vector<double>> x {{0.498531, 0.622145, 0.746551, 0.899687,
      0.995019, 1.24803, 1.49695, 1.7464, 1.86737, 1.92478, 2.07206, 2.12789,
      2.23212}};
  std::vector<double> y {17.0676, 20.0356, 24.0914, 27.5963, 28.9598, 31.9736,
      34.6866, 33.7931, 31.9415, 30.6897, 28.3853, 23.9687, 18.5146};

  // Create the model and cost function
  Unfit::Examples::ParabolicModel parabola;
  Unfit::GenericModelCostFunction parabola_cost(parabola, x, y);

  // Choose an initial guess for our parameters (c0, c1, c2)
  std::vector<double> c {1.0, 1.0, 1.0};
  const int number_of_unknowns =c.size();

  // Set bounds on our parameters
  std::vector<double> lower_bnd {-100.0, -100.0, -100.0};
  std::vector<double> upper_bnd {100.0, 100.0, 100.0};

  // Run multiple instances of the Differential Evolution optimiser, each with
  // a different random seed, to provide wide coverage of the parameter space.
  // The optimiser we want to use in this first level, (which is multi-start) is
  // passed in as a template parameter. There is also an alternate interface
  // which gives a lot more control over what the multi-start does.
  c = Unfit::GenerateInitialGuess<Unfit::DifferentialEvolution>(parabola_cost,
      lower_bnd, upper_bnd, number_of_unknowns);

  // The vector c now contains the best solution across all of the Differential
  // Evolution instances that were run. We can then use this as a starting point
  // for a local optimiser, e.g. Nelder Mead.
  Unfit::NelderMead nm;
  auto rc = nm.FindMin(parabola_cost, c);
  // Check the optimiser came back with a converged solution
  CHECK_EQUAL(0, rc);
}

// Finding a good initial population with multi-level, multi-start - note MultiThreaded versions
TEST(MultiLevelMultiStartInitialPopulation)
{
  // Set the experimental data
  std::vector<std::vector<double>> x {{0.498531, 0.622145, 0.746551, 0.899687,
      0.995019, 1.24803, 1.49695, 1.7464, 1.86737, 1.92478, 2.07206, 2.12789,
      2.23212}};
  std::vector<double> y {17.0676, 20.0356, 24.0914, 27.5963, 28.9598, 31.9736,
      34.6866, 33.7931, 31.9415, 30.6897, 28.3853, 23.9687, 18.5146};

  // Create the model and cost function
  Unfit::Examples::ParabolicModel parabola;
  Unfit::GenericModelCostFunction parabola_cost(parabola, x, y);

  // Choose an initial guess for our parameters (c0, c1, c2)
  std::vector<double> c {1.0, 1.0, 1.0};
  const int number_of_unknowns =c.size();

  // Set bounds on our parameters
  std::vector<double> lower_bnd {-100.0, -100.0, -100.0};
  std::vector<double> upper_bnd {100.0, 100.0, 100.0};

  // We want to use a multi-start approach to generate an initial population
  // for the Differential Evolution algorithm
  Unfit::DifferentialEvolution de_opt;
  // The initial Differential Evolution population is set to be the result of
  // a multi-start Particle Swarm approach. Note that we could use Differential
  // evolution for this too if we wanted to. There is also an alternate
  // interface which gives a lot more control over what the multi-start does.
  // Here the best member from each swarm is added to the initial Differential
  // Evolution population.
  de_opt.SetPopulation(Unfit::GenerateInitialPopulation<Unfit::ParticleSwarm>(
      parabola_cost, lower_bnd, upper_bnd, number_of_unknowns));
  // Run the Differential Evolution optimisation and return the result in c
  auto rc = de_opt.FindMin(parabola_cost, c);
  // Check the optimiser came back with a converged solution
  CHECK_EQUAL(0, rc);
}

// a fast escape from a cost function if something goes wrong ?? Need to be careful, modify the GenericModelCostFunction


// This test is just a copy of what is in the main.cpp file that is run when you
// use the non-test version of Unfit, and as such, should be kept up to date if
// main.cpp changes. It is added as a test here so we make sure we don't break
// non-test Unfit. It is also a good example to see how we can switch easily
// between optimizers for the the same model.
TEST(MainCPPUnfit)
{
  // A tolerance for testing the solutions against
  const double tol = 1.0e-3;

  // First we set/get our data points. You could also read them in; we have
  // a DataFileReader class to help with that.
  std::vector<std::vector<double>> x {{0.498531, 0.622145, 0.746551, 0.899687,
      0.995019, 1.24803, 1.49695, 1.7464, 1.86737, 1.92478, 2.07206, 2.12789,
      2.23212}};
  std::vector<double> y {17.0676, 20.0356, 24.0914, 27.5963, 28.9598, 31.9736,
      34.6866, 33.7931, 31.9415, 30.6897, 28.3853, 23.9687, 18.5146};

  // Make ourselves an object of the desired model type, in this case a parabola
  Unfit::Examples::ParabolicModel parabola;
  // We want to fit the data by optimizing our parameters, so we need a function
  // to calculate the cost of our model given this data set (x, y)
  Unfit::GenericModelCostFunction parabola_cost(parabola, x, y);
  // Choose an initial guess for (c0, c1, c2)
  std::vector<double> c {1.0, 1.0, 1.0};

  // First, let's try the Nelder-Mead simplex optimisation algorithm
  Unfit::NelderMead nm_opt;
  // Pass in the cost function and our initial guess for c
  auto rc = nm_opt.FindMin(parabola_cost, c);
  // ...and voila', your optimized parameter set is returned in c
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(-21.115, c[0], tol);
  CHECK_CLOSE(61.0367, c[1], tol);
  CHECK_CLOSE(-9.61497, c[2], tol);
  CHECK_CLOSE(24.9151, nm_opt.GetCost(), tol);

  // If that is enough for you - great! Go and start constructing your own
  // models. If not, below we solve the same problem again with all of the
  // different optimisation methods implemented in Unfit.

  // Now try the Levenberg-Marquardt optimisation algorithm
  Unfit::LevenbergMarquardt lm_opt;
  // Use the same initial guess
  c = {1.0, 1.0, 1.0};
  // Pass in the cost function and our initial guess for c
  rc = lm_opt.FindMin(parabola_cost, c);
  // ...and voila' again, your optimized parameter set
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(-21.115, c[0], tol);
  CHECK_CLOSE(61.0367, c[1], tol);
  CHECK_CLOSE(-9.61497, c[2], tol);
  CHECK_CLOSE(24.9151, lm_opt.GetCost(), tol);

  // Now try a Genetic Algorithm approach
  Unfit::GeneticAlgorithm ga_opt;
  // Here the initial guess just needs to have the correct number of entries
  c = {1.0, 1.0, 1.0};
  // For the Genetic Algorithm, we need to set some sensible bounds
  ga_opt.bounds.SetBounds(0, -100.0, 100.0);  // Bounds on c0
  ga_opt.bounds.SetBounds(1, -100.0, 100.0);  // Bounds on c1
  ga_opt.bounds.SetBounds(2, -100.0, 100.0);  // Bounds on c2
  // Pass in the cost function and our initial guess for c
  rc = ga_opt.FindMin(parabola_cost, c);
  // ...and voila' again, your optimized parameter set
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(-19.7735, c[0], tol);
  CHECK_CLOSE(57.214, c[1], tol);
  CHECK_CLOSE(-7.3312, c[2], tol);
  CHECK_CLOSE(26.4066, ga_opt.GetCost(), tol);

  // Now try the Differential Evolution algorithm
  Unfit::DifferentialEvolution de_opt;
  // Here the initial guess just needs to have the correct number of entries
  c = {1.0, 1.0, 1.0};
  // For Differential Evolution, we need to set some sensible bounds
  de_opt.bounds.SetBounds(0, -100.0, 100.0);  // Bounds on c0
  de_opt.bounds.SetBounds(1, -100.0, 100.0);  // Bounds on c1
  de_opt.bounds.SetBounds(2, -100.0, 100.0);  // Bounds on c2
  // Pass in the cost function and our initial guess for c
  rc = de_opt.FindMin(parabola_cost, c);
  // ...and voila' again, your optimized parameter set
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(-21.115, c[0], tol);
  CHECK_CLOSE(61.0367, c[1], tol);
  CHECK_CLOSE(-9.61497, c[2], tol);
  CHECK_CLOSE(24.9151, de_opt.GetCost(), tol);

  // Now try the Particle Swarm algorithm
  Unfit::ParticleSwarm ps_opt;
  // Here the initial guess just needs to have the correct number of entries
  c = {1.0, 1.0, 1.0};
  // For Particle Swarm, we need to set some sensible bounds
  ps_opt.bounds.SetBounds(0, -100.0, 0.0);  // Bounds on c0
  ps_opt.bounds.SetBounds(1, 0.0, 100.0);   // Bounds on c1
  ps_opt.bounds.SetBounds(2, -10.0, 10.0);  // Bounds on c2
  // Pass in the cost function and our initial guess for c
  rc = ps_opt.FindMin(parabola_cost, c);
  // ...and voila' again, your optimized parameter set
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(-20.2212, c[0], tol);
  CHECK_CLOSE(58.4493, c[1], tol);
  CHECK_CLOSE(-8.02656, c[2], tol);
  CHECK_CLOSE(25.6234, ps_opt.GetCost(), tol);

  // Now try the Simulated Annealing algorithm
  Unfit::SimulatedAnnealing sa_opt;
  // Here the initial guess just needs to have the correct number of entries
  std::vector<double> init_guess_sa {1.0, 1.0, 1.0};
  // For Simulated Annealing, we need to set some sensible bounds
  sa_opt.bounds.SetBounds(0, -50.0, 0.0);   // Bounds on c0
  sa_opt.bounds.SetBounds(1, 0.0, 100.0);   // Bounds on c1
  sa_opt.bounds.SetBounds(2, -10.0, 10.0);  // Bounds on c2
  // Pass in the cost function and our initial guess for c
  rc = sa_opt.FindMin(parabola_cost, c);
  // ...and voila' again, your optimized parameter set
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(-21.115, c[0], tol);
  CHECK_CLOSE(61.0367, c[1], tol);
  CHECK_CLOSE(-9.61495, c[2], tol);
  CHECK_CLOSE(24.9151, sa_opt.GetCost(), tol);
}
}  // suite UnfitExamples
