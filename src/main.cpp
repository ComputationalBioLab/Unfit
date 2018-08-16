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
#include <string>
#include <vector>
#include "GenericModelCostFunction.hpp"
#include "ParabolicModel.hpp"
#include "Unfit.hpp"

/**
 * This file provides a simple example of how to use Unfit with a simple
 * cost function. After going through it, take a look at how to write your
 * own cost function so you can perform your own optimizations/fits.
 *
 * In this problem, we have some experimental data, y, at observations, x.
 * We want to fit a parabolic equation to the data of the form:
 *
 *   f(x) = c0*x*x + c1*x + c2
 *
 * such that f(x) best approximates the experimental data. Here we have three
 * parameters to find, c0, c1, and c2. How do we choose values for these
 * parameters? That is where Unfit comes in, for this and many other problems.
 *
 * In terms of the code, the above parabolic equation is stored as a model in
 * the header file ParabolicModel.hpp. Given x and c, it will return f(x). To
 * perform an optimisation, we need to calculate the distance between f(x) and
 * our data, y, which we call the cost, and hence the function to calculate
 * the cost we call a cost function. You could write that yourself, but we have
 * a nice GenericNDCostFunction that is able to convert most models to cost
 * functions without any work on your part.
 *
 * Under the hood, the cost function calculates a residual vector. Consider one
 * data point, we will call it x_i. If we evaluate the model we will get f(x_i),
 * and then we can compare it with the data point y_i. The residual for this
 * data point is calculated as:
 *
 *   r_i = y_i - f(x_i)
 *
 * Do that for each data point and you have the residual vector that our
 * optimisers use to find your parameter values. Here we first solve the problem
 * using the Nelder-Mead Simplex optimisation method. If you are interested, we
 * then go on to solve the same problem using each of the five other algorithms
 * implemented in Unfit: Levenberg-Marquardt, Genetic Algorithm,Differential
 * Evolution, Particle Swarm and Simulated Annealing.
 *
 * Both Nelder-Mead and Levenberg-Marquardt are considered local minimizers and
 * can quickly converge if given a nice cost function and/or a good initial
 * guess. The remaining algorithms are global minimizers and thus, they can be
 * slower, plus they require bounds on your parameters. For most problems
 * Differential Evolution is our choice, followed by either Nelder-Mead or
 * Levenberg-Marquardt. If your problem is a tricky one, have a look at our
 * multi-start, multi-level optimisation in TestExamplesGenerateInitialGuess.cpp
 *
 * FYI, the solution to this problem is approximately
 *   c0 = -21.15
 *   c1 = 61.04
 *   c2 = -9.61
 */
int main()
{
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
  nm_opt.FindMin(parabola_cost, c);
  // ...and voila', your optimized parameter set is returned in c
  std::cout << "NM Result: " << c[0] << " " << c[1] << " " << c[2] << std::endl;
  std::cout << "NM Cost (SSE):   " << nm_opt.GetCost() << std::endl;
  std::cout << std::endl;

  // If that is enough for you - great! Go and start constructing your own
  // models. If not, below we solve the same problem again with all of the
  // different optimisation methods implemented in Unfit.

  // Now try the Levenberg-Marquardt optimisation algorithm
  Unfit::LevenbergMarquardt lm_opt;
  // Use the same initial guess
  c = {1.0, 1.0, 1.0};
  // Pass in the cost function and our initial guess for c
  lm_opt.FindMin(parabola_cost, c);
  // ...and voila' again, your optimized parameter set
  std::cout << "LM Result: " << c[0] << " " << c[1] << " " << c[2] << std::endl;
  std::cout << "LM Cost (SSE):   " << lm_opt.GetCost() << std::endl;
  std::cout << std::endl;

  // Now try a Genetic Algorithm approach
  Unfit::GeneticAlgorithm ga_opt;
  // Here the initial guess just needs to have the correct number of entries
  c = {1.0, 1.0, 1.0};
  // For the Genetic Algorithm, we need to set some sensible bounds
  ga_opt.bounds.SetBounds(0, -100.0, 100.0);  // Bounds on c0
  ga_opt.bounds.SetBounds(1, -100.0, 100.0);  // Bounds on c1
  ga_opt.bounds.SetBounds(2, -100.0, 100.0);  // Bounds on c2
  // Pass in the cost function and our initial guess for c
  ga_opt.FindMin(parabola_cost, c);
  // ...and voila' again, your optimized parameter set
  std::cout << "GA Result: " << c[0] << " " << c[1] << " " << c[2] << std::endl;
  std::cout << "GA Cost (SSE):   " << ga_opt.GetCost() << std::endl;
  std::cout << std::endl;

  // Now try the Differential Evolution algorithm
  Unfit::DifferentialEvolution de_opt;
  // Here the initial guess just needs to have the correct number of entries
  c = {1.0, 1.0, 1.0};
  // For Differential Evolution, we need to set some sensible bounds
  de_opt.bounds.SetBounds(0, -100.0, 100.0);  // Bounds on c0
  de_opt.bounds.SetBounds(1, -100.0, 100.0);  // Bounds on c1
  de_opt.bounds.SetBounds(2, -100.0, 100.0);  // Bounds on c2
  // Pass in the cost function and our initial guess for c
  de_opt.FindMin(parabola_cost, c);
  // ...and voila' again, your optimized parameter set
  std::cout << "DE Result: " << c[0] << " " << c[1] << " " << c[2] << std::endl;
  std::cout << "DE Cost (SSE):   " << de_opt.GetCost() << std::endl;
  std::cout << std::endl;

  // Now try the Particle Swarm algorithm
  Unfit::ParticleSwarm ps_opt;
  // Here the initial guess just needs to have the correct number of entries
  c = {1.0, 1.0, 1.0};
  // For Particle Swarm, we need to set some sensible bounds
  ps_opt.bounds.SetBounds(0, -100.0, 0.0);  // Bounds on c0
  ps_opt.bounds.SetBounds(1, 0.0, 100.0);   // Bounds on c1
  ps_opt.bounds.SetBounds(2, -10.0, 10.0);  // Bounds on c2
  // Pass in the cost function and our initial guess for c
  ps_opt.FindMin(parabola_cost, c);
  // ...and voila' again, your optimized parameter set
  std::cout << "PS Result: " << c[0] << " " << c[1] << " " << c[2] << std::endl;
  std::cout << "PS Cost (SSE):   " << ps_opt.GetCost() << std::endl;
  std::cout << std::endl;

  // Now try the Simulated Annealing algorithm
  Unfit::SimulatedAnnealing sa_opt;
  // Here the initial guess just needs to have the correct number of entries
  std::vector<double> init_guess_sa {1.0, 1.0, 1.0};
  // For Simulated Annealing, we need to set some sensible bounds
  sa_opt.bounds.SetBounds(0, -50.0, 0.0);   // Bounds on c0
  sa_opt.bounds.SetBounds(1, 0.0, 100.0);   // Bounds on c1
  sa_opt.bounds.SetBounds(2, -10.0, 10.0);  // Bounds on c2
  // Pass in the cost function and our initial guess for c
  sa_opt.FindMin(parabola_cost, c);
  // ...and voila' again, your optimized parameter set
  std::cout << "SA Result: " << c[0] << " " << c[1] << " " << c[2] << std::endl;
  std::cout << "SA Cost (SSE):   " << sa_opt.GetCost() << std::endl;
  std::cout << std::endl;

  return 0;
}
