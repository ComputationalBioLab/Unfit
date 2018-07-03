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
#ifndef UNFIT_NOTFORRELEASE_GENERATEINITIALGUESS_HPP_
#define UNFIT_NOTFORRELEASE_GENERATEINITIALGUESS_HPP_

#include <algorithm>
#include <future>
#include <random>
#include <vector>

#include <iostream>

namespace  // file scope
{
template <class T>
std::vector<double> RunOneTrial(Unfit::GenericCostFunction &cost_func,
    const std::vector<double> &lb, const std::vector<double> &ub,
    int number_of_unknowns, int trial, int population_per_trial,
    double cost_tol = 1e-4, double geom_tol = 1e-2)
{
  T opt;
  for (auto i = 0; i < number_of_unknowns; ++i) {
    opt.bounds.SetBounds(i, lb[i], ub[i]);
  }
  opt.options.SetPopulationSize(population_per_trial);
  // Reduced tolerances as we don't need an exact solution here
  opt.options.SetCostTolerance(cost_tol);
  opt.options.SetGeometricTolerance(geom_tol);
  opt.options.SetRandomSeed(10*trial);
  // Initial guess
  std::vector<double> coordinates(number_of_unknowns, 0.0);
  if (!opt.GetIsPopulationBased()) {
    // If the algorithm does not generate a population internally, we will have
    // to generate a random guess here or we will always be starting our
    // optimization from the origin
    std::mt19937 random_engine(opt.options.GetRandomSeed());
    for (auto i = 0; i < number_of_unknowns; ++i) {
      auto dist = std::uniform_real_distribution<double>(lb[i], ub[i]);
      coordinates[i] = dist(random_engine);
    }
  }
  // Minimise
  opt.FindMin(cost_func, coordinates);
  // Append the cost to the coordinates as needed by SetPopulation
  auto residuals = cost_func(coordinates);
  auto cost = std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0);
  coordinates.push_back(cost);
  return coordinates;
}
}  // namespace

namespace Unfit
{
template <class T>
std::vector<std::vector<double>> GenerateInitialPopulation(
    GenericCostFunction &cost_func, const std::vector<double> &lb,
    const std::vector<double> &ub, int number_of_unknowns, int number_of_trials,
    int population_per_trial, double cost_tol = 1e-4,
    double geom_tol = 1e-2)
{
  std::vector<std::vector<double>> population(number_of_trials);
  // Serial implementation
//  for (auto trial = 0; trial < number_of_trials; ++trial) {
//    population[trial] = RunOneTrial<T>(cost_func, lb, ub, number_of_unknowns,
//        trial, population_per_trial, cost_tol, geom_tol);
//  }
  // Threaded (parallel) implementation
  std::vector<std::future<std::vector<double>>>population_futures;
  for (auto trial = 0; trial < number_of_trials; ++trial) {
    population_futures.push_back(std::async(std::launch::async,
        &RunOneTrial<T>, std::ref(cost_func), lb, ub, number_of_unknowns, trial,
        population_per_trial, cost_tol, geom_tol));
  }
  for (auto trial = 0; trial < number_of_trials; ++trial) {
    population[trial] = population_futures[trial].get();
  }
  return population;
}

template <class T>
std::vector<std::vector<double>> GenerateInitialPopulation(
    GenericCostFunction &cost_func, const std::vector<double> &lb,
    const std::vector<double> &ub, int number_of_unknowns)
{
  const int number_of_trials {16};
  const int population_per_trial {16};
  const double cost_tol {1e-4};
  const double geom_tol {1e-2};

  return GenerateInitialPopulation<T>(cost_func, lb, ub, number_of_unknowns,
    number_of_trials, population_per_trial, cost_tol, geom_tol);
}

template <class T>
std::vector<double> GenerateInitialGuess(GenericCostFunction &cost_func,
    const std::vector<double> &lb, const std::vector<double> &ub,
    int number_of_unknowns, int number_of_trials, int population_per_trial,
    double cost_tol = 1e-4, double geom_tol = 1e-2)
{
  std::vector<std::vector<double>> population = GenerateInitialPopulation<T>(
      cost_func, lb, ub, number_of_unknowns, number_of_trials,
      population_per_trial, cost_tol, geom_tol);
  // Sort the population by cost low -> high
  std::sort(begin(population), end(population),
      [](const std::vector<double> &x, const std::vector<double> &y) {
        return y.back() > x.back();
      });
  // Return the best member, with the cost removed
  population[0].pop_back();
  return population[0];
}

template <class T>
std::vector<double> GenerateInitialGuess(GenericCostFunction &cost_func,
    const std::vector<double> &lb, const std::vector<double> &ub,
    int number_of_unknowns)
{
  const int number_of_trials {16};
  const int population_per_trial {10*number_of_unknowns};
  const double cost_tol {1e-4};
  const double geom_tol {1e-2};

  return GenerateInitialGuess<T>(cost_func, lb, ub, number_of_unknowns,
      number_of_trials, population_per_trial, cost_tol, geom_tol);
}
}  // namespace Unfit

#endif
