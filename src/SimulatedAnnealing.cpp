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
#include <cmath>
#include <random>
#include <vector>
#include "Bounds.hpp"
#include "Options.hpp"
#include "SimulatedAnnealing.hpp"

namespace Unfit
{
// Default constructor: keep Reset in sync with this if you change it
SimulatedAnnealing::SimulatedAnnealing()
  : cost_ {0}
  , dimensions_ {0}
  , previous_best_cost_ {0.0}
  , step_sizes_ {}
  , acceptance_ratios_ {}
  , generator_ {0}
  , uniform_dist_(0.0, 1.0)
{
  is_population_based_ = false;
}

void SimulatedAnnealing::Reset()
{
  ResetGenericOptimizer();
  cost_ = 0;
  dimensions_ = 0;
  previous_best_cost_ = 0.0;
  step_sizes_.clear();
  acceptance_ratios_.clear();
  generator_.seed(0);
  uniform_dist_.reset();
  // Remember not to reset is_population_based_
}

void SimulatedAnnealing::InitialiseParameters()
{
  const auto initial_step_size = 1.0;
  step_sizes_.assign(dimensions_, initial_step_size);
  const auto initial_acceptance_ratio = 1.0;
  acceptance_ratios_.assign(dimensions_, initial_acceptance_ratio);
  // Set the correct number of default bounds (does not overwrite existing ones)
  bounds.SetNumberOfBounds(dimensions_);
  // Seed the random number generator with the specified seed
  generator_.seed(options.GetRandomSeed());
}

void SimulatedAnnealing::GenerateTrialPoint(std::vector<double> &trial_point,
    int i)
{
  // Take a simple step of random length in the selected direction. If we go
  // out of bounds just clamp it within the bounds.
  const auto distance = 2.0 * uniform_dist_(generator_) - 1.0;
  trial_point[i] += distance * step_sizes_[i];
  bounds.ClampWithinBounds(trial_point);

//  // Alternatively we could choose a random point within the bounds.
//  const auto distance = 2.0 * uniform_dist_(generator_) - 1.0;
//  trial_point[i] += distance * step_sizes_[i];
//  if (!bounds.IsWithinBounds(i, trial_point[i])) {
//    const auto lb = bounds.GetLowerBound(i);
//    const auto ub = bounds.GetUpperBound(i);
//    trial_point[i] = lb + uniform_dist_(generator_) * (ub - lb);
//  }

//  // One option is to keep looking for a candidate point until we find one
//  // within our bounds, but if our bounds are tight or our step size is large
//  // then this can stall and increase solution time.
//  //
//  const auto lb = bounds.GetLowerBound(i);
//  const auto ub = bounds.GetUpperBound(i);
//  do {
//    const auto distance = 2.0 * uniform_dist_(generator_) - 1.0;
//    trial_point[i] += distance * step_sizes_[i];
//  } while (trial_point[i] > ub || trial_point[i] < lb);

//  // This is a more complicated method of selecting a trial point, suggested
//  // by Inger (Very Fast Simulated Reanealing), but not really needed at the
//  // moment. Leave here in case we need it.
//  //
//  // The signum function it uses can be easily implemented as:
//  // template <typename T> int signum(T val) {
//  //   return (T(0) < val) - (val < T(0));
//  // }
//  //
//  const auto lb = bounds.GetLowerBound(i);
//  const auto ub = bounds.GetUpperBound(i);
//  do {
//    const auto y = signum(uniform_dist_(generator_) - 0.5) * temperature *
//        (std::pow((1.0 + 1.0 / temperature), std::fabs(2.0 *
//        uniform_dist_(generator_) - 1.0)) - 1.0);
//    trial_point[i] += y * (ub - lb);
//  } while (trial_point[i] > ub || trial_point[i] < lb);
}

void SimulatedAnnealing::UpdateStepSizes() noexcept
{
  for (auto i = 0u; i < dimensions_; ++i) {
    // These are Equation 7.22 from the Belegundu book, page 316
    if (acceptance_ratios_[i] > 0.6) {
      step_sizes_[i] *= (1.0 + 2.0 * ((acceptance_ratios_[i] - 0.6) / 0.4));
    }
    else if (acceptance_ratios_[i] < 0.4) {
      step_sizes_[i] /= (1.0 + 2.0 * ((0.4 - acceptance_ratios_[i]) / 0.4));
    }
    const auto lb = bounds.GetLowerBound(i);
    const auto ub = bounds.GetUpperBound(i);
    if (step_sizes_[i] > (ub - lb)) {
      step_sizes_[i] = ub - lb;
    }
    acceptance_ratios_[i] = 1.0;
  }
}

void SimulatedAnnealing::ResetStepSizes(double step_size) noexcept
{
  for (auto i = 0u; i < dimensions_; ++i) {
    const auto lb = bounds.GetLowerBound(i);
    const auto ub = bounds.GetUpperBound(i);
    step_sizes_[i] = step_size * (ub - lb);
  }
}

int SimulatedAnnealing::FindMin(GenericCostFunction &CostFunction,
    std::vector<double> &coordinates)
{
  // Check the initial guess is a non-zero size
  if (coordinates.empty()) return -1;
  dimensions_ = coordinates.size();
  cost_ = dimensions_;
  // Increase coordinate vector size by 1 as we append the cost to the end
  coordinates.push_back(0.0);
  // Now we have the size of the problem we can set the member variables
  InitialiseParameters();
  // If we have not asked to use the supplied guess, generate a random one
  if (!options.GetAddInitialToPopulation()) {
    do {
      for (auto i = 0u; i < dimensions_; ++i) {
        const auto lb = bounds.GetLowerBound(i);
        const auto ub = bounds.GetUpperBound(i);
        coordinates[i] = lb + uniform_dist_(generator_) * (ub - lb);
      }
    } while (!CalculateCost(CostFunction, coordinates));
  }
  else {
    // Check the supplied coordinates are within the specified bounds
    if (!bounds.IsWithinBounds(coordinates)) return -2;
    // Check if the cost of the start points is finite
    if (!CalculateCost(CostFunction, coordinates)) return -3;
  }
  // Store the best cost for convergence checking
  previous_best_cost_ = coordinates[cost_];

  auto rc = ProcessFindMin(CostFunction, coordinates);

  // Add the best point to the population so we can use the print and output
  // functions from GenericOptimizer.
  population_.clear();
  population_.push_back(coordinates);
  PrintFinalOutput();
  // Return the best points (with the cost removed)
  coordinates.pop_back();
  return rc;
}

int SimulatedAnnealing::ProcessFindMin(GenericCostFunction &CostFunction,
    std::vector<double> &coordinates)
{
  // Get the simulation control parameters
  const auto max_iters = options.GetMaxIterations();
  const auto max_func_evals = options.GetMaxFunctionEvaluations();
  const auto cost_tolerance = options.GetCostTolerance();
  const auto temperature_reduction_factor =
      options.GetTemperatureReductionFactor();
  const auto step_reduction_factor = options.GetStepReductionFactor();
  const auto number_of_cycles = options.GetNumberOfCycles();
  const auto number_of_temperature_loops =
      options.GetNumberOfTemperatureLoops();

  // Initialization
  auto best_point = coordinates;
  auto trial_point = coordinates;
  auto temperature = options.GetTemperature();
  double step_size = 1.0;
  PrintInitialOutput(best_point[cost_]);

  while (iterations_ < max_iters && function_evaluations_ < max_func_evals) {
    ++iterations_;
    for (auto t_loop = 0; t_loop < number_of_temperature_loops; ++t_loop) {
      for (auto cycle = 0; cycle < number_of_cycles; ++cycle) {
        // Loop over each of the coordinate directions (dimensions)
        for (auto i = 0u; i < dimensions_; ++i) {
          GenerateTrialPoint(trial_point, i);
          // If the cost of the trial point is no good, skip this loop so we can
          // generate another trial point
          if (!CalculateCost(CostFunction, trial_point)) continue;
          // Update our coordinates if the trial has a better cost. Store it
          // if it is the best cost we have seen so far.
          if (trial_point[cost_] <= coordinates[cost_]) {
            coordinates[i] = trial_point[i];
            if (trial_point[cost_] < best_point[cost_]) {
              best_point = trial_point;
            }
          }
          else {
            // If the cost of the trial point is worse than our current point,
            // apply the Metropolis Criterion to see if we want to accept it
            const auto boltzmann_probability = std::exp((coordinates[cost_] -
                trial_point[cost_]) / temperature);
            // If random number is less than Boltzmann_probability accept the
            // trial point. Otherwise, reject it and update acceptance ratios.
            if (uniform_dist_(generator_) < boltzmann_probability) {
              coordinates[i] = trial_point[i];
            }
            else {
              trial_point[i] = coordinates[i];
              acceptance_ratios_[i] -= 1.0 /
                  static_cast<double>(number_of_cycles);
            }
          }  // metropolis
        }  // i
      }  // cycles
      UpdateStepSizes();
    }  // temperature loops

    // Check for convergence: Cost has converged
    if (best_point[cost_] < cost_tolerance) break;
    // Check for convergence: Cost is not changing
    const auto converge_tol = cost_tolerance * (best_point[cost_]);
    if ((best_point[cost_] < previous_best_cost_) &&
        (previous_best_cost_ - best_point[cost_]) < converge_tol) break;

    // Reset step size
    ResetStepSizes(step_size);
    // Reduce step factor
    step_size *= step_reduction_factor;
    // Reduce temperature as per the annealing/cooling schedule
    temperature *= temperature_reduction_factor;
//    // Below is a more complex schedule suggested by Inger for his very fast
//    // annealing, assuming all dimensions have the same schedule. However,
//    // testing has shown that it is better at getting the right answer
//    // more consistently, but it is very slow.
//    //
//    const double n = std::log(static_cast<double>(max_iters));
//    const double m = -std::log(cost_tolerance);
//    const double c = m * std::exp(-n / dimensions_);
//    temperature = initial_temperature * std::exp(-c *
//        std::pow(iterations_, 1.0 / dimensions_));

    coordinates = best_point;
    trial_point = best_point;
    previous_best_cost_ = best_point[cost_];
    PrintIterationOutput(best_point[cost_]);
  }  // main loop

  if (function_evaluations_ >= max_func_evals) return 1;
  if (iterations_ >= max_iters) return 2;
  return 0;
}

}  // namespace Unfit
