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
#include <algorithm>
#include <future>
#include <iostream>
#include <vector>
#include "Bounds.hpp"
#include "DifferentialEvolution.hpp"
#include "Options.hpp"

namespace Unfit
{
// Default construtor: keep Reset in sync with this if you change it
DifferentialEvolution::DifferentialEvolution()
  : new_population_ {}
  , best_member_ {}
  , dimensions_ {0}
  , cost_ {0}
{}

DifferentialEvolution::~DifferentialEvolution() = default;

void DifferentialEvolution::Reset()
{
  ResetGenericOptimizer();
  new_population_.clear();
  best_member_.clear();
  dimensions_ = 0;
  cost_ = 0;
}

int DifferentialEvolution::FindMin(GenericCostFunction &CostFunction,
    std::vector<double> &coordinates)
{
  if (coordinates.empty()) return -1;
  dimensions_ = static_cast<unsigned>(coordinates.size());
  // the value of the size of dimensions_ is the same as the value of
  // of the index for accessing the cost in a given vector i.e. cost is
  // always the last vector element. So, we assign dimensions_ to cost_.
  cost_ = dimensions_;

  // Check whether the user has given us a population already, otherwise
  // generate a uniform random population
  if (options.GetUserSetPopulation()) {
    // Check the population is not empty
    if (population_.empty() || population_[0].empty()) return -2;
    // Check the minimum size requirement is met
    if (options.GetPopulationSize() < 6u) return -2;
    // Check the initial guess is consistent with the population
    // +1 because cost is appended to the end of member vector
    for (auto member : population_) {
      if (member.size() != dimensions_ + 1u) return -2;
    }
    GenerateRandomEngines();
  }
  else {
    // Set the population size to be 10x dimensions unless the user has set this
    if (!options.GetUserSetPopulationSize()) {
      options.SetPopulationSize(10u * dimensions_);
    }
    // If the user has set the size, make sure it is >=6 as some strategies
    // require this.
    else if (options.GetPopulationSize() < 6u) {
      options.SetPopulationSize(6u);
    }
    GeneratePopulation(CostFunction, dimensions_);
  }

  // If requested and possible, add the initial coordinates to the population
  if (options.GetAddInitialToPopulation()) {
    auto member = coordinates;
    member.push_back(0.0);  // to store the cost
    if (CalculateCost(CostFunction, member)) {
      // Initial coordinates result in a valid cost, overwrite the first member
      population_[0] = member;
    }
    else {
      std::cout << "WARNING: Failed to add the initial coordinates to the "
          << "population (DE)" << std::endl;
    }
  }

  // Find and store the best member and its cost
  auto it = std::min_element(begin(population_), end(population_),
      [](const std::vector<double> &x, const std::vector<double> &y) {
        return y.back() > x.back();
      });
  best_member_ = *it;
  new_population_ = population_;
  ++iterations_;  // The initial population is counted as one iteration

  int rc = ProcessFindMin(CostFunction);
  // Return the best member of the population (with the cost removed)
  coordinates = best_member_;
  coordinates.pop_back();
  return rc;
}

int DifferentialEvolution::ProcessFindMin(GenericCostFunction &CostFunction)
{
  const auto max_iters = options.GetMaxIterations();
  const auto max_func_evals = options.GetMaxFunctionEvaluations();
  const auto population_size = options.GetPopulationSize();
  const auto use_multi_threaded = options.GetUseMultiThreaded();
  PrintInitialOutput(best_member_[cost_]);

  while (iterations_ < max_iters && function_evaluations_ < max_func_evals) {
    if (!use_multi_threaded) {
      // Serial implementation
      for (auto i = 0u; i < population_size; ++i) {
        new_population_[i] = NewPopulationMember(CostFunction, i);
      }
    }
    else {
      // Threaded (parallel) implementation
      std::vector<std::future<std::vector<double>>> new_population_futures;
      for (auto i = 0u; i < population_size; ++i) {
        new_population_futures.push_back(std::async(std::launch::async,
            &DifferentialEvolution::NewPopulationMember, this,
            std::ref(CostFunction), i));
      }
      for (auto i = 0u; i < population_size; ++i) {
        new_population_[i] = new_population_futures[i].get();
      }
    }

    population_.swap(new_population_);
    auto it = std::min_element(begin(population_), end(population_),
        [](const std::vector<double> &x, const std::vector<double> &y) {
          return y.back() > x.back();
        });
    best_member_ = *it;
    ++iterations_;
    PrintIterationOutput(best_member_[cost_]);
    if (IsConverged(best_member_)) break;
  }
  PrintFinalOutput();
  if (function_evaluations_ >= max_func_evals) return 1;
  if (iterations_ >= max_iters) return 2;
  return 0;
}

std::vector<double> DifferentialEvolution::NewPopulationMember(
    GenericCostFunction &CostFunction, unsigned member)
{
  auto trial = GenerateTrialMember(member);
  if (options.GetUseHardBounds()) bounds.ClampWithinBounds(trial);
  if (!CalculateCost(CostFunction, trial)) {
    return population_[member];
  }
  else if (trial[cost_] <= population_[member][cost_]) {
    return trial;
  }
  else {
    return population_[member];
  }
}

std::vector<double> DifferentialEvolution::GenerateTrialMember(unsigned mem)
{
  const auto population_size = options.GetPopulationSize();
  std::uniform_int_distribution<unsigned> random_member(0u, population_size-1);
  std::uniform_int_distribution<unsigned> random_parameter(0u, dimensions_-1);
  std::uniform_real_distribution<double> random_probability(0.0, 1.0);
  const auto strategy = options.GetStrategy();
  const auto F = options.GetWeightingFactor();
  const auto CR = options.GetCrossOver();

  auto r1 = 0u;  // Indices for randomly selected members
  auto r2 = 0u;
  auto r3 = 0u;
  auto r4 = 0u;
  auto r5 = 0u;
  do {
    r1 = random_member(random_engines_[mem]);
  } while (r1 == mem);
  do {
    r2 = random_member(random_engines_[mem]);
  } while ((r2 == mem) || (r2 == r1));
  do {
    r3 = random_member(random_engines_[mem]);
  } while ((r3 == mem) || (r3 == r1) || (r3 == r2));
  do {
    r4 = random_member(random_engines_[mem]);
  } while ((r4 == mem) || (r4 == r1) || (r4 == r2) || (r4 == r3));
  do {
    r5 = random_member(random_engines_[mem]);
  } while ((r5 == mem) || (r5 == r1) || (r5 == r2) || (r5 == r3) || (r5 == r4));

  auto trial = population_[mem];
  auto n = random_parameter(random_engines_[mem]);
  if (strategy == 1u) {  // DE/best/1/exp
    auto j = 0u;
    do {
      trial[n] = best_member_[n] + F*(population_[r2][n] - population_[r3][n]);
      n = (n + 1) % dimensions_;
    } while ((random_probability(random_engines_[mem]) < CR) &&
        (++j < dimensions_));
  }
  else if (strategy == 2u) {  // DE/rand/1/exp
    auto j = 0u;
    do {
      trial[n] =population_[r1][n] + F*(population_[r2][n] -
          population_[r3][n]);
      n = (n + 1) % dimensions_;
    } while ((random_probability(random_engines_[mem]) < CR) &&
        (++j < dimensions_));
  }
  else if (strategy == 3u) {  // DE/rand-to-best/1/exp
    auto j = 0u;
    do {
      trial[n] = trial[n] + F*(best_member_[n] - trial[n]) +
          F*(population_[r1][n] - population_[r2][n]);
      n = (n + 1) % dimensions_;
    } while ((random_probability(random_engines_[mem]) < CR) &&
        (++j < dimensions_));
  }
  else if (strategy == 4u) {  // DE/best/2/exp
    auto j = 0u;
    do {
      trial[n] = best_member_[n] + F*(population_[r1][n] +
          population_[r2][n] - population_[r3][n] -
          population_[r4][n]);
      n = (n + 1) % dimensions_;
    } while ((random_probability(random_engines_[mem]) < CR) &&
        (++j < dimensions_));
  }
  else if (strategy == 5u) {  // DE/rand/2/exp
    auto j = 0u;
    do {
      trial[n] = population_[r5][n] + F*(population_[r1][n] +
          population_[r2][n] - population_[r3][n] - population_[r4][n]);
      n = (n + 1) % dimensions_;
    } while ((random_probability(random_engines_[mem]) < CR) &&
        (++j < dimensions_));
  }
  else if (strategy == 6u) {  // DE/best/1/bin
    for (auto j = 0u; j < dimensions_; ++j) {
      if (random_probability(random_engines_[mem]) < CR) {
        trial[j] = best_member_[j] + F*(population_[r2][j] -
            population_[r3][j]);
      }
    }
    trial[n] = best_member_[n] + F*(population_[r2][n] - population_[r3][n]);
  }
  else if (strategy == 7u) {  // DE/rand/1/bin
    for (auto j = 0u; j < dimensions_; ++j) {
      if (random_probability(random_engines_[mem]) < CR) {
        trial[j] = population_[r1][j] + F*(population_[r2][j] -
            population_[r3][j]);
      }
    }
    trial[n] = population_[r1][n] + F*(population_[r2][n] - population_[r3][n]);
  }
  else if (strategy == 8u) {  // DE/rand-to-best/1/bin
    for (auto j = 0u; j < dimensions_; ++j) {
      if (random_probability(random_engines_[mem]) < CR) {
        trial[j] = trial[j] + F*(best_member_[j] - trial[j]) +
            F*(population_[r1][j] - population_[r2][j]);
      }
    }
    trial[n] = trial[n] + F*(best_member_[n] - trial[n]) +
        F*(population_[r1][n] - population_[r2][n]);
  }
  else if (strategy == 9u) {  // DE/best/2/bin
    for (auto j = 0u; j < dimensions_; ++j) {
      if (random_probability(random_engines_[mem]) < CR) {
        trial[j] = best_member_[j] + F*(population_[r1][j] +
            population_[r2][j] - population_[r3][j] - population_[r4][j]);
      }
    }
    trial[n] = best_member_[n] + F*(population_[r1][n] +
        population_[r2][n] - population_[r3][n] - population_[r4][n]);
  }
  else {  // DE/rand/2/bin
    for (auto j = 0u; j < dimensions_; ++j) {
      if (random_probability(random_engines_[mem]) < CR) {
        trial[j] = population_[r5][j] + F*(
            population_[r1][j] + population_[r2][j] -
            population_[r3][j] - population_[r4][j]);
      }
    }
    trial[n] = population_[r5][n] + F*(population_[r1][n] +
        population_[r2][n] - population_[r3][n] - population_[r4][n]);
  }
  return trial;
}

}  // namespace Unfit
