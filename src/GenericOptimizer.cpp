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
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <vector>
#include "Bounds.hpp"
#include "Options.hpp"
#include "GenericCostFunction.hpp"
#include "GenericOptimizer.hpp"

namespace Unfit
{
GenericOptimizer::GenericOptimizer()
    : bounds {}
    , options {}
    , population_ {}
    , random_engines_ {}
    , function_evaluations_ {0u}
    , iterations_ {0u}
    , is_population_based_ {true}
{}

GenericOptimizer::~GenericOptimizer()
{}

void GenericOptimizer::ResetGenericOptimizer()
{
  bounds.ResetBounds();
  options.ResetOptions();
  population_.clear();
  random_engines_.clear();
  function_evaluations_ = 0u;
  iterations_ = 0u;
  // NOTE: Do not reset is_population_based_
}

double GenericOptimizer::GetCost(std::size_t index) const noexcept
{
  if (index < population_.size()) return population_[index].back();
  return std::numeric_limits<double>::infinity();
}

bool GenericOptimizer::GetIsPopulationBased() const noexcept
{
  return is_population_based_;
}

std::size_t GenericOptimizer::GetNumberOfIterations() const noexcept
{
  return iterations_;
}

std::size_t GenericOptimizer::GetNumberOfFunctionEvaluations() const noexcept
{
  return function_evaluations_;
}

std::vector<std::vector<double>> GenericOptimizer::GetPopulation() const
{
  return population_;
}

std::vector<double> GenericOptimizer::GetSolution(std::size_t index) const
{
  if (index < population_.size()) {
    auto member = population_[index];
    member.pop_back();  // remove the cost
    return member;
  }
  return {};
}

void GenericOptimizer::SetPopulation(
    const std::vector<std::vector<double>> &population)
{
  options.SetUserSetPopulation(true);
  options.SetPopulationSize(population.size());
  population_ = population;
}

bool GenericOptimizer::CalculateCost(GenericCostFunction &CostFunction,
    std::vector<double> &x)
{
  const auto residuals = CostFunction(x);
  ++function_evaluations_;
  // First check to see if the cost function failed before calculating
  if (!std::isfinite(residuals[0])) {
    x.back() = std::numeric_limits<double>::infinity();
    return false;
  }
  x.back() = std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0);
  return std::isfinite(x.back());
}

void GenericOptimizer::GeneratePopulation(GenericCostFunction &CostFunction,
    std::size_t dimensions)
{
  const auto population_size = options.GetPopulationSize();
  if (dimensions == 0u || population_size == 0u) {
    options.SetPopulationSize(0u);
    population_.clear();
    random_engines_.clear();
    return;
  }
  if (random_engines_.size() != population_size) GenerateRandomEngines();

  // Create a uniform distribution for each variable
  bounds.SetNumberOfBounds(dimensions);
  std::vector<std::uniform_real_distribution<double>> distributions;
    for (auto i = 0u; i < dimensions; ++i) {
    const auto lb = bounds.GetLowerBound(i);
    const auto ub = bounds.GetUpperBound(i);
    distributions.emplace_back(std::uniform_real_distribution<double>(lb, ub));
  }
  // Generate the population, in parallel where available
  std::vector<double> member(dimensions + 1);
  population_.assign(population_size, member);
  const auto use_multi_threaded = options.GetUseMultiThreaded();
  if (!use_multi_threaded) {
    // Serial implementation
    for (auto mem = 0u; mem < population_size; ++mem) {
      population_[mem] = GeneratePopulationMember(CostFunction, distributions,
          dimensions, mem);
    }
  }
  else {
    // Threaded (parallel) implementation
    std::vector<std::future<std::vector<double>>> population_futures;
    for (auto mem = 0u; mem < population_size; ++mem) {
      population_futures.push_back(std::async(std::launch::async,
          &GenericOptimizer::GeneratePopulationMember, this,
          std::ref(CostFunction), std::ref(distributions), dimensions, mem));
    }
    for (auto mem = 0u; mem < population_size; ++mem) {
      population_[mem] = population_futures[mem].get();
    }
  }
}

void GenericOptimizer::GenerateRandomEngines()
{
  const auto random_seed = options.GetRandomSeed();
  std::mt19937 engine_seeder(random_seed);
  std::uniform_int_distribution<int> seed_distribution(0,
      std::numeric_limits<int>::max());
  random_engines_.assign(options.GetPopulationSize(), std::mt19937 {});
  for (auto &random_engine : random_engines_) {
    random_engine.seed(seed_distribution(engine_seeder));
  }
}

bool GenericOptimizer::IsConverged(const std::vector<double> &best_member,
    std::size_t truncated_index) const
{
  // Cost convergence - always do first as geometric tolerance is more work
  auto cost_tolerance = options.GetCostTolerance();
  const auto best_cost = best_member.back();
  // If the best cost is good enough we consider it converged
  if (fabs(best_cost) < cost_tolerance) return true;

  // Check the worst to see if all points have converged to the same cost
  double worst_cost = 0.0;
  if (truncated_index >= population_.size()) truncated_index = 0u;
  if (truncated_index == 0u) {
    const auto it = std::max_element(begin(population_), end(population_),
        [](const std::vector<double> &x, const std::vector<double> &y) {
          return y.back() > x.back();
        });
    worst_cost = (*it).back();
  }
  else {
    worst_cost = population_[truncated_index].back();
  }
  // Switch to a relative tolerance if the best cost is > 1
  if (fabs(best_cost) > 1.0) cost_tolerance *= fabs(best_cost);
  if (fabs(best_cost - worst_cost) > cost_tolerance) return false;

  // Compare the distance between the best member and each of the other members
  // for geometric convergence
  const auto dimensions = best_member.size() - 1;
  std::vector<double> geom_tol(dimensions, options.GetGeometricTolerance());
  // Switch to a relative tolerance if the best coordinate is > 1
  for (auto i = 0u; i < dimensions; ++i) {
    if (fabs(best_member[i]) > 1.0) geom_tol[i] *= fabs(best_member[i]);
  }
  // Now do the comparisons
  auto last_member = population_.size() - 1;
  if (truncated_index != 0u) last_member = truncated_index;
  for (auto n = 0u; n <= last_member; ++n) {
    for (auto i = 0u; i < dimensions; ++i) {
      if (fabs(best_member[i] - population_[n][i]) > geom_tol[i]) return false;
    }
  }
  // If both the cost and geometric tests pass, we have a converged solution
  return true;
}

void GenericOptimizer::PrintInitialOutput(double best_cost) const
{
  const auto output_level = options.GetOutputLevel();
  if (output_level > 0) {
    std::cout << "Iterations";
    if (output_level > 1) {
      std::cout << "  Functions          Cost";
    }
    std::cout << std::endl;
    PrintIterationOutput(best_cost);
  }
}

void GenericOptimizer::PrintIterationOutput(double best_cost) const
{
  const auto output_level = options.GetOutputLevel();
  if (output_level > 0) {
    std::cout.width(10);
    std::cout << iterations_ << " ";
    if (output_level > 1) {
      std::cout.width(10);
      std::cout << function_evaluations_ << " ";
      std::cout.width(13);
      std::cout << std::setprecision(6) << std::scientific << best_cost;
    }
    std::cout << std::endl;
  }
}

void GenericOptimizer::PrintFinalOutput() const
{
  const auto output_level = options.GetOutputLevel();
  if (output_level > 2) {
    const auto population_size = population_.size();
    const auto dimensions = population_[0].size() - 1;
    std::cout << std::endl << "Final Results:" << std::endl;
    std::cout << std::setprecision(6);
    std::cout << std::scientific;
    for (auto k = 0u; k < population_size; ++k) {
      std::cout << " x[";
      std::cout.width(3);
      std::cout << std::right << k + 1 << "] = (";
      for (auto l = 0u; l < dimensions - 1; ++l) {
        std::cout.width(13);
        std::cout << std::right << population_[k][l] << " , ";
      }
      std::cout.width(13);
      std::cout << std::right << population_[k][dimensions - 1]
          << ") = " << std::right << population_[k].back()
          << std::endl;
    }
  }
}

void GenericOptimizer::SortPopulation() noexcept
{
  std::sort(begin(population_), end(population_),
      [](const std::vector<double> &x, const std::vector<double> &y) {
        return y.back() > x.back();
      });
}

std::vector<double> GenericOptimizer::GeneratePopulationMember(
    GenericCostFunction &CostFunction,
    std::vector<std::uniform_real_distribution<double>> &distributions,
    std::size_t dimensions, std::size_t mem)
{
  std::vector<double> member(dimensions + 1);
  do {
    for (auto dimension = 0u; dimension < dimensions; ++dimension) {
      member[dimension] = distributions[dimension](random_engines_[mem]);
    }
  } while (!CalculateCost(CostFunction, member));
  return member;
}

}  // namespace Unfit


