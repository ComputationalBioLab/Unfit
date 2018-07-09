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
#include <random>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

#include "Bounds.hpp"
#include "GeneticAlgorithm.hpp"
#include "Options.hpp"

// NOTE: May need to provide some sort of iteration counter for our do-while
//       loops to prevent infinite loops in the case of e.g. an always nan cost
//       function

namespace Unfit
{
// Default construtor: keep Reset in sync with this if you change it
GeneticAlgorithm::GeneticAlgorithm()
  : n_keep_ {0}
  , dimensions_ {0}
  , ranks_ {}
  , generator_ {0}
  , distributions_ {}
{
  options.SetGamma(0.2);
}

void GeneticAlgorithm::Reset()
{
  ResetGenericOptimizer();
  n_keep_ = 0;
  dimensions_ = 0;
  ranks_.clear();
  generator_.seed(0);
  distributions_.clear();
  options.SetGamma(0.2);
}

void GeneticAlgorithm::InitialiseBounds()
{
  // Make sure we have bounds on each variable
  bounds.SetNumberOfBounds(dimensions_);
  // Create a distribution for each variable
  for (auto i = 0u; i < dimensions_; ++i) {
    const auto lb = bounds.GetLowerBound(i);
    const auto ub = bounds.GetUpperBound(i);
    distributions_.emplace_back(std::uniform_real_distribution<double>(lb, ub));
  }
}

void GeneticAlgorithm::MutateGenes(GenericCostFunction &CostFunction)
{
  const auto population_size = options.GetPopulationSize();
  const auto elitism = options.GetElitism();
  if (elitism >= population_size) return;  // All elite = no mutations
  std::uniform_int_distribution<unsigned> eligible_chromosomes(elitism,
      population_size - 1u);
  std::uniform_int_distribution<unsigned> eligible_genes(0u, dimensions_ - 1u);
  const auto gamma = options.GetGamma();
  const auto number_of_mutations = static_cast<unsigned>(gamma *
      population_size * dimensions_);
  for (auto i = 0u; i < number_of_mutations; ++i) {
    auto chromosome = eligible_chromosomes(generator_);
    auto gene = eligible_genes(generator_);
    do {
      population_[chromosome][gene] = distributions_[gene](generator_);
    } while (!CalculateCost(CostFunction, population_[chromosome]));
  }
}

void GeneticAlgorithm::CalculateRanks()
{
  const auto population_size = options.GetPopulationSize();
  const auto survival_rate = options.GetSurvivalRate();
  n_keep_ = static_cast<unsigned>(population_size * survival_rate);
  if (n_keep_ >= population_size) return;  // Keep all the parents
  if (n_keep_ < 2u) n_keep_ = 2u;  // Keep one pair
  ranks_.resize(n_keep_);
  for (auto i = 0u; i < n_keep_; ++i) ranks_[i] = n_keep_ - i;
  std::partial_sum(begin(ranks_), end(ranks_), begin(ranks_));
  for (auto &r : ranks_) r /= ranks_.back();
}

std::pair<unsigned, unsigned> GeneticAlgorithm::GetMatingPair()
{
  const auto population_size = options.GetPopulationSize();
  if (n_keep_ == population_size) return {0u, 0u};
  std::uniform_real_distribution<double> weights(0.0, 1.0);
  // Select the first parent
  unsigned parent_1 = 0u;
  auto random_num = weights(generator_);
  while (random_num > ranks_[parent_1]) ++parent_1;
  // Select the second parent
  unsigned parent_2 = 0u;
  do {
    random_num = weights(generator_);
    while (random_num > ranks_[parent_2]) ++parent_2;
    if (parent_1 == (ranks_.size()-1)) parent_2 = 0u;  // Prevent inf loop
  } while (parent_1 == parent_2);
  return std::make_pair(parent_1, parent_2);
}

void GeneticAlgorithm::CrossOver(const std::vector<double> &parent_1,
    const std::vector<double> &parent_2, std::vector<double> &offspring_1,
    std::vector<double> &offspring_2)
{
  std::uniform_real_distribution<double> weights(0.0, 1.0);
  offspring_1 = parent_1;
  offspring_2 = parent_2;
  const auto cross_over_point = static_cast<unsigned>(weights(generator_) *
      (dimensions_ - 1) + 0.5);  // Indices 0-> dim-1; Round, not truncate
  const auto blend = weights(generator_) * (parent_1[cross_over_point] -
      parent_2[cross_over_point]);
  offspring_1[cross_over_point] -= blend;
  offspring_2[cross_over_point] += blend;
  if (cross_over_point != (dimensions_ - 1)) {
    // Swap genes after the crossover point
    for (auto i = cross_over_point + 1; i < dimensions_; ++i) {
      offspring_1[i] = parent_2[i];
      offspring_2[i] = parent_1[i];
    }
  }
  else {
    // Instead of swapping the genes before the crossover point, just swap the
    // last (blended) entry. These are equivalent because the order of the
    // offspring does not matter
    std::swap(offspring_1[dimensions_ - 1], offspring_2[dimensions_ - 1]);
  }
}

void GeneticAlgorithm::GeneratePopulation(GenericCostFunction &CostFunction)
{
  const auto population_size = options.GetPopulationSize();
  if (dimensions_ == 0u || population_size == 0u) {
    population_.clear();
    return;
  }
  std::vector<double> chromosome(dimensions_ + 1);
  population_.assign(population_size, chromosome);
  for (auto &individual : population_) {
    do {
      for (auto gene = 0u; gene < dimensions_; ++gene) {
        individual[gene] = distributions_[gene](generator_);
      }
    } while (!CalculateCost(CostFunction, individual));
  }
}

void GeneticAlgorithm::Reproduce(GenericCostFunction &CostFunction)
{
  const auto population_size = options.GetPopulationSize();
  for (auto i = n_keep_; i < population_size; i += 2) {
    const auto parents = GetMatingPair();
    std::vector<double> offspring_1;
    std::vector<double> offspring_2;
    do {
      CrossOver(population_[parents.first], population_[parents.second],
          offspring_1, offspring_2);
    } while (!CalculateCost(CostFunction, offspring_1) ||
        !CalculateCost(CostFunction, offspring_2));
    population_[i] = std::move(offspring_1);
    if ((i + 1) < population_size) population_[i+1] = std::move(offspring_2);
  }
}

int GeneticAlgorithm::FindMin(GenericCostFunction &CostFunction,
    std::vector<double> &coordinates)
{
  if (coordinates.empty()) return -1;
  dimensions_ = static_cast<unsigned>(coordinates.size());
  generator_.seed(options.GetRandomSeed());
  // Check whether the user has given us a population already, otherwise
  // generate a uniform random population
  if (options.GetUserSetPopulation()) {
    // Check the population is not empty
    if (population_.empty() || population_[0].empty()) return -2;
    // Check the minimum size requirement is met
    if (options.GetPopulationSize() < 3u) return -2;
    // Check the initial guess is consistent with the population
    // +1 because cost is appended to the end of member vector
    for (auto member : population_) {
      if (member.size() != dimensions_ + 1u) return -2;
    }
    InitialiseBounds();
    GenerateRandomEngines();
  }
  else {
    // Set the population size to be 10x dimensions unless the user has set it.
    // If the user has set the size, make sure it is >=3 as a minimum
    if (!options.GetUserSetPopulationSize()) {
      options.SetPopulationSize(10u * dimensions_);
    }
    else if (options.GetPopulationSize() < 3u) {
      options.SetPopulationSize(3u);
    }
    InitialiseBounds();
    GeneratePopulation(CostFunction);
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
          << "population (GA)" << std::endl;
    }
  }

  SortPopulation();
  CalculateRanks();
  int rc = ProcessFindMin(CostFunction);
  // Return the best member of the population (with the cost removed)
  coordinates = population_[0];
  coordinates.pop_back();
  return rc;
}

int GeneticAlgorithm::ProcessFindMin(GenericCostFunction &CostFunction)
{
  ++iterations_;  // The initial population is counted as one iteration
  PrintInitialOutput(population_[0].back());
  const auto max_iters = options.GetMaxIterations();
  const auto max_func_evals = options.GetMaxFunctionEvaluations();
  while (iterations_ < max_iters && function_evaluations_ < max_func_evals) {
    ++iterations_;
    Reproduce(CostFunction);
    MutateGenes(CostFunction);
    SortPopulation();
    PrintIterationOutput(population_[0].back());
    if (IsConverged(population_[0], n_keep_ - 1)) break;
  }
  PrintFinalOutput();
  if (function_evaluations_ >= max_func_evals) return 1;
  if (iterations_ >= max_iters) return 2;
  return 0;
}

}  // namespace Unfit
