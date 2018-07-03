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
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <numeric>
#include <limits>
#include "NelderMead.hpp"

namespace Unfit
{
// Default construtor: keep Reset in sync with this if you change it
NelderMead::NelderMead()
  : contract_(),
    centroid_(),
    reflect_(),
    expand_(),
    dimensions_(0),
    best_vertex_(0),
    worst_vertex_(0),
    next_worst_vertex_(0),
    cost_(0),
    restart_best_cost_(0.0),
    process_(Restarted)
{
  is_population_based_ = false;
}

void NelderMead::Reset()
{
  ResetGenericOptimizer();
  contract_.clear();
  centroid_.clear();
  reflect_.clear();
  expand_.clear();
  dimensions_ = 0;
  best_vertex_ = 0;
  worst_vertex_ = 0;
  next_worst_vertex_ = 0;
  cost_ = 0;
  restart_best_cost_ = 0.0;
  process_ = Restarted;
}

int NelderMead::GeneratePopulation(GenericCostFunction &CostFunction,
    const std::vector<double> &initial_point)
{
  // Add the initial point to the simplex and set the population size
  auto vertex = initial_point;
  vertex.push_back(0.0);  // Default cost
  population_.assign(dimensions_ + 1, vertex);

  // Stop and flag an error code of 1 if the initial point is invalid
  if (!CalculateCost(CostFunction, population_[0]) ||
      !bounds.IsWithinBounds(population_[0])) return 1;

  const double nonzero_scale(0.05);
  const double zero_scale(0.00025);
  const double zero_tolerance = 1e-16;
  // Generating the simplex vertices by perturbing the coordinates
  for (auto j = 0u; j < dimensions_; ++j) {
    // Always start with the initial point
    std::vector<double> new_vertex(population_[0]);
    double increment = 0.0;
    const double coord = initial_point[j];
    if (bounds.IsWithinBounds(j, coord-zero_scale)) increment = -zero_scale;
    if (bounds.IsWithinBounds(j, coord+zero_scale)) increment = zero_scale;
    // If the initial coordinate is not zero
    if (fabs(new_vertex[j]) >= zero_tolerance) {
      if (bounds.IsWithinBounds(j, coord*(1.0-nonzero_scale))) {
        increment = -coord*nonzero_scale;
      }
      if (bounds.IsWithinBounds(j, coord*(1.0+nonzero_scale))) {
        increment = coord*nonzero_scale;
      }
    }
    // If we couldn't fit the simplex within the bounds
    if (fabs(increment) < zero_tolerance) return 2;
    // Update the new vertex with the coordinate we have found
    new_vertex[j] = coord + increment;
    // Check we have a valid cost for our vertex
    if (!CalculateCost(CostFunction, new_vertex)) return 3;
    // Add the new valid vertex to the simplex
    population_[j+1] = new_vertex;
  }
  return 0;
}

// MB Think the loops should be in the other order so we are accessing the right
// index in the inner loop
void NelderMead::ComputeCentroid()
{
  centroid_.assign(dimensions_+1, 0);
  // i refers to the element within the vertex
  for (auto j = 0u; j < dimensions_; ++j) {
    // j refers to the vertex
    double centroid_j = 0.0;
    for (auto i = 0u; i < dimensions_+1; ++i) {
      if (i == worst_vertex_) continue;
      centroid_j += population_[i][j];
    }
    centroid_[j] = centroid_j / dimensions_;
  }
}

bool NelderMead::Reflect(GenericCostFunction &CostFunction)
{
  const auto alpha = options.GetAlpha();
  for (auto j = 0u; j < dimensions_; ++j) {
    reflect_[j] = centroid_[j] + alpha*(centroid_[j] -
        population_[worst_vertex_][j]);
  }
  if (!bounds.IsWithinBounds(reflect_)) return false;
  if (!CalculateCost(CostFunction, reflect_)) return false;
  return true;
}

bool NelderMead::Expand(GenericCostFunction &CostFunction)
{
  const auto beta = options.GetBeta();
  for (auto j = 0u; j < dimensions_; ++j) {
    expand_[j] = centroid_[j] + beta*(reflect_[j] - centroid_[j]);
  }
  if (!bounds.IsWithinBounds(expand_)) return false;
  if (!CalculateCost(CostFunction, expand_)) return false;
  return true;
}

bool NelderMead::ContractOutside(GenericCostFunction &CostFunction)
{
  const auto alpha = options.GetAlpha();
  const auto gamma = options.GetGamma();
  // Contractions can't be based on reflect_ as reflect may have failed
  for (auto j = 0u; j < dimensions_; ++j) {
    contract_[j] = centroid_[j] + gamma*alpha*(centroid_[j] -
        population_[worst_vertex_][j]);
  }
  if (!bounds.IsWithinBounds(contract_)) return false;
  if (!CalculateCost(CostFunction, contract_)) return false;
  return true;
}

bool NelderMead::ContractInside(GenericCostFunction &CostFunction)
{
  const auto alpha = options.GetAlpha();
  const auto gamma = options.GetGamma();
  // Contractions can't be based on reflect_ as reflect may have failed
  for (auto j = 0u; j < dimensions_; ++j) {
    contract_[j] = centroid_[j] - gamma*alpha*(centroid_[j] -
        population_[worst_vertex_][j]);
  }
  // There is no need to check the domain here as we are inside the simplex
  if (!CalculateCost(CostFunction, contract_)) return false;
  return true;
}

bool NelderMead::Shrink(GenericCostFunction &CostFunction)
{
  const auto delta = options.GetDelta();
  for (auto i = 0u; i < dimensions_+1; ++i) {
    if (i == best_vertex_) continue;
    for (auto j = 0u; j < dimensions_; ++j) {
      population_[i][j] = population_[best_vertex_][j] +
          delta*(population_[i][j] - population_[best_vertex_][j]);
    }
    // There is no need to check the domain here as we are inside the simplex
    if (!CalculateCost(CostFunction, population_[i])) return false;
  }
  return true;
}

bool NelderMead::IsDegenerate()
{
  // Check to see if any of the vertex coordinates is further away from the
  // centroid than the tolerance for degeneracy
  const auto degenerate_tol = options.GetDegenerateTolerance();
  for (auto j = 0u; j < dimensions_; ++j) {
    if (fabs(population_[worst_vertex_][j]-centroid_[j]) > degenerate_tol)
        return false;
  }
  // If all of the coordinates are too close to the centroid the simplex is
  // degenerate
  return true;
}

int NelderMead::FindMin(GenericCostFunction &CostFunction,
    std::vector<double> &coordinates)
{
  // Check that the initialpoint is not zero
  if (coordinates.empty()) return -4;
  // Initialise if this is the first time
  if (dimensions_ == 0) {
    dimensions_ = coordinates.size();
    InitialiseVectors();
    cost_ = dimensions_;
  }
  if (options.GetUseAdaptiveParameters()) UseAdaptiveParameters();

  // Check whether the user has given us a population already, otherwise
  // generate a uniform random population
  if (options.GetUserSetPopulation()) {
    // Check the population is not empty
    if (population_.empty()) return -1;
    // Check the population members are not empty
    if (population_[0].empty()) return -2;
    // Check the minimum size requirement is met
    if (population_.size() != population_[0].size()) return -3;
    // Check the initial guess is consistent with the population
    for (auto member : population_) {
      if (member.size() != dimensions_ + 1u) return -4;
    }
    // Check the cost of each vertex
    for (auto member : population_) {
      if (!std::isfinite(member.back())) return -5;
    }
  }
  else {
    // Generate the initial simplex
    auto rc = GeneratePopulation(CostFunction, coordinates);
    if (rc != 0) return -rc;
  }
  restart_best_cost_ = population_[best_vertex_][cost_];
  // Perform the optimisation
  auto rc = ProcessFindMin(CostFunction);
  // Return the best point, without the cost -> end()-1
  std::vector<double>::iterator first = population_[best_vertex_].begin();
  std::vector<double>::iterator last = population_[best_vertex_].end() - 1;
  coordinates.assign(first, last);

  return rc;
}

int NelderMead::RegeneratePopulation(GenericCostFunction &CostFunction)
{
  auto best = population_[best_vertex_];
  best.pop_back();  // Remove the cost as GeneratePopulation will push it back
  return GeneratePopulation(CostFunction, best);
}

int NelderMead::ProcessFindMin(GenericCostFunction &CostFunction)
{
  // The initial simplex counts as one iteration in Matlab, so we will do the
  // same thing here.
  ++iterations_;
  PrintInitialOutput(population_[best_vertex_][cost_]);
  SortPopulation();
  best_vertex_ = 0;
  worst_vertex_ = dimensions_;
  next_worst_vertex_ = dimensions_ - 1;
  // Make sure we only call these functions once
  const auto cost_tol = options.GetCostTolerance();
  const auto max_iters = options.GetMaxIterations();
  const auto max_func_evals = options.GetMaxFunctionEvaluations();

  while (iterations_ < max_iters && function_evaluations_ < max_func_evals) {
    ++iterations_;

    ComputeCentroid();
    // Calculate the reflection point and its cost
    if (!Reflect(CostFunction)) {
      // Unable to find a valid reflect point, so set the cost of a
      // reflection such that it forces a contraction and will not be adopted
      reflect_[cost_] = population_[worst_vertex_][cost_] + cost_tol;
    }
    // If the cost of the reflect point is better than the 2nd worst cost
    if (reflect_[cost_] < population_[next_worst_vertex_][cost_]) {
      // If the cost of the reflect point is not as good as the best cost,
      // replace the worst vertex with the reflect point
      if (population_[best_vertex_][cost_] < reflect_[cost_]) {
        population_[worst_vertex_] = reflect_;
        process_ = Reflected;
      }
      else {
        // If the reflect cost is better than the best cost, calculate the
        // position & cost of the expansion point
        if (!Expand(CostFunction)) {
          // Unable to find a valid expansion point, so set the cost of an
          // expansion such that it forces a reflect
          expand_[cost_] = reflect_[cost_] + cost_tol;
        }
        // If the expand cost is better than the reflect cost, replace the
        // vertex with the worst cost with the expand point, otherwise replace
        // the vertex with the worst cost with the reflect point.
        if (expand_[cost_] < reflect_[cost_]) {
          population_[worst_vertex_] = expand_;
          process_ = Expanded;
        }
        else {
          population_[worst_vertex_] = reflect_;
          process_ = Reflected;
        }
      }
    }
    else {
      // If the cost of the reflect point is better than the worst cost then
      // perform an outside contraction, otherwise perform an inside contraction
      if (reflect_[cost_] < population_[worst_vertex_][cost_]) {
        if (!ContractOutside(CostFunction)) {
          // Unable to find a valid contraction point, so set the cost of a
          // contraction such that it forces a shrink
          contract_[cost_] = population_[worst_vertex_][cost_] + cost_tol;
        }
        // If the cost of the outside contract point is better than the worst
        // cost, replace the worst point with the outside contract point
        if (contract_[cost_] < population_[worst_vertex_][cost_]) {
          population_[worst_vertex_] = contract_;
          process_ = ContractedOut;
        }
        else {
          // We can't find a replacement for the worst vertex, so shrink
          if (!Shrink(CostFunction)) return 3;  // Unable to find valid simplex
          process_ = Shrunk;
        }
      }
      else {
        if (!ContractInside(CostFunction)) {
          // Unable to find a valid contraction point, so set the cost of a
          // contraction such that it forces a shrink
          contract_[cost_] = population_[worst_vertex_][cost_] + cost_tol;
        }
        // If we contracted inside and the cost is better than the reflect cost,
        // replace the vertex with the worst cost with the contract point
        if (contract_[cost_] <= reflect_[cost_]) {
          population_[worst_vertex_] = contract_;
          process_ = ContractedIn;
        }
        else {
          // We can't find a replacement for the worst vertex, so shrink
          if (!Shrink(CostFunction)) return 3;  // Unable to find valid simplex
          process_ = Shrunk;
        }
      }
    }
    SortPopulation();
    if (IsConverged(population_[best_vertex_])) break;
    if (process_ == ContractedIn || process_ == ContractedOut ||
        process_ == Shrunk) {
      if (IsDegenerate()) {
        // If we made progress after the last restart it is worth trying again
        if (population_[best_vertex_][cost_] < restart_best_cost_) {
          // Generate a new simplex about the best point but return an error
          // if we can't regenerate the simplex
          if (RegeneratePopulation(CostFunction) != 0) return 3;
          SortPopulation();
          restart_best_cost_ = population_[best_vertex_][cost_];
          process_ = Restarted;
        }
        else {
          // If we have already restarted and didn't make any progress then
          // there is no point in restarting again as we will just get into
          // an infinite loop. Just flag this.
          return 4;
        }
      }
    }
    PrintIterationOutput(population_[best_vertex_][cost_]);
  }
  PrintFinalOutput();
  if (function_evaluations_ >= max_func_evals) return 1;
  if (iterations_ >= max_iters) return 2;
  return 0;
}

void NelderMead::UseAdaptiveParameters()
{
  double alpha = 1.0;
  double beta = 1.0 + (2.0/dimensions_);
  double gamma = 0.75 - (1.0/(2.0*dimensions_));
  double delta = 1.0 - (1.0/dimensions_);
  options.SetNelderMeadStepSizes(alpha, beta, delta, gamma);
}

void NelderMead::InitialiseVectors()
{
  // A vector to store the centroid. All of these are dimensions_ + 1 because we
  // are currently storing the cost with the coordinates to make sorting easy.
  centroid_.assign(dimensions_ + 1, 0);
  // A vector to store the contract point
  contract_.assign(dimensions_ + 1, 0);
  // A vector to store the expand point
  expand_.assign(dimensions_ + 1, 0);
  // A vector to store the reflect point
  reflect_.assign(dimensions_ + 1, 0);
  // Set the default bounds (does not overwrite those that exist)
  bounds.SetNumberOfBounds(dimensions_);
}

void NelderMead::PrintInitialOutput(double best_cost) const
{
  const auto output_level = options.GetOutputLevel();
  if (output_level > 0) {
    std::cout << "Iterations";
    if (output_level > 1) {
      std::cout << "  Functions          Cost  Operation";
    }
    std::cout << std::endl;
    PrintIterationOutput(best_cost);
  }
}

void NelderMead::PrintIterationOutput(double best_cost) const
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
      switch (process_) {
        case Reflected:
          std::cout << "  reflect";
          break;
        case Expanded:
          std::cout << "  expand";
          break;
        case ContractedIn:
          std::cout << "  contract inside";
          break;
        case ContractedOut:
          std::cout << "  contract outside";
          break;
        case Shrunk:
          std::cout << "  shrink";
          break;
        case Restarted:
          std::cout << "  restart";
          break;
        default:
          std::cout << "  WARNING: Operation is undefined";
      }
    }
    std::cout << std::endl;
  }
}

}  // namespace Unfit
