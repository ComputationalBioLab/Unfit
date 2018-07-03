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
#include <cmath>
#include <future>
#include <iostream>
#include <vector>
#include "ParticleSwarm.hpp"

namespace Unfit
{
ParticleSwarm::ParticleSwarm()
  : best_particle_ {}
  , dimensions_ {0}
  , cost_ {0}
{
  options.SetAlpha(1.0);
  options.SetBeta(0.5);
  options.SetDelta(0.99);
}

ParticleSwarm::~ParticleSwarm() = default;

void ParticleSwarm::Reset()
{
  ResetGenericOptimizer();
  best_particle_ = {};
  dimensions_ = 0;
  cost_ = 0;
  options.SetAlpha(1.0);
  options.SetBeta(0.5);
  options.SetDelta(0.99);
}

int ParticleSwarm::FindMin(GenericCostFunction &CostFunction,
    std::vector<double> &coordinates)
{
  if (coordinates.empty()) return -1;
  dimensions_ = coordinates.size();
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
    if (options.GetPopulationSize() < 2u) return -2;
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
    // Set the minimum population size if the user has set it too low
    else if (options.GetPopulationSize() < 2u) {
      options.SetPopulationSize(2u);
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
          << "population (PS)" << std::endl;
    }
  }

  // Find and store the best member and its cost
  auto it = std::min_element(begin(population_), end(population_),
      [](const std::vector<double> &x, const std::vector<double> &y) {
        return y.back() > x.back();
      });
  best_particle_ = *it;
  ++iterations_;  // The initial population is counted as one iteration

  int rc = ProcessFindMin(CostFunction);
  // Return the best member of the population (with the cost removed)
  coordinates = best_particle_;
  coordinates.pop_back();
  return rc;
}

int ParticleSwarm::ProcessFindMin(GenericCostFunction &CostFunction)
{
  const auto max_iters = options.GetMaxIterations();
  const auto max_func_evals = options.GetMaxFunctionEvaluations();
  const auto population_size = options.GetPopulationSize();
  const auto delta = options.GetDelta();
  const auto use_chaos_enhancement = options.GetUseAdaptiveParameters();
  const auto enhancement_strategy = options.GetStrategy();
  const auto use_multi_threaded = options.GetUseMultiThreaded();
  PrintInitialOutput(best_particle_[cost_]);

  while (iterations_ < (max_iters) && function_evaluations_ < max_func_evals) {
    if (!use_multi_threaded) {
      // Serial implementation
      for (auto i = 0u; i < population_size; ++i) {
        UpdatePopulationMember(CostFunction, i);
      }
    }
    else {
      // Threaded (parallel) implementation
      std::vector<std::future<void>> population_futures;
      for (auto i = 0u; i < population_size; ++i) {
        population_futures.push_back(std::async(std::launch::async,
            &ParticleSwarm::UpdatePopulationMember, this,
            std::ref(CostFunction), i));
      }
      for (auto i = 0u; i < population_size; ++i) {
        population_futures[i].get();
      }
    }

    // Find the best member in the current population
    auto it = std::min_element(begin(population_), end(population_),
        [](const std::vector<double> &x, const std::vector<double> &y) {
          return y.back() > x.back();
        });
    // Update our saved particle if it is the best of all time
    if ((*it)[cost_] < best_particle_[cost_]) best_particle_ = *it;
    ++iterations_;
    // After each iteration, alpha is reduced by a factor delta
    options.SetAlpha(options.GetAlpha()*delta);
    // After each iteration, beta is altered if an adaptive strategy is used
    if (use_chaos_enhancement) ChaosEnhancement(enhancement_strategy);
    PrintIterationOutput(best_particle_[cost_]);
    if (IsConverged(best_particle_)) break;
  }
  PrintFinalOutput();
  if (function_evaluations_ >= max_func_evals) return 1;
  if (iterations_ >= max_iters) return 2;
  return 0;
}

void ParticleSwarm::ChaosEnhancement(unsigned enhancement_strategy)
{
  auto beta = options.GetBeta();
  while (beta > 1.0) beta -= 1.0;
  while (beta < 0.0) beta += 1.0;
  const double beta_zero_tol = 1e-6;
  const double my_pi = acos(-1.0);
  const auto two_pi = 2.0*my_pi;
  switch (enhancement_strategy) {
    case 1: {  // Chebyshev Map
      // Beta cannot be 0.0, 0.5 or 1.0
      if (fabs(beta) < beta_zero_tol) beta += beta_zero_tol;
      if (fabs(beta - 0.5) < beta_zero_tol) beta -= beta_zero_tol;
      if (fabs(beta - 1.0) < beta_zero_tol) beta -= beta_zero_tol;
      auto new_beta = cos(iterations_ * acos(beta));
      options.SetBeta(fabs(new_beta));
      break;
    }
    case 2: {  // Circle Map
      const double a = 0.5;
      const double b = 0.2;
      auto temp_beta = beta + (b - (a/two_pi) * sin(two_pi * beta));
      auto new_beta = fmod(temp_beta, 1.0);
      options.SetBeta(fabs(new_beta));
      break;
    }
    case 3: {  // Gauss/Mouse Map
      // Beta cannot be 0.0
      if (fabs(beta) < beta_zero_tol) {
        options.SetBeta(beta_zero_tol);
      }
      else {
        options.SetBeta(fmod(1.0/beta, 1.0));
      }
      break;
    }
    case 4: {  // Intermittency Map
      // Beta cannot be 1.0
      if (fabs(beta - 1.0) < beta_zero_tol) beta -= beta_zero_tol;
      const double p = 0.5;
      const double epsilon = 0.001;
      const double c = (1.0 - epsilon - p) / (p * p);
      if (beta <= p) {
        options.SetBeta(epsilon + beta + c*beta*beta);
      }
      else {
        options.SetBeta((beta - p)/(1.0 - p));
      }
      break;
    }
    case 5: {  // Iterative Map
      // Beta cannot be 0.0
      if (fabs(beta) < beta_zero_tol) beta += beta_zero_tol;
      const double a = 0.111;  // Should be between 0 and 1
      options.SetBeta(sin(a * my_pi / beta));
      break;
    }
    case 6: {  // Liebovitch Map
      // Beta cannot be 0.0, 0.5 or 1.0
      if (fabs(beta) < beta_zero_tol) beta += beta_zero_tol;
      if (fabs(beta - 0.5) < beta_zero_tol) beta -= beta_zero_tol;
      if (fabs(beta - 1.0) < beta_zero_tol) beta -= beta_zero_tol;
      const double p1 = 1.0 / 3.0;
      const double p2 = 2.0 / 3.0;
      if (beta <= p1) {
        const double alpha = (p2 / p1) * (1.0 - (p2 - p1));
        options.SetBeta(alpha * beta);
      }
      else if (beta <= p2) {
        options.SetBeta((p2 - beta) / (p2 - p1));
      }
      else {
        const double delta = (1.0 / (p2 - 1.0)) * ((p2 - 1) - p1*(p2 - p1));
        options.SetBeta(1.0 - delta*(1.0 - beta));
      }
      break;
    }
    case 7: {  // Logistic Map
      const double a = 4.0;
      // Beta cannot be 0.25, 0.5, 0.75 or 1.0
      if (fabs(beta) < beta_zero_tol) beta += beta_zero_tol;
      if (fabs(beta - 0.25) < beta_zero_tol) beta += beta_zero_tol;
      if (fabs(beta - 0.5) < beta_zero_tol) beta += beta_zero_tol;
      if (fabs(beta - 0.75) < beta_zero_tol) beta -= beta_zero_tol;
      if (fabs(beta - 1.0) < beta_zero_tol) beta -= beta_zero_tol;
      options.SetBeta(a * beta * (1 - beta));
      break;
    }
    case 8: {  // Piecewise Map
      // Beta cannot be 0, 0.5 or 1
      if (fabs(beta) < beta_zero_tol) beta += beta_zero_tol;
      if (fabs(beta - 0.5) < beta_zero_tol) beta += beta_zero_tol;
      if (fabs(beta - 1.0) < beta_zero_tol) beta -= beta_zero_tol;
      const double p = 0.24;
      if (beta < p) {
        options.SetBeta(fabs(beta / p));
      }
      else if (beta < 0.5) {
        options.SetBeta(fabs((beta - p) / (0.5 - p)));
      }
      else if (beta < (1.0 - p)) {
        options.SetBeta(fabs((1.0 - p - beta) / (0.5 - p)));
      }
      else {
        options.SetBeta(fabs((1.0 - beta) / p));
      }
      break;
    }
    case 9: {  // Sine Map
      // Beta cannot be 0
      if (fabs(beta) < beta_zero_tol) beta += beta_zero_tol;
      const double a = 2.2;  // Can be 0 < a <= 4
      options.SetBeta((a / 4.0) * sin(my_pi * beta));
      break;
    }
    case 10: {  // Singer Map
      // Beta cannot be 0
      if (fabs(beta) < beta_zero_tol) beta += beta_zero_tol;
      const double mu = 1.0;  // Can be betwen 0.9 and 1.08
      options.SetBeta(mu * (7.86*beta - 23.31*beta*beta + 28.75*beta*beta*beta -
          13.3*beta*beta*beta*beta));
      break;
    }
    case 11: {  // Sinusoidal Map
      // Beta cannot be 0
      if (fabs(beta) < beta_zero_tol) beta += beta_zero_tol;
      const double a = 2.3;
      options.SetBeta(fabs(a * beta * beta * sin(my_pi * beta)));
      break;
    }
    case 12: {  // Tent Map
      // Beta cannot be 0
      if (fabs(beta) < beta_zero_tol) beta += beta_zero_tol;
      const double p = 0.7;
      if (beta < p) {
        options.SetBeta(beta / p);
      }
      else {
        options.SetBeta((10.0 / 3.0) * (1.0 - beta));
      }
      break;
    }
    default: {
      options.SetBeta(beta);
    }
  }
}

std::vector<double> ParticleSwarm::GenerateTrialParticle(std::size_t member)
{
  const auto alpha = options.GetAlpha();
  const auto beta = options.GetBeta();
  const auto one_minus_beta = 1.0 - beta;
  std::normal_distribution<double> norm_rand_num(0.0, 1.0);
  std::vector<double> trial_particle(dimensions_ + 1);
  // The equation for accelerated particle swarm (APSO) is rearranged here:
  // (1-beta)*p + beta*best + alpha*p*rand = p*(1-beta+(alpha*rand)) + beta*best
  for (auto j = 0u; j < dimensions_; ++j) {
    trial_particle[j] = population_[member][j] * (one_minus_beta +
        alpha * norm_rand_num(random_engines_[member])) +
        (beta * best_particle_[j]);
  }
  return trial_particle;
}

void ParticleSwarm::UpdatePopulationMember(GenericCostFunction &CostFunction,
    std::size_t member)
{
  auto trial_particle = GenerateTrialParticle(member);
  if (options.GetUseHardBounds()) bounds.ClampWithinBounds(trial_particle);
  // Only include the new particle if it has a finite cost
  if (CalculateCost(CostFunction, trial_particle)) {
    population_[member] = trial_particle;
  }
}

}  // namespace unfit
