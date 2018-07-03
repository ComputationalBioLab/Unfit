// Unfit: Data fitting and optimization software
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
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
#include "Bounds.hpp"
#include "ConjugateGradient.hpp"
#include "GenericCostFunction.hpp"
#include "LevenbergMarquardt.hpp"
#include "Matrix.hpp"
#include "Options.hpp"

namespace Unfit
{
LevenbergMarquardt::LevenbergMarquardt()
  : jacobian_ {},
    dimensions_ {0},
    observation_size_ {0},
    stored_cost_ {0.0}
{
  options.SetDelta(1e-10);
  is_population_based_ = false;
}

void LevenbergMarquardt::Reset()
{
  ResetGenericOptimizer();
  jacobian_.Assign(0, 0, 0.0);
  dimensions_ = 0;
  observation_size_ = 0;
  stored_cost_ = 0.0;
  options.SetDelta(1e-10);
}

std::vector<double> LevenbergMarquardt::CalculateResiduals(
    GenericCostFunction &cost_function, const std::vector<double> &coordinates)
{
  ++function_evaluations_;
  return cost_function(coordinates);
}

int LevenbergMarquardt::FindMin(GenericCostFunction &cost_function,
    std::vector<double> &coordinates)
{
  // Check the inputs are in order before starting our iterations
  if (coordinates.empty()) return 1;  // No unknown parameters to find
  if (!bounds.IsWithinBounds(coordinates)) return 2;  // Guess out of bounds
  auto residuals = CalculateResiduals(cost_function, coordinates);
  auto cost = SumOfSquares(residuals);
  population_.clear();
  population_.push_back(coordinates);
  population_[0].push_back(cost);
  if (!std::isfinite(cost)) return 6;  // Initial cost is not finite
  dimensions_ = coordinates.size();
  observation_size_ = residuals.size();
  if (observation_size_ < dimensions_) return 3;  // Not enough observations
  const auto cost_tol = options.GetCostTolerance();
  if (cost < cost_tol) {  // Initial guess is close enough
    PrintInitialOutput(cost);
    PrintFinalOutput();
    stored_cost_ = cost;
    return 0;
  }
  PrintInitialOutput(cost);

  const auto epsilon = options.GetEpsilon();
  const auto max_iters = options.GetMaxIterations();
  const auto max_func_evals = options.GetMaxFunctionEvaluations();
  const auto using_broyden_updates = options.GetUseBroydenUpdates();
  const auto recalculate_jacobian_frequency = std::max(static_cast<size_t>(10),
      dimensions_);
  std::vector<double> gradient(observation_size_);
  std::vector<double> hlm(dimensions_);
  Matrix jacobian_sq;
  auto mu = options.GetTau();
  double nu = 2.0;
  bool have_updated_solution = true;
  bool have_new_jacobian = false;
  bool calculate_new_jacobian = false;
  auto num_broyden_updates = recalculate_jacobian_frequency;

  while (iterations_ < max_iters && function_evaluations_ < max_func_evals) {
    ++iterations_;
    // Check for solution breakdown
    if (nu > std::numeric_limits<int>::max()) return 7;

    if (using_broyden_updates) {
      if ((num_broyden_updates == recalculate_jacobian_frequency) ||
          (have_updated_solution && nu > 16.0)) {
        calculate_new_jacobian = true;
        num_broyden_updates = 0;
      }
    }
    else {
      if (have_updated_solution) calculate_new_jacobian = true;
    }

    // Recalculate the jacobian using finite differences
    if (calculate_new_jacobian) {
      FindJacobian(cost_function, residuals, coordinates);
      calculate_new_jacobian = false;
      have_updated_solution = false;
      have_new_jacobian = true;
    }

    // Update JTJ and JTf after a broyden update or jacobian recalculation
    if (have_new_jacobian) {
      have_new_jacobian = false;
      gradient = InnerProduct(jacobian_, residuals);  // JTf
      // If the L_inf norm of gradient is less than epsilon we have converged
      if (std::all_of(begin(gradient), end(gradient), [&](double g) {
          return fabs(g) < epsilon;
          })) break;
      jacobian_sq = InnerProduct(jacobian_, true);  // JTJ
      if (iterations_ == 1) mu *= MaxDiagonalMat(jacobian_sq);  // First iter
      nu = 2.0;
    }

    // Solve JTJ*hlm = JTf to calculate the next step in the solution, hlm
    jacobian_sq.AddConstantToDiagonal(mu);
    ConjugateGradient(hlm, jacobian_sq, Scale(gradient, -1));
    jacobian_sq.AddConstantToDiagonal(-1.0*mu);  // Restore the diagonals
    // If the step size is small then we have converged
    if (FindL2NormOfVector(hlm) <= epsilon*(FindL2NormOfVector(coordinates) +
        epsilon)) break;

    // Create a new solution using the step and calculate the residuals & cost
    auto coordinates_new = LoopUntilWithinBounds(coordinates, hlm);
    auto residuals_new = CalculateResiduals(cost_function, coordinates_new);
    const auto cost_new = SumOfSquares(residuals_new);
    if (!std::isfinite(cost_new)) {  // Cost is not finite, reject the solution
      mu *= nu;
      nu *= 2.0;
      continue;
    }

    if (using_broyden_updates) {
      // If the new cost is an improvement, do a Broyden rank one update
      if ((cost_new < cost) || have_updated_solution) {
        BroydenUpdate(residuals_new, residuals, hlm, jacobian_);
        ++num_broyden_updates;
        have_new_jacobian = true;
      }
    }

    // Check the gain ratio to see if we will accept or reject the new solution
    const auto gain_ratio = (cost - cost_new) / (0.5*mu*InnerProduct(hlm, hlm) -
        0.5*InnerProduct(hlm, gradient));
    if (gain_ratio > 0) {  // New solution is accepted (alt: if cost_new < cost)
      coordinates = std::move(coordinates_new);
      residuals = std::move(residuals_new);
      cost = cost_new;
      if (cost < cost_tol) break;  // Converged if the cost is small enough
      const auto gr_temp = 2.0*gain_ratio - 1.0;
      mu *= std::max(1.0/3.0, 1.0 - (gr_temp*gr_temp*gr_temp));
      nu = 2.0;
      have_updated_solution = true;
    }
    else {  // New solution is rejected
      mu *= nu;
      nu *= 2.0;
    }
    PrintIterationOutput(cost);
  }
  population_[0] = coordinates;
  population_[0].push_back(cost);
  PrintFinalOutput();
  stored_cost_ = cost;
  if (function_evaluations_ >= max_func_evals) return 4;
  if (iterations_ >= max_iters) return 5;
  return 0;
}

void LevenbergMarquardt::FindJacobian(GenericCostFunction &cost_function,
    const std::vector<double> &residuals, const std::vector<double> &variables,
    bool one_sided_difference)
{
  // The passed in residuals vector gives us our residuals at the central point
  // 'x'. Here we just need to calculate the residual at "x+dx" to get an
  // approximation of the Jacobian. Here, we actually create and store J^T
  // instead of J as it makes a lot of calculations access memory in a linear
  // fashion instead of jumping back and forwards.
  //
  //  J^T [#row = dimensions_][#col = observations]
  const auto delta_x = options.GetDelta();
  if (jacobian_.Size() != dimensions_*observation_size_) {
    jacobian_.Assign(dimensions_, observation_size_, 0);
  }
  if (one_sided_difference) {
    auto variables_delta = variables;  // For steps in plus delta_x_
    for (auto i = 0u; i < dimensions_; ++i) {
      variables_delta[i] += delta_x;
      auto residuals_delta = CalculateResiduals(cost_function, variables_delta);
      for (auto j = 0u; j < observation_size_; ++j) {
        jacobian_.values_[i*observation_size_ + j] =
            (residuals_delta[j] - residuals[j]) / delta_x;
      }
      variables_delta[i] = variables[i];
    }
  }
  else {
    auto variables_delta_p = variables;  // For steps in plus delta_x_
    auto variables_delta_m = variables;  // For steps in minus delta_x_
    for (auto i = 0u; i < dimensions_; ++i) {
      variables_delta_p[i] += delta_x;
      variables_delta_m[i] -= delta_x;
      auto residuals_delta_p = CalculateResiduals(cost_function,
          variables_delta_p);
      auto residuals_delta_m = CalculateResiduals(cost_function,
          variables_delta_m);
      for (auto j = 0u; j < observation_size_; ++j) {
        jacobian_.values_[i*observation_size_ + j] =
            (residuals_delta_p[j] - residuals_delta_m[j]) / (2.0 * delta_x);
      }
      variables_delta_p[i] = variables[i];
      variables_delta_m[i] = variables[i];
    }
  }
}

std::vector<double> LevenbergMarquardt::LoopUntilWithinBounds(
    const std::vector<double> &point, std::vector<double> &hlm)
{
  std::vector<double> result {0.0};
  std::vector<double> upper {0.0};
  std::vector<double> lower {0.0};
  bounds.GetBounds(lower, upper);
  while (!bounds.IsWithinBounds(
      result = AddSubtractVectors(point, hlm, true))) {
    for (unsigned i = 0; i < hlm.size(); ++i) {
      if (!bounds.IsWithinBounds(i, result[i])) hlm[i] *= 0.5;
    }
  }
  return result;
//// MB: One possible alternative
//  auto result = AddSubtractVectors(point, hlm, true);
//  bounds.ClampWithinBounds(result);
//  hlm = AddSubtractVectors(result, point, false);
//  return result;
//// MB: Could also modify this to reduce hlm until a valid cost is found
}

void LevenbergMarquardt::BroydenUpdate(const std::vector<double> &residuals_new,
    const std::vector<double> &residuals, const std::vector<double> &hlm,
    Matrix &jacobianT)
{
  // NOTE: This all assumes we have JT (J transpose) passed in, not J
  const auto one_over_hth = 1.0 / std::inner_product(begin(hlm), end(hlm),
      begin(hlm), 0.0);
  // u = f(xnew)
  auto u_vec = residuals_new;
  // u = f(xnew) - f(x)
  std::transform(begin(u_vec), end(u_vec), begin(residuals), begin(u_vec),
      std::minus<double>());
  // u = f(xnew) - f(x) - J*h
  auto it_j = begin(jacobianT.values_);
  for (auto h : hlm) {
    for (auto &u : u_vec) u -= (*it_j++) * h;
  }
  // JT += (u*hT)T, where u = (1/hTh)*(f(x_new)-f(x)-Jh)
  it_j = begin(jacobianT.values_);
  for (auto hT : hlm) {
    for (auto u : u_vec) (*it_j++) += u * hT * one_over_hth;
  }
}

}  // namespace Unfit
