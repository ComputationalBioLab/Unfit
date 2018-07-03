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
#include "Options.hpp"

namespace Unfit
{
Options::Options()
  : max_function_evaluations_ {default_max_function_evaluations_}
  , max_iterations_ {default_max_iterations_}
  , output_level_ {default_output_level_}
  , cost_tolerance_ {default_cost_tolerance_}
  , degenerate_tolerance_ {default_degenerate_tolerance_}
  , geometric_tolerance_ {default_geometric_tolerance_}
  , alpha_ {default_alpha_}
  , beta_ {default_beta_}
  , delta_ {default_delta_}
  , gamma_ {default_gamma_}
  , epsilon_ {default_epsilon_}
  , tau_ {default_tau_}
  , use_adaptive_ {default_use_adaptive_}
  , population_size_ {default_population_size_}
  , seed_ {default_seed_}
  , strategy_ {default_strategy_}
  , weighting_factor_ {default_weighting_factor_}
  , cross_over_ {default_cross_over_}
  , elitism_ {default_elitism_}
  , survival_rate_ {default_survival_rate_}
  , user_has_set_population_size_ {default_user_has_set_population_size_}
  , user_has_set_population_ {default_user_has_set_population_}
  , use_broyden_updates_ {default_use_broyden_updates_}
  , add_initial_to_population_ {default_add_initial_to_population_}
  , use_hard_bounds_ {default_use_hard_bounds_}
  , use_multi_threaded_ {default_use_multi_threaded_}
  , temperature_ {default_temperature_}
  , step_reduction_factor_ {default_step_reduction_factor_}
  , temperature_reduction_factor_ {default_temperature_reduction_factor_}
  , num_cycles_ {default_num_cycles_}
  , num_temperature_loops_ {default_num_temperature_loops_}
{}

void Options::ResetOptions()
{
  max_function_evaluations_ = default_max_function_evaluations_;
  max_iterations_ = default_max_iterations_;
  output_level_ = default_output_level_;
  cost_tolerance_ = default_cost_tolerance_;
  degenerate_tolerance_ = default_degenerate_tolerance_;
  geometric_tolerance_ = default_geometric_tolerance_;
  alpha_ = default_alpha_;
  beta_ = default_beta_;
  delta_ = default_delta_;
  gamma_ = default_gamma_;
  epsilon_ = default_epsilon_;
  tau_ = default_tau_;
  use_adaptive_ = default_use_adaptive_;
  population_size_ = default_population_size_;
  seed_ = default_seed_;
  strategy_ = default_strategy_;
  weighting_factor_ = default_weighting_factor_;
  cross_over_ = default_cross_over_;
  elitism_ = default_elitism_;
  survival_rate_ = default_survival_rate_;
  user_has_set_population_size_ = default_user_has_set_population_size_;
  user_has_set_population_ = default_user_has_set_population_;
  use_broyden_updates_ = default_use_broyden_updates_;
  add_initial_to_population_ = default_add_initial_to_population_;
  use_hard_bounds_ = default_use_hard_bounds_;
  use_multi_threaded_ = default_use_multi_threaded_;
  temperature_ = default_temperature_;
  step_reduction_factor_ = default_step_reduction_factor_;
  temperature_reduction_factor_ = default_temperature_reduction_factor_;
  num_cycles_ = default_num_cycles_;
  num_temperature_loops_ = default_num_temperature_loops_;
}

double Options::GetCostTolerance() const noexcept
{
  return cost_tolerance_;
}

void Options::SetCostTolerance(double tolerance)
{
  if (tolerance <= 0.0) return;
  cost_tolerance_ = tolerance;
}

double Options::GetDegenerateTolerance() const noexcept
{
  return degenerate_tolerance_;
}

void Options::SetDegenerateTolerance(double tolerance)
{
  if (tolerance <= 0.0) return;
  degenerate_tolerance_ = tolerance;
}

double Options::GetGeometricTolerance() const noexcept
{
  return geometric_tolerance_;
}

void Options::SetGeometricTolerance(double tolerance)
{
  if (tolerance <= 0.0) return;
  geometric_tolerance_ = tolerance;
}

double Options::GetAlpha() const noexcept
{
  return alpha_;
}

void Options::SetAlpha(double alpha)
{
  alpha_ = alpha;
}

double Options::GetBeta() const noexcept
{
  return beta_;
}

void Options::SetBeta(double beta)
{
  beta_ = beta;
}

double Options::GetDelta() const noexcept
{
  return delta_;
}

void Options::SetDelta(double delta)
{
  delta_ = delta;
}

double Options::GetGamma() const noexcept
{
  return gamma_;
}

void Options::SetGamma(double gamma)
{
  gamma_ = gamma;
}

double Options::GetEpsilon() const noexcept
{
  return epsilon_;
}

void Options::SetEpsilon(double epsilon)
{
  epsilon_ = epsilon;
}

void Options::SetTau(double tau)
{
  if (tau < 0.0) return;
  tau_ = tau;
}

double Options::GetTau() const noexcept
{
  return tau_;
}

unsigned Options::GetMaxFunctionEvaluations() const noexcept
{
  return max_function_evaluations_;
}

void Options::SetMaxFunctionEvaluations(unsigned max_func_evals)
{
  max_function_evaluations_ = max_func_evals;
}

unsigned Options::GetMaxIterations() const noexcept
{
  return max_iterations_;
}

void Options::SetMaxIterations(unsigned max_iters)
{
  max_iterations_ = max_iters;
}

void Options::GetNelderMeadStepSizes(double &alpha, double &beta, double &delta,
    double &gamma)
{
  alpha = alpha_;
  beta = beta_;
  delta = delta_;
  gamma = gamma_;
}

void Options::SetNelderMeadStepSizes(double alpha, double beta, double delta,
    double gamma)
{
  if (gamma <= 0.0 || gamma >= 1.0) return;
  if (delta <= 0.0 || delta >= 1.0) return;
  if (gamma >= alpha) return;
  if (alpha >= beta) return;
  alpha_ = alpha;
  beta_ = beta;
  delta_ = delta;
  gamma_ = gamma;
}

unsigned Options::GetOutputLevel() const noexcept
{
  return output_level_;
}

void Options::SetOutputLevel(unsigned output_level)
{
  output_level_ = output_level;
}

bool Options::GetUseAdaptiveParameters() const noexcept
{
  return use_adaptive_;
}

void Options::SetUseAdaptiveParameters(bool adaptive)
{
  use_adaptive_ = adaptive;
}

unsigned Options::GetPopulationSize() const noexcept
{
  return population_size_;
}

void Options::SetPopulationSize(unsigned pop_size)
{
  user_has_set_population_size_ = true;
  population_size_ = pop_size;
}

unsigned Options::GetRandomSeed() const noexcept
{
  return seed_;
}

void Options::SetRandomSeed(unsigned seed)
{
  seed_ = seed;
}

unsigned Options::GetStrategy() const noexcept
{
  return strategy_;
}

void Options::SetStrategy(unsigned strategy)
{
  strategy_ = strategy;
}

double Options::GetCrossOver() const noexcept
{
  return cross_over_;
}

void Options::SetCrossOver(double cross_over)
{
  cross_over_ = cross_over;
  if (cross_over < 0.0) cross_over_ = 0.0;
  if (cross_over > 1.0) cross_over_ = 1.0;
}

double Options::GetWeightingFactor() const noexcept
{
  return weighting_factor_;
}

void Options::SetWeightingFactor(double weighting_factor)
{
  weighting_factor_ = weighting_factor;
}

bool Options::GetUserSetPopulationSize() const noexcept
{
  return user_has_set_population_size_;
}

unsigned Options::GetElitism() const noexcept
{
  return elitism_;
}

void Options::SetElitism(unsigned elite)
{
  if (elite <= GetPopulationSize()) elitism_ = elite;
}

double Options::GetSurvivalRate() const noexcept
{
  return survival_rate_;
}

void Options::SetSurvivalRate(double rate)
{
  survival_rate_ = rate;
  if (rate < 0.0) survival_rate_ = 0.0;
  if (rate > 1.0) survival_rate_ = 1.0;
}

bool Options::GetUseBroydenUpdates() const noexcept
{
  return use_broyden_updates_;
}

void Options::SetUseBroydenUpdates(bool use_broyden)
{
  use_broyden_updates_ = use_broyden;
}

bool Options::GetAddInitialToPopulation() const noexcept
{
  return add_initial_to_population_;
}

void Options::SetAddInitialToPopulation(bool add_initial)
{
  add_initial_to_population_ = add_initial;
}

bool Options::GetUseHardBounds() const noexcept
{
  return use_hard_bounds_;
}

void Options::SetUseHardBounds(bool use_hard_bounds)
{
  use_hard_bounds_ = use_hard_bounds;
}

bool Options::GetUserSetPopulation() const noexcept
{
  return user_has_set_population_;
}

void Options::SetUserSetPopulation(bool has_set_population)
{
  user_has_set_population_ = has_set_population;
}

bool Options::GetUseMultiThreaded() const noexcept
{
  return use_multi_threaded_;
}

void Options::SetUseMultiThreaded(bool use_multi_threaded)
{
  use_multi_threaded_ = use_multi_threaded;
}

double Options::GetTemperature() const noexcept
{
  return temperature_;
}

void Options::SetTemperature(double temperature)
{
  temperature_ = temperature;
}

double Options::GetStepReductionFactor() const noexcept
{
  return step_reduction_factor_;
}

void Options::SetStepReductionFactor(double step_factor)
{
  step_reduction_factor_ = step_factor;
}

double Options::GetTemperatureReductionFactor() const noexcept
{
  return temperature_reduction_factor_;
}

void Options::SetTemperatureReductionFactor(double temperature_factor)
{
  temperature_reduction_factor_ = temperature_factor;
}

int Options::GetNumberOfCycles() const noexcept
{
  return num_cycles_;
}

void Options::SetNumberOfCycles(int num_cycles)
{
  num_cycles_ = num_cycles;
}

int Options::GetNumberOfTemperatureLoops() const noexcept
{
  return num_temperature_loops_;
}

void Options::SetNumberOfTemperatureLoops(int num_temperature_loops)
{
  num_temperature_loops_ = num_temperature_loops;
}
}  // namespace Unfit
