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
#ifndef UNFIT_INCLUDE_OPTIONS_HPP_
#define UNFIT_INCLUDE_OPTIONS_HPP_

namespace Unfit
{
/**
 * This class is designed to handle all of the optimization options for Unfit.
 * The model is that all options should be here, so as to provide a uniform
 * interface to setting an option from any optimization method. The drawback is
 * that some options may not be valid for some methods, but if someone sets an
 * option that is not relevant, it will have no effect. In the future it may be
 * a good idea to warn the user if an option will do nothing for their chosen
 * optimization algorithm.
 */
class Options
{
 public:
  /**
   * Set all of the optimization options to their default values
   */
  Options();

  /**
   * Resets all of the optimization options to their default values
   */
  void ResetOptions();

  /**
   * Get the value of the cost tolerance which controls the termination of
   * the optimization based on the value of the user-supplied cost function.
   *
   * \return The value of cost tolerance
   */
  double GetCostTolerance() const noexcept;

  /**
   * Set the value of the cost tolerance which controls the termination of
   * the optimization based on the value of the user-supplied cost function. The
   * tolerance must be positive, and anything else will be ignored.
   *
   * \param tolerance The value of cost tolerance
   */
  void SetCostTolerance(double tolerance);

  /**
   * Get the value of the degenerate tolerance. At the moment this is only used
   * in the Nelder-Mead algorithm to check the size of the simplex. If all of
   * the points become colinear then we need to restart with a new simplex.
   * This parameter controls how colinear the vertices can be before a restart
   * is required.
   *
   * \return The value of degenerate tolerance
   */
  double GetDegenerateTolerance() const noexcept;

  /**
   * Set the value of the degenerate tolerance. At the moment this is only used
   * in the Nelder-Mead algorithm to check the size of the simplex. If all of
   * the points become colinear then we need to restart with a new simplex.
   * This parameter controls how colinear the vertices can be before a restart
   * is required. The tolerance must be positive, and anything else will be
   * ignored.
   *
   * \param tolerance The value of degenerate tolerance
   */
  void SetDegenerateTolerance(double tolerance);

  /**
   * Get the value of the geometric tolerance which controls the termination of
   * the optimization based on a metric related to the size of the problem
   * space. It will be algorithm specific as to whether this is used and as to
   * what it actually checks.
   *
   * \return The value of geometric tolerance
   */
  double GetGeometricTolerance() const noexcept;

  /**
   * Set the value of the geometric tolerance which controls the termination of
   * the optimization based on a metric related to the size of the problem
   * space. It will be algorithm specific as to whether this is used and as to
   * what it actually checks. The tolerance must be positive, and anything else
   * will be ignored.
   *
   * \param tolerance The value of geometric tolerance
   */
  void SetGeometricTolerance(double tolerance);

  /**
   * Get the value of the alpha parameter. It will be algorithm specific as to
   * whether this is used and as to what it actually checks. In the NelderMead
   * algorithm this scales the reflect operation.
   *
   * \return The value of alpha
   */
  double GetAlpha() const noexcept;

  /**
   * Set the value of the alpha parameter. It will be algorithm specific as to
   * whether this is used and as to what it actually checks. In the NelderMead
   * algorithm this scales the reflect operation. If you are using NelderMead
   * you should use SetNelderMeadStepSizes as this does not check the input
   * is in keeping with the NelderMead constraints.
   *
   * \param alpha The value of alpha
   */
  void SetAlpha(double alpha);

  /**
   * Get the value of the beta parameter. It will be algorithm specific as to
   * whether this is used and as to what it actually checks. In the NelderMead
   * algorithm this scales the expand operation.
   *
   * \return The value of beta
   */
  double GetBeta() const noexcept;

  /**
   * Set the value of the beta parameter. It will be algorithm specific as to
   * whether this is used and as to what it actually checks. In the NelderMead
   * algorithm this scales the expand operation. If you are using NelderMead
   * you should use SetNelderMeadStepSizes as this does not check the input
   * is in keeping with the NelderMead constraints.
   *
   * \param beta The value of beta
   */
  void SetBeta(double beta);

  /**
   * Get the value of the delta parameter. It will be algorithm specific as to
   * whether this is used and as to what it actually checks. In the NelderMead
   * algorithm this scales the shrink operation. This is the
   * size of the finite difference step used to approximate the jacobian
   * matrix in the LevenbergMarquardt algorithm.
   *
   * \return The value of delta
   */
  double GetDelta() const noexcept;

  /**
   * Set the value of the delta parameter. It will be algorithm specific as to
   * whether this is used and as to what it actually checks. In the NelderMead
   * algorithm this scales the shrink operation. If you are using NelderMead
   * you should use SetNelderMeadStepSizes as this does not check the input
   * is in keeping with the NelderMead constraints. This is the
   * size of the finite difference step used to approximate the jacobian
   * matrix in the LevenbergMarquardt algorithm.
   *
   * \param delta The value of delta
   */
  void SetDelta(double delta);

  /**
   * Get the value of the gamma parameter. It will be algorithm specific as to
   * whether this is used and as to what it actually checks. In the NelderMead
   * algorithm this scales the contract operations.
   *
   * \return The value of gamma
   */
  double GetGamma() const noexcept;

  /**
   * Set the value of the gamma parameter. It will be algorithm specific as to
   * whether this is used and as to what it actually checks. In the NelderMead
   * algorithm this scales the contract operations. If you are using NelderMead
   * you should use SetNelderMeadStepSizes as this does not check the input
   * is in keeping with the NelderMead constraints.
   *
   * \param gamma The value of gamma
   */
  void SetGamma(double gamma);

  /**
   * This method returns the value of epsilon, which is used in the Levenberg
   * Marquardt algorithm as a convergence measure for the gradients at each
   * iteration. If none of the gradients are larger than epsilon we have
   * convergence.
   *
   * \return the value of epsilon
   */
  double GetEpsilon() const noexcept;

  /**
   * This method sets the value of epsilon, which is used in the Levenberg
   * Marquardt algorithm as a convergence measure for the gradients at each
   * iteration. If none of the gradients are larger than epsilon we have
   * convergence.
   *
   * \param epsilon the intended value of epsilon
   */
  void SetEpsilon(double epsilon);

  /**
   * Set the value of the tau parameter. It will be algorithm specific as to
   * whether this is used and as to what it actually checks. In the Levenberg
   * Marquardt algorithm tau scales the initial damping parameter mu. Here tau
   * must be positive and any negative input will be ignored.
   *
   * \param tau The value of tau
   */
  void SetTau(double tau);

  /**
   * Set the value of the tau parameter. It will be algorithm specific as to
   * whether this is used and as to what it actually checks. In the Levenberg
   * Marquardt algorithm tau scales the initial damping parameter mu.
   *
   * \return The value of tau
   */
  double GetTau() const noexcept;

  /**
   * Get the maximum number of function evaluations (calls to the cost function)
   * for the optimization problem. The optimization will terminate if the
   * maximum number of function evaluations is reached before convergence.
   *
   * \return Maximum number of function evaluations
   */
  unsigned GetMaxFunctionEvaluations() const noexcept;

  /**
   * Set the maximum number of function evaluations (calls to the cost function)
   * for the optimization problem. The optimization will terminate if the
   * maximum number of function evaluations is reached before convergence.
   *
   * \param max_func_evals Maximum number of function evaluations
   */
  void SetMaxFunctionEvaluations(unsigned max_func_evals);

  /**
   * Get the maximum number of iterations for the optimization problem. The
   * optimization will terminate if the maximum number of iterations is reached
   * before convergence.
   *
   * \return Maximum number of iterations
   */
  unsigned GetMaxIterations() const noexcept;

  /**
   * Set the maximum number of iterations for the optimization problem. The
   * optimization will terminate if the maximum number of iterations is reached
   * before convergence.
   *
   * \param max_iters Maximum number of iterations
   */
  void SetMaxIterations(unsigned max_iters);

  /**
   * Get the step sizes for the NelderMead algorithm. The alpha parameter scales
   * the reflect operation, the beta parameter scales the expand operation, the
   * delta parameter scales the shrink operation and the gamma parameter scales
   * the contract operations. The design of the interface is such that you can
   * get the parameters, change one or more, then pass them back in to set them
   * without having to expose the user to the added complexity of something like
   * a tuple.
   *
   * \param alpha Scales the reflect operation (default 1.0)
   * \param beta  Scales the expand operation (default 2.0)
   * \param delta Scales the shrink operation (default 0.5)
   * \param gamma Scales the contract operation (default 0.5)
   */
  void GetNelderMeadStepSizes(double &alpha, double &beta, double &delta,
    double &gamma);

  /**
   * Set the step sizes for the NelderMead algorithm. If you are not using the
   * NelderMead algorithm then setting these options will have no effect. The
   * alpha parameter scales the reflect operation, the beta parameter scales
   * the expand operation, the delta parameter scales the shrink operation and
   * the gamma parameter scales the contract operations. To ensure normal
   * operations, it is enforced that 0 < gamma, delta < 1, and that gamma <
   * alpha < beta. If this is not the case then the set operation will fail and
   * the parameters will remain unchanged.
   *
   * \param alpha Scales the reflect operation (default 1.0)
   * \param beta  Scales the expand operation (default 2.0)
   * \param delta Scales the shrink operation (default 0.5)
   * \param gamma Scales the contract operation (default 0.5)
   */
  void SetNelderMeadStepSizes(double alpha, double beta, double delta,
    double gamma);

  /**
   * Gets the current level of output. This may be algorithm specific, but in
   * general 0 = no output; 1 = iteration counter only; 2 = iteration by
   * iteration information; 3+ = same as 2, plus additional output (e.g. final
   * result).
   *
   * \return The current level of output
   */
  unsigned GetOutputLevel() const noexcept;

  /**
   * Set the desired level of output. This may be algorithm specific, but in
   * general 0 = no output; 1 = iteration counter only; 2 = iteration by
   * iteration information; 3+ = same as 2, plus additional output (e.g. final
   * result).
   *
   * \param output_level The selected level of output
   */
  void SetOutputLevel(unsigned output_level);

  /**
   * Get a boolean which states whether or not the algorithm is currently using
   * adaptive parameters.
   *
   * \return True if adaptive parameters are used, otherwise false
   */
  bool GetUseAdaptiveParameters() const noexcept;

  /**
   * Set the maximum number of iterations for the optimization problem. The
   * optimization will terminate if the maximum number of iterations is reached
   * before convergence.
   *
   * \param adaptive Set to true if adaptive parameters are wanted, or turn
   *   them off by passing false
   */
  void SetUseAdaptiveParameters(bool adaptive);

  /**
   * Returns the size of the population used in the GeneticAlgorithm or
   * DifferentialEvolution solution method.
   *
   * \return The size of the population
   */
  unsigned GetPopulationSize() const noexcept;

  /**
   * Set the size of the population used by GeneticAlgorithm and
   * DifferentialEvolution. The minimum population size is three for GA and six
   * for DE, and if you try to set a number less than this these minimums will
   * be imposed. The default is to use 10x the number of unknowns.
   *
   * \param pop_size The intended size of the population
   */
  void SetPopulationSize(unsigned pop_size);

  /**
   * For both GeneticAlgorithm and DifferentialEvolution, the default population
   * size is determined by the code. If a user wants to override this then we
   * need to set a flag. Call this method to determine if the user has set their
   * own population size.
   *
   * \return true if the user has specified a population size, otherwise false
   */
  bool GetUserSetPopulationSize() const noexcept;

  /**
   * Gets the seed to be used by the random number engine used by
   * Genetic Algorithm and Differential Evolution.
   *
   * \return The seed for the random number generator
   */
  unsigned GetRandomSeed() const noexcept;

  /**
   * Sets the seed to be used by the random number engine used by
   * Genetic Algorithm and Differential Evolution. Note that value of 0 and 1
   * generates the same set of values.
   *
   * \param seed The seed for the random number generator
   */
  void SetRandomSeed(unsigned seed);

  /**
   * Get the current strategy used by the differential evolution algorithm
   * (range = 1-10).
   *
   * \return the strategy to adopt for differential evolution
   */
  unsigned GetStrategy() const noexcept;

  /**
   * Choose amongst the 10 strategies (numbered 1-10) for the vector
   * combinations at each iteration for the differential evolution algorithm.
   * Different strategies can work better or worse for different problems. The
   * default is 1, which seems to work for our test problems.
   *
   * \param strategy the strategy to adopt for differential evolution
   */
  void SetStrategy(unsigned strategy);

  /**
   * Get the cross over probability (CR) used by the differential evolution
   * algorithm
   *
   * \return the cross over probability
   */
  double GetCrossOver() const noexcept;

  /**
   * Set the cross over probability (CR) to be used by the differential
   * evolution algorithm. This determines how likely it is that a population
   * member will evolve.
   *
   * \param cross_over the cross over probability
   */
  void SetCrossOver(double cross_over);

  /**
   * Get the differential weighting factor (F) used by the differential
   * evolution algorithm
   *
   * \return the differential weighting factor
   */
  double GetWeightingFactor() const noexcept;

  /**
   * Set the differential weighting factor (F) to be used by the differential
   * evolution algorithm. This is used as a weight when combining different
   * population vectors.
   *
   * \param weighting_factor the differential weighting factor
   */
  void SetWeightingFactor(double weighting_factor);

  /**
   * Returns whether or not the initial coordinates should form part of the
   * initial population for Differential Evolution. The default is no (false).
   *
   * \return true if the passed in coordinates should be added to the population
   */
  bool GetAddInitialToPopulation() const noexcept;

  /**
   * Choose whether or not the initial coordinates should form part of the
   * initial population for Differential Evolution. The default is no (false).
   *
   * \param add_initial true if the passed in coordinates should be included in
   *        the population
   */
  void SetAddInitialToPopulation(bool add_initial);

  /**
   * Returns whether or not we are using Broyden rank one updates for the
   * Levenberg Marquart algorithm. The default is yes (true).
   *
   * \return true if using Broyden rank one updates, otherwise false
   */
  bool GetUseBroydenUpdates() const noexcept;

  /**
   * Choose whether or not we are using Broyden rank one updates for the
   * Levenberg Marquart algorithm. The default is yes (true).
   *
   * \param use_broyden true if using Broyden rank one updates, otherwise false
   */
  void SetUseBroydenUpdates(bool use_broyden);

  /**
   * Elite chromosomes are immune from the mutation process associated with the
   * Genetic Algorithm optimisation technique. This returns the number of elite
   * chromosomes.
   *
   * \return The number of elite chromosomes
   */
  unsigned GetElitism() const noexcept;

  /**
   * Elite chromosomes are immune from the mutation process associated with the
   * Genetic Algorithm optimisation technique. This sets the number of elite
   * chromosomes.
   *
   * \param elite The number of elite chromosomes
   */
  void SetElitism(unsigned elite);

  /**
   * For Genetic Algorithms, in each generation, part of the population survives
   * and part is replaced. This method gets the proportion of chromosomes that
   * survive at each generation.
   *
   * \return The proportion of the population that survives in each generation
   */
  double GetSurvivalRate() const noexcept;

  /**
   * For Genetic Algorithms, in each generation, part of the population survives
   * and part is replaced. This method sets the proportion of chromosomes that
   * survive at each generation. This should be between zero and one, with the
   * caveat that a minimum of two chromosomes must survive in any given
   * generation, to produce offspring for the next. Therefore, if you set the
   * value to zero, two will still survive. If you set it to one, then the only
   * change between one generation and the next will be the mutation component.
   *
   * \param rate The fraction of the population that survives in a generation
   */
  void SetSurvivalRate(double rate);

  /**
   * For Differential Evolution, this returns whether or not the specified
   * bounds are strictly enforced. If false, the solution can be outside the
   * specified bounds.
   *
   * \return Whether or not hard bounds are in use
   */
  bool GetUseHardBounds() const noexcept;

  /**
   * For Differential Evolution, this sets whether or not the specified
   * bounds are strictly enforced. If false, the solution can be outside the
   * specified bounds (bounds are then essentially only used to generate the
   * initial populations). If true, the solution must lie within the same bounds
   * that were used to generate the population.
   *
   * \param use_hard_bounds Whether or not hard bounds are in use
   */
  void SetUseHardBounds(bool use_hard_bounds);

  /**
   * Returns whether or not the user has provided a population to the optimizer.
   *
   * \return True if a population has been given, false otherwise
   */
  bool GetUserSetPopulation() const noexcept;

  /**
   * Sets whether or not the user has provided a population to the optimizer.
   *
   * \param has_set_population True if a population is given, false otherwise
   */
  void SetUserSetPopulation(bool has_set_population);

  /**
   * Returns whether or not the algorithms in Unfit run multi-threaded or on
   * a single thread (parallel or serial execution).
   *
   * \return True if using multi-threaded code, false otherwise
   */
  bool GetUseMultiThreaded() const noexcept;

  /**
   * Sets whether or not the algorithms in Unfit run multi-threaded or on
   * a single thread (parallel or serial execution).
   *
   * \param use_multi_threaded True if using multi-threaded code, false
   *        otherwise
   */
  void SetUseMultiThreaded(bool use_multi_threaded);

  /**
   * Returns the initial temperature (in Kelvin) to be used in the first
   * iteration for Simulated Annealing algorithm. The default initial
   * temperature is 1000.
   *
   * \return The initial annealing temperature
   */
  double GetTemperature() const noexcept;

  /**
   * Sets the initial temperature (in Kelvin) to be used in the first iteration
   * for Simulated Annealing algorithm.
   *
   * \param temperature The initial temperature (>0)
   */
  void SetTemperature(double temperature);

  /**
   * Returns the step reduction factor to be used to reduce the step size
   * after each temperature loop in Simulated Annealing. The default value
   * is 0.9
   *
   * \return The step reduction factor
   */
  double GetStepReductionFactor() const noexcept;

  /**
   * Sets the step reduction factor to be used to reduce the step size
   * after each temperature loop in Simulated Annealing.
   *
   * \param step_factor The factor to reduce step size
   */
  void SetStepReductionFactor(double step_factor);

  /**
   * Returns the temperature reduction factor to be used to reduce temperature
   * after each temperature loop in Simulated Annealing. The default value
   * is 0.5
   *
   * \return The temperature reduction factor
   */
  double GetTemperatureReductionFactor() const noexcept;

  /**
   * Sets the temperature reduction factor to be used to reduce temperature
   * after each temperature loop in Simulated Annealing.
   *
   * \param temperature_factor The factor to reduce temperature
   */
  void SetTemperatureReductionFactor(double temperature_factor);

  /**
   * Returns the number of cycles to be used in Simulated Annealing.
   * The default value is 20
   *
   * \return The number of cycles
   */
  int GetNumberOfCycles() const noexcept;

  /**
   * Sets the number of cycles to be used in Simulated Annealing.
   *
   * \param num_cycles The number of cycles
   */
  void SetNumberOfCycles(int num_cycles);

  /**
   * Returns the number of temperature loops to be used in Simulated
   * Annealing. The default value is 5
   *
   * \return The number of temperature loops
   */
  int GetNumberOfTemperatureLoops() const noexcept;

  /**
   * Sets the number of temperature loops to be used in Simulated Annealing.
   *
   * \param num_temperature_loops The number of temperature loops
   */
  void SetNumberOfTemperatureLoops(int num_temperature_loops);

 private:
  /** Maximum number of function evaluations permitted before termination*/
  unsigned max_function_evaluations_;
  /** Maximum number of iterations permitted before termination*/
  unsigned max_iterations_;
  /** Specify the amount of output to write out during the optimization*/
  unsigned output_level_;
  /**
   * The tolerance for the convergence of the cost function. This represents
   * the epsilon parameter in the Levenberg-Marquardt algorithm.
   */
  double cost_tolerance_;
  /** A tolerance for Nelder-Mead to check the colinearity of the vertices*/
  double degenerate_tolerance_;
  /** A tolerance for the convergence of the solution space*/
  double geometric_tolerance_;
  /**
   * Used as a scale for the reflect operation for the NelderMead algorithm.
   * Used as initial guess for the step size in the steepest descent algorithm
   */
  double alpha_;
  /**
   * Used as a scale for the expand operation in the NelderMead algorithm.
   * Used as Armijo rule factor in the steepest descent algorithm.
   */
  double beta_;
  /**
   * Scales the shrink operation for the NelderMead algorithm. This is the
   * size of the finite difference step used to approximate the jacobian
   * matrix in the LevenbergMarquardt algorithm.
   */
  double delta_;
  /** Scales the contract operation for the NelderMead algorithm*/
  double gamma_;
  /** Sets the convergence tolerance for the increment vector in LM*/
  double epsilon_;
  /** Scales the initial damping parameter (mu) in LevenbergMarquardt*/
  double tau_;
  /** Sets whether or not adaptive parameters are to be used*/
  bool use_adaptive_;
  /** Sets the size of the population for the GA and DE algorithms*/
  unsigned population_size_;
  /** Sets the seed for the random number generators in GA and DE*/
  unsigned seed_;
  /** Sets the strategy for the DE algorithm (1-10)*/
  unsigned strategy_;
  /** Sets the weighting factor for the DE algorithm (sometimes called F)*/
  double weighting_factor_;
  /** Sets the cross over for the DE algorithm (sometimes called CR)*/
  double cross_over_;
  /** Sets the number of chromosomes immune from mutation*/
  unsigned elitism_;
  /** Sets the proportion of chromosomes that survive each generation*/
  double survival_rate_;
  /** True if the user has set the population size themselves, otherwise false*/
  bool user_has_set_population_size_;
  /** True if the user has set the population themselves, otherwise false*/
  bool user_has_set_population_;
  /**
   * If true, the Levenberg Marquardt algorithm will use Broyden rank one
   * updates to minimise cost-function calls. If false gradients are always
   * calculated using finite differences.
   */
  bool use_broyden_updates_;
  /**
   * If true, the initial coordinates passed in to FindMin are included in the
   * initial population. If false, all population members are randomly
   * generated. (Currently only implemented for Differential Evolution)
   */
  bool add_initial_to_population_;
  /**
   * Whether or not the solution can lie outside the specified bounds.
   * Currently only implemented for Differential Evolution.
   */
  bool use_hard_bounds_;
  /**
   * Whether or not to use a multi-threaded implementation of the algorithms,
   * where available.
   */
  bool use_multi_threaded_;
  /** Sets the initial temperature for the SA algorithm*/
  double temperature_;
  /** Sets the step reduction factor for the SA algorithm*/
  double step_reduction_factor_;
  /** Sets the temperature reduction factor for the SA algorithm*/
  double temperature_reduction_factor_;
  /** Sets the number of cycles for the SA algorithm*/
  int num_cycles_;
  /** Sets the number of temperature loops for the SA algorithm*/
  int num_temperature_loops_;

  /** Default maximum number of function evaluations before termination*/
  static constexpr unsigned default_max_function_evaluations_ = 100000;
  /** Default maximum number of iterations before termination*/
  static constexpr unsigned default_max_iterations_ = 10000;
  /** Default output level for the optimization*/
  static constexpr unsigned default_output_level_ = 0;
  /** Default tolerance for the convergence of the cost function*/
  static constexpr double default_cost_tolerance_ = 1e-12;
  /** Default tolerance for Nelder-Mead to check vertex colinearity*/
  static constexpr double default_degenerate_tolerance_ = 1e-8;
  /** Default tolerance for the convergence of the solution space*/
  static constexpr double default_geometric_tolerance_ = 1e-4;
  /** Default value of the parameter alpha*/
  static constexpr double default_alpha_ = 1.0;
  /** Default value of the parameter beta*/
  static constexpr double default_beta_ = 2.0;
  /** Default scale the shrink operation for the NelderMead algorithm*/
  static constexpr double default_delta_ = 0.5;
  /** Default convergence tolerance for the increment vector in LM*/
  static constexpr double default_epsilon_ = 1e-12;
  /** Default scale the contract operation for the NelderMead algorithm*/
  static constexpr double default_gamma_ = 0.5;
  /** Default scale the initial damping parameter (mu) in LevenbergMarquardt*/
  static constexpr double default_tau_ = 1e-3;
  /** Default is to not use adaptive parameters*/
  static constexpr bool default_use_adaptive_ = false;
  /** Default size of the population for the GA and DE algorithms*/
  static constexpr unsigned default_population_size_ = 20;
  /** Default seed for the random number generators in GA and DE*/
  static constexpr unsigned default_seed_ = 0;
  /** Default strategy for the DE algorithm*/
  static constexpr unsigned default_strategy_ = 1;
  /** Default weighting factor for the DE algorithm (sometimes called F)*/
  static constexpr double default_weighting_factor_ = 0.8;
  /** Default cross over for the DE algorithm (sometimes called CR)*/
  static constexpr double default_cross_over_ = 0.9;
  /** Sets the number of chromosomes immune from mutation*/
  static constexpr unsigned default_elitism_ = 1;
  /** Sets the proportion of chromosomes that survive each generation*/
  static constexpr double default_survival_rate_ = 0.5;
  /** By default, the program will try to choose a population size*/
  static constexpr bool default_user_has_set_population_size_ = false;
  /** By default, the program will generate a population*/
  static constexpr bool default_user_has_set_population_ = false;
  /** By default, the program will use Broyden rank one updates for LM*/
  static constexpr bool default_use_broyden_updates_ = true;
  /** By default, the program will generate the whole population for DE*/
  static constexpr bool default_add_initial_to_population_ = false;
  /** By default, the program allows solutions outside the set bounds for DE*/
  static constexpr bool default_use_hard_bounds_ = false;
  /** By default, the program will run everything on one thread*/
  static constexpr bool default_use_multi_threaded_ = false;
  /** Default initial temperature for the SA algorithm*/
  static constexpr double default_temperature_ = 1000.0;
  /** Default step reduction factor for the SA algorithm*/
  static constexpr double default_step_reduction_factor_ = 0.9;
  /** Default temperature reduction factor for the SA algorithm*/
  static constexpr double default_temperature_reduction_factor_ = 0.5;
  /** Default number of cycles for the SA algorithm*/
  static constexpr int default_num_cycles_ = 20;
  /** Default number of temperature loops for the SA algorithm*/
  static constexpr int default_num_temperature_loops_ = 5;
};

}  // namespace Unfit

#endif
