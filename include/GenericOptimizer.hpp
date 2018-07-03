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
#ifndef UNFIT_INCLUDE_GENERICOPTIMIZER_HPP_
#define UNFIT_INCLUDE_GENERICOPTIMIZER_HPP_

#include <atomic>
#include <random>
#include <vector>
#include "Bounds.hpp"
#include "Options.hpp"

namespace Unfit
{
class GenericCostFunction;  // Forward declaration, see GenericCostFunction.hpp

/**
 * A generic base class for optimization methods to derive from. This defines
 * the FindMin interface through which the optimization is performed, but is not
 * a pure interface class as it may also contains common features for the
 * optimization algorithms.
 */
class GenericOptimizer
{
 public:
  /**
   * Default constructor. Here we override the generated one as we need to
   * initialize the member variables to something sensible.
   */
  GenericOptimizer();

  /**
   * As we are deriving from this class, the destructor should be virtual. In
   * this case (at present) we have nothing that will not be deleted when the
   * class goes out of scope so an empty destructor method is fine.
   */
  virtual ~GenericOptimizer();

  /**
   * All derived classes must implement the FindMin method. This is the main
   * method that should be called to perform an optimization.
   *
   * \param cost_function The cost function to be minimized
   * \param coordinates On entry this is the initial guess of the coordinates
   *        that need fitting. On exit this is overwritten by the final result
   *        from the fitting procedure.
   * \return An integer error code, which should be zero to indicate success
   */
  virtual int FindMin(GenericCostFunction &cost_function,
      std::vector<double> &coordinates) = 0;

  /**
   * Each optimization class must implement a Reset method to set all of their
   * parameters back to their default values.
   */
  virtual void Reset() = 0;

  /**
   * This method resets all of the GenericOptimizer parameters back to their
   * default values and should therefore be called from the Reset method of all
   * optimizers deriving from this class.
   */
  void ResetGenericOptimizer();

  /**
   * This method which returns the cost of a particular member of the current
   * solution. For optimizers where more than one cost is stored the optional
   * parameter can be used to choose which cost to get. The costs are not
   * guaranteed to be sorted. If an invalid member index is passed in as an
   * argument the returned cost will be infinity.
   *
   * DifferentialEvolution: Members are unsorted \n
   * GeneticAlgorithm: Members are sorted (best = 0, worst = n) \n
   * LevenbergMarquardt: Only one cost is stored (index = 0) \n
   * NelderMead: Members are sorted (best = 0, worst = n)
   *
   * \param index The index of the member whose cost will be returned
   * \return The cost of the solution of the specified member
   */
  virtual double GetCost(std::size_t index = 0) const noexcept;

  /**
   * This method returns whether or not the current optimization method is
   * population based (e.g. Differential Evolution) or not (e.g Nelder-Mead).
   *
   * \return True if the optimizer is population based, otherwise false
   */
  virtual bool GetIsPopulationBased() const noexcept;

  /**
   * This method returns the number of iterations the optimizer has gone through
   * up to this point.
   *
   * \return The number of iterations
   */
  virtual std::size_t GetNumberOfIterations() const noexcept;

  /**
   * This method returns the number of function evaluations the optimizer has
   * performed up to this point.
   *
   * \return The number of function evaluations
   */
  virtual std::size_t GetNumberOfFunctionEvaluations() const noexcept;

  /**
   * This method returns the stored population with the costs appended to each
   * member. There is no guarantee that the population will be sorted.
   *
   * \return The population
   */
  virtual std::vector<std::vector<double>> GetPopulation() const;

  /**
   * This method returns a member of the current solution. For optimizers where
   * more than one solution is stored (e.g. a population) the optional parameter
   * can be used to choose which solution to get. The solutions are not
   * guaranteed to be sorted. If an invalid member index is passed in as an
   * argument the returned member will be empty. The cost is not appended to the
   * end of the member that is returned.
   *
   * DifferentialEvolution: Members are unsorted \n
   * GeneticAlgorithm: Members are sorted (best = 0, worst = n) \n
   * LevenbergMarquardt: Only one member is stored (index = 0) \n
   * NelderMead: Members are sorted (best = 0, worst = n)
   *
   * \param index The index of the member to be returned
   * \return The specified member (sans cost)
   */
  virtual std::vector<double> GetSolution(std::size_t index = 0) const;

  /**
   * This method allows the user to set the population directly, rather than
   * having it generated internally. Different methods will have a different
   * minimum population size. In NelderMead, for example, the population size
   * must be exactly one more than your number of dimensions (to get a simplex).
   * Here the population that is passed in is not tested. Each optimizer will
   * implement its own checks to make sure the population that is passed in is
   * valid. Note also that the population should be passed in with the costs
   * appended.
   *
   * \param population The population to be used (with costs appended)
   */
  virtual void SetPopulation(
      const std::vector<std::vector<double>> &population);

  /**
   * Use this object to set and manipulate bounds for your problem, if any.
   */
  Unfit::Bounds bounds;

  /**
   * Use this object to set and manipulate the optimization options, such as
   * the maximum number of iterations, tolerances and output levels.
   */
  Unfit::Options options;

 protected:
  /**
   * This method calls a cost function that has been derived from
   * GenericCostFunction with a specified set of parameters. GenericCostFunction
   * returns a vector of residuals (as is needed by the LevenbergMarquardt
   * algorithm) and this function uses a sum-of-squared residuals approach to
   * convert the vector of residuals into a single cost value. This method also
   * increments the number of function evaluations.
   *
   * \param CostFunction The cost function to be evaluated
   * \param x The parameters with which the cost function will be evaluated
   * \return True if the cost is finite, otherwise false
   */
  virtual bool CalculateCost(GenericCostFunction &CostFunction,
      std::vector<double> &x);

  /**
   * This method creates a population that follows a uniform random distribution
   * between the specified bounds. It is recommended that you set sensible
   * bounds as the defaults are +/- max double, which can make population
   * generation take a very long time. Costs are calculated via the provided
   * cost function and are appended to the end of each population member. Any
   * member that does not compute a valid cost is automatically discarded and
   * regenerated. The population is not sorted within this method.
   *
   * \param CostFunction The function used to evaluate the cost of each member
   * \param dimensions The number of dimensions for each member
   */
  void GeneratePopulation(GenericCostFunction &CostFunction,
      std::size_t dimensions);

  /**
   * This method generates set of random number engines (one for each member of
   * the population) based on a random seed. The random seed provided in the
   * options.SetRandomSeed() is then used to generate random seeds for each of
   * the engines. There needs to be one per member so all of the member
   * calculations can be done in parallel.
   */
  void GenerateRandomEngines();

  /**
   * \brief Checks to see if the population has converged
   *
   * This method checks to see if the population has converged. In particular,
   * it checks three things - the best cost, the spread of costs in the
   * population, and the spread of the population coordinates in each dimension.
   * In order to be considered converged, either the best cost has to be less
   * than the cost tolerance, or BOTH the cost and the geometry have to be
   * converged to the same point. Details on the second case are given below.
   *
   * Cost: Uses options.GetCostTolerance(). If the best cost is < 1.0 then the
   * difference between the best and worst costs must be less than the cost
   * tolerance. If the best cost is > 1.0 then the difference between the best
   * and worst costs must be less than the best_cost*cost_tolerance.
   *
   * Geometry: Uses options.GetGeometricTolerance(). In each dimension, if the
   * position with the best cost is < 1.0 then the difference between all
   * member locations and the best member location must be less than the
   * geometric tolerance. If the position with the best cost (again separately
   * for each dimension) is > 1.0 then the difference between the best location
   * and all other member locations must be less than
   * best_location*geometric_tolerance.
   *
   * \param best_member The best member of the current population (with cost)
   * \param truncated_index This parameter should be used if you only want to
   *        check part of the population for convergence. An example of this is
   *        GeneticAlgorithm which only checks the surviving members of the
   *        population for convergence. If you have a population 0..7, and you
   *        want to check only members 0,1,2,3, then truncated_index will 3. A
   *        truncated_index of zero will check the whole population.
   * \return true if both convergence criteria are met, otherwise false
   */
  virtual bool IsConverged(const std::vector<double> &best_member,
      std::size_t truncated_index = 0) const;

  /**
   * Prints the column headers relating to iteration by iteration output to the
   * screen if the output level has been set to be >= 1.
   *
   * \param best_cost The best cost after the initial setup
   */
  virtual void PrintInitialOutput(double best_cost) const;

  /**
   * At each iteration, prints nothing if the output level is zero, the
   * iteration number if the output level is one, and the iteration number plus
   * the number of function evaluations and the best cost if the output level is
   * greater than one.
   *
   * \param best_cost The best cost after the latest iteration
   */
  virtual void PrintIterationOutput(double best_cost) const;

  /**
   * Upon completion, this method prints the final population if the output
   * level has set to be greater than two.
   */
  virtual void PrintFinalOutput() const;

  /**
   * This method sorts the members of the population into ascending order by
   * their cost, such that the first member (index = 0) contains the lowest cost
   * and the last member contains the highest cost.
   */
  virtual void SortPopulation() noexcept;

  /**
   * A vector of population member vectors. The cost of a given member is
   * appended to that member, so the size will be the number of dimensions in
   * the problem + 1.
   */
  std::vector<std::vector<double>> population_;
  /** A vector containing one random engine for each member of the population */
  std::vector<std::mt19937> random_engines_;
  /** Variable to store the number of function evaluations*/
  std::atomic<std::size_t> function_evaluations_;
  /** Variable to store the current number of iterations*/
  std::atomic<std::size_t> iterations_;
  /** Variable to store if the optimizer is population based*/
  bool is_population_based_;

 private:
  /**
   * This method generates a single population member for the GeneratePopulation
   * method. Each parameter in the member is generated between the specified
   * bounds (with the default bounds being +/- max(double), and each has its
   * own random number generator (to allow parallel execution)
   *
   * \param CostFunction The cost function to be optimised
   * \param distributions A uniform distribution with the correct bounds for
   *        each dimension
   * \param dimensions The number of dimensions in the generated member
   * \param mem The index of the population member to be generated
   * \return A vector containing the new population member
   */
  std::vector<double> GeneratePopulationMember(
      GenericCostFunction &CostFunction,
      std::vector<std::uniform_real_distribution<double>> &distributions,
      std::size_t dimensions, std::size_t mem);
};

}  // namespace Unfit

#endif

