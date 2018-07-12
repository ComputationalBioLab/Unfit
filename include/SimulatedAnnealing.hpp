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
#ifndef UNFIT_INCLUDE_SIMULATEDANNEALING_HPP_
#define UNFIT_INCLUDE_SIMULATEDANNEALING_HPP_

#include <vector>
#include "GenericCostFunction.hpp"
#include "GenericOptimizer.hpp"

namespace Unfit
{
/**
 * \brief A class to implement the Simulated Annealing optimization method
 *
 * This class implements the Simulated Annealing optimization method. It
 * requires a cost function (written by the user) and an initial guess of the
 * unknown parameters. The initial guess is only needed to obtain the number
 * of unknowns, so make sure it is of the correct length. The algorithm will
 * proceed to try finding a set of parameters that minimizes the cost function.
 *
 * The inspiration for this implementation was derived from the book:
 * Belegundu, A. D., & Chandrupatla, T. R. (2011). Optimization concepts and
 * applications in engineering. Cambridge University Press.
 * In particular, look at Chapter 7 for their recommendations of the Simulated
 * Annealing algorithm implementation.
 */
class SimulatedAnnealing : public GenericOptimizer
{
  /** For unit testing purposes only*/
  friend class TestSimulatedAnnealing;
 public:
  /**
   * A constructor that sets all parameters to their default values.
   */
  SimulatedAnnealing();

  /**
   * As we could derive from this class, the destructor should be virtual. In
   * this case, the default destructor is just fine as we are putting everything
   * into std::vectors, which would take care of themselves.
   */
  virtual ~SimulatedAnnealing() = default;

  /**
   * \brief A method to find a minimum point of a function using a Simulated
   * Annealing approach.
   *
   * This is the main method to call so as to perform a minimization of the
   * provided cost function. Here, a Simulated Annealing optimization technique
   * is adopted. There are a few points worth noting before starting the
   * optimization. First, make sure that the length of the coordinates vector
   * passed in to this function matches the number of unknowns you have in your
   * cost function. As the cost function is supplied by user, there is no way of
   * checking if the two match up. By default, the vector you supply is not used
   * (other than to see it's length, and a random initial guess is generated
   * within the specified bounds. However, if you want to start from your guess
   * simply call:
   *
   *   sa.options.SetAddInitialToPopulation(true);
   *
   * before you call FindMin, and we will use yours. Next, it is highly
   * advisable to add bounds to each of the variables in the fit. The algorithm
   * generates possible answers between the provided bounds, and if no bounds
   * are provided the bounds are +/- the maximum double the computer can
   * represent (read - slow). Chances are you know enough about the problem to
   * make a more sensible choice. Lastly, the final result (answer with the
   * best fit) is returned in the coordinates vector.
   *
   * Intended use:
   *   int rc = object.FindMin(cost_function, coordinates)
   *
   * \param CostFunction The user-supplied function to minimize
   * \param coordinates On input, it must be the same size as the number of
   *        unknown parameters in the cost function. On exit, it returns the
   *        values of the unknown parameters that provided the minimum cost.
   * \return Zero upon success, non-zero upon failure. \n
   *        -1 = Empty coordinate vector \n
   *        -2 = Initial user supplied coordinates are out of bounds \n
   *        -3 = Cost of the initial coordinates is inf or nan \n
   *         1 = Maximum number of function evaluations was reached \n
   *         2 = Maximum number of iterations was reached
   */
  int FindMin(GenericCostFunction &CostFunction,
      std::vector<double> &coordinates) override;

  /**
   * Resets the Simulated Annealing algorithm back to its default state. This is
   * equivalent to destroying the object and creating a new one.
   */
  void Reset() override;

 private:
  /**
   * This method is called by FindMin and it contains the main iteration loop
   * where the actual optimization takes place.
   *
   * \param CostFunction The user-supplied function to minimize
   * \param coordinates On input, it must be the same size as the number of
   *        unknown parameters in the cost function. On exit, it returns the
   *        values of the unknown parameters that provided the minimum cost.
   * \return Zero upon success, non-zero upon failure. \n
   *         1 = Maximum number of function evaluations was reached \n
   *         2 = Maximum number of iterations was reached
   */
  int ProcessFindMin(GenericCostFunction &CostFunction,
      std::vector<double> &coordinates);

  /**
   * Once we know the dimensions of the problem, this method is called to
   * initialize the relevant private member variables.
   */
  void InitialiseParameters();

  /**
   * This method generates a trial point from the current vector of coordinates
   * by perturbing one of the values. The trial point, perturbing direction i,
   * is generated by:
   *
   *   trial_point[i] += random_number*step_size[i]
   *
   * where random_number is a uniform random number in the range of -1 to 1.
   *
   * \param trial_point The coordinate vector to be perturbed
   * \param i The coordinate/direction to perturb
   */
  void GenerateTrialPoint(std::vector<double> &trial_point, int i);

  /**
   * This method updates the step sizes after a cycle has been completed. The
   * value of the acceptance ratios are used to update the step sizes. After the
   * step sizes have been adjusted, the acceptance ratios are reset to 1.
   */
  void UpdateStepSizes() noexcept;

  /**
   * The step size may take a low values eventually and may not allow a step to
   * accept a point with higher cost even at moderate temperatures. To overcome
   * this issue, this method resets each step size at the end of each
   * temperature loop.
   *
   * \param step_size A scaling factor for setting the step sizes. The actual
   *        updated size will be this factor multiplied by the difference
   *        between the upper and lower bounds for each variable.
   */
  void ResetStepSizes(double step_size) noexcept;

  /** The index of the cost of each point*/
  std::size_t cost_;
  /** Variable to store the number of dimensions*/
  std::size_t dimensions_;
  /** The cost from the previous iteration*/
  double previous_best_cost_;
  /** A vector of step sizes*/
  std::vector<double> step_sizes_;
  /**A vector of acceptance ratios*/
  std::vector<double> acceptance_ratios_;
  /** A random number engine*/
  std::mt19937 generator_;
  /** A uniform distribution over the interval [0, 1]*/
  std::uniform_real_distribution<double> uniform_dist_;
};

}  // namespace Unfit

#endif
