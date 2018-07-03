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
#ifndef UNFIT_INCLUDE_DIFFERENTIALEVOLUTION_HPP_
#define UNFIT_INCLUDE_DIFFERENTIALEVOLUTION_HPP_

#include <vector>
#include "GenericCostFunction.hpp"
#include "GenericOptimizer.hpp"

namespace Unfit
{
/**
 * \brief A class to implement the Differential Evolution optimization method
 *
 * This class implements the Differential Evolution optimization method. It
 * requires a cost function (written by the user) and an initial guess of the
 * unknown parameters. The initial guess is only needed to obtain the number
 * of unknowns, so make sure it is of the correct length. The algorithm will
 * proceed to try to find a set of parameters that minimize the cost function.
 *
 * The idea behind Differential Evolution comes from Kenneth Price and
 * Rainer Storn. Inspiration for this implementation was derived from their
 * work, in particular see: http://www1.icsi.berkeley.edu/~storn/code.html
 */
class DifferentialEvolution : public GenericOptimizer
{
  friend class TestDifferentialEvolution;
 public:
  /**
   * The constructor sets all of the parameters to their default values.
   */
  DifferentialEvolution();

  /**
   * As we could derive from this class, the destructor should be virtual. In
   * this case the default destructor is just fine as we are putting everything
   * into std::vectors which take care of themselves.
   */
  virtual ~DifferentialEvolution();

  /**
   * \brief A method to find the minimum of a model/function using a
   * Differential Evolution approach.
   *
   * This is the main method to call to perform a minimization of the provided
   * cost function. Here a Differential Evolution optimization technique is
   * adopted. There are a few points worth noting before starting the
   * optimization. First, make sure that the length of the coordinates vector
   * that is passed in to this function matches the number of unknowns you have
   * in your cost function. As the cost function is user supplied there is no
   * way of checking these match up. The values you supply are not important as
   * the initial population is randomly generated. The final result (coordinates
   * with best fit) is returned in the vector. Second, it is advisable to add
   * bounds to each of the variables in the fit. The algorithm generates
   * possible answers between the provided bounds, and if no bounds are provided
   * the bounds are +/- the maximum double the computer can represent. Chances
   * are you know enough about the problem to make more sensible choices.
   * Convergence with +/- 1e8, for example, is still much faster than the
   * default bounds.
   *
   * Intended use:
   *   int rc = object.FindMin(CostFunction, coordinates)
   *
   * \param CostFunction The user-supplied function to minimise
   * \param coordinates On input must be the same size as the number of
   *        parameters in the cost function. On exit returns the unknowns
   *        that provided the minimum cost.
   * \return Zero upon success, non-zero upon failure. \n
   *        -1 = Empty coordinate vector \n
   *         1 = Maximum number of function evaluations was reached \n
   *         2 = Maximum number of iterations was reached
   */
  int FindMin(GenericCostFunction &CostFunction,
      std::vector<double> &coordinates);

  /**
   * Resets DifferentialEvolution back to its default state. This is equivalent
   * to destroying the object and creating a new one.
   */
  void Reset();

 private:
  /**
   * This method calls GenerateTrialMember to get a proposed addition to the
   * population. If the proposed member has a better cost it is returned to be
   * added to the population, otherwise the existing member is returned such
   * that this member remains unchanged at this iteration.
   *
   * \param CostFunction The cost function to be optimised
   * \param member The index of the population member to be chosen
   * \return A vector containing the new population member
   */
  std::vector<double> NewPopulationMember(GenericCostFunction &CostFunction,
      unsigned member);

  /**
   * This method generates a trial member. Once a trial member has been obtained
   * the trial member is compared to the relevant member of the current
   * population to decide if it is a worthy replacement (done in
   * ProcessFindMin). There are a number of strategies that can be employed to
   * generate a trial member. Following Price & Storn, those implemented here
   * are (in order, see SetStrategy):
   *
   *   - best/1/exp
   *   - rand/1/exp
   *   - rand-to-best/1/exp
   *   - best/2/exp
   *   - rand/2/exp
   *   - best/1/bin
   *   - rand/1/bin
   *   - rand-to-best/1/bin
   *   - best/2/bin
   *   - rand/2/bin
   *
   * The leading best or rand or rand-to-best is the base member to which a
   * weighted average of other members is added. The number refers to
   * how many pairs (1 or 2) of population members are used to calculate the
   * increment fo the base member. Finally, the trailing exp or bin refers
   * to an exponential or binary cross over method.
   *
   * \param i the index of the population member for which a trial member will
   *        be generated
   * \return a vector containing the trial member
   */
  std::vector<double> GenerateTrialMember(unsigned i);

  /**
   * This method is called by FindMin and contains the main iteration loop. Here
   * is where the actual optimisation takes place, the population is mutated by
   * vector combinations and convergence is determined.
   *
   * \param CostFunction Returns the residuals of the model
   * \return Zero upon success, non-zero upon failure \n
   *         1 = Maximum number of function evaluations was reached \n
   *         2 = Maximum number of iterations was reached
   */
  int ProcessFindMin(GenericCostFunction &CostFunction);

  /** The next population, constructed then swapped with population_*/
  std::vector<std::vector<double>> new_population_;
  /** A vector that contains the best solution found so far*/
  std::vector<double> best_member_;
  /**variable to store the number of dimensions*/
  unsigned dimensions_;
  /**the index of the cost of each vertex*/
  unsigned cost_;
};

}  // namespace Unfit

#endif
