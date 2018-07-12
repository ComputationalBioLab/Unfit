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
#ifndef UNFIT_INCLUDE_LEVENBERGMARQUARDT_HPP_
#define UNFIT_INCLUDE_LEVENBERGMARQUARDT_HPP_

#include <vector>
#include "Bounds.hpp"
#include "GenericCostFunction.hpp"
#include "GenericOptimizer.hpp"
#include "Matrix.hpp"
#include "Options.hpp"

namespace Unfit
{
/**
 * \brief A class to implement the Levenberg Marquardt optimization method
 *
 * This class implements the Levenberg Marquardt optimization method. It
 * requires a cost function (written by the user) and an initial guess of the
 * unknown parameters. It then proceeds to try to find a set of parameters that
 * minimize the cost function.
 *
 * The inspiration for this implementation was derived from the book:
 *   Methods for Non-linear Least Squares Problems,
 *   by K. Madsen, H.B. Nielsen, O. Tingleff.
 * At the time of writing, this was freely available online. There is also a
 * sample matlab implementation on their website that was useful for getting the
 * Broyden rank one updates to work.
 */
class LevenbergMarquardt : public GenericOptimizer
{
  /** For unit testing purposes only*/
  friend class TestLevenbergMarquardt;
 public:
  /**
   * The constructor sets all of the parameters to their defaults.
   */
  LevenbergMarquardt();

  /**
   * As we are deriving from this class, the destructor should be virtual. In
   * this case an empty destructor is just fine as we are putting everything
   * into std::vectors which take care of themselves.
   */
  virtual ~LevenbergMarquardt() {}

  /**
   * \brief Implements the Levenberd-Marquardt optimization method.
   *
   * This method takes in a cost function (user supplied) and an initial guess
   * for the unknown parameters and uses the Levenberd-Marquardt algorithm to
   * try to find the set of parameters that provides a minimum cost. This is
   * implemented as a nonlinear least squares problem, so the cost will always
   * be >= 0. This is known as a local optimizer, so the results you get will
   * depend on how good the initial guess is.
   *
   * Intended use:
   *   auto rc = object.FindMin(cost_function, coordinates)
   *
   * \param cost_function residuals to be tested
   * \param coordinates vector of initial guesses
   * \return 0 if successful \n
   *         1 if empty coordinates are passed in \n
   *         2 if the initial coordinates are out of bounds \n
   *         3 if the the data_points < number of variables \n
   *         4 if the maximum number of function evaluations was exceeded \n
   *         5 if the maximum number of iterations has been exceeded \n
   *         6 if the initial cost is not finite \n
   *         7 if the solution breaks down (nu gets very large)
   */
  int FindMin(GenericCostFunction &cost_function,
      std::vector<double> &coordinates) override;

  /**
   * Resets everything back to their default values. This is equivalent to
   * destroying the object and creating a new one.
   */
  void Reset() override;

 private:
  /**
   * This method calculates the new coordinates based on the old coordinates
   * and the step size, hlm. If the new coordinates fall within the bounds of
   * the domain that is all it does. If the new coordinates are outside the
   * bounds of the domain, the step size is iteratively reduced by a factor of
   * two until it is inside the domain.
   *
   * \param point current coordinates
   * \param hlm step to be taken
   * \return the new coordinates (point + hlm)
   */
  std::vector<double> LoopUntilWithinBounds(
      const std::vector<double> &point, std::vector<double> &hlm);

  /**
   * Calculate the residuals (costs) associated with the the current estimate
   * of the unknown parameters. NOTE that you should never call the cost
   * function directly, call this method instead. This method also tracks the
   * number of function evaluations.
   *
   * \param cost_function the function to be used to calculate the residuals
   * \param coordinates the estimate of the unknown parameters
   * \return a vector containing the residuals associated with the current
   *         coordinates
   */
  std::vector<double> CalculateResiduals(GenericCostFunction &cost_function,
      const std::vector<double> &coordinates);

  /**
   * This method computes the first derivatives of the cost function with
   * respect to each of the unknown parameters and stores the values in the
   * jacobian matrix. Note that it is more efficient to store JT
   * (J-transpose) so that is what is done here.
   *
   * \param cost_function cost functions to calculate the residuals
   * \param residuals a vector of residuals at the current guess. Only used for
   *        one-sided differencing; can be empty for two-sided differencing.
   * \param coordinates vector of the coefficients of the function
   * \param one_sided_difference select between one sided (forward) and two
   *        sided finite different approximations for the Jacobian calculation
   */
  void FindJacobian(GenericCostFunction &cost_function,
      const std::vector<double> &residuals,
      const std::vector<double> &coordinates, bool one_sided_difference = true);

  /**
   * Rather than calculate the jacobian matrix via finite differences at each
   * iteration, for expensive cost functions in particular, it is worth using
   * Brodyen Rank One Updates. These given an estimate of the step size and
   * do not require any function evaluations. See the text quoted in the class
   * documentation for details.
   *
   * \param residuals_new the residuals at the proposed parameter estimate
   * \param residuals the residuals at the current parameter estimate
   * \param hlm the step size between the current and new parameter estimates
   * \param jacobianT the transpose of the jacobian matrix (updated)
   */
  void BroydenUpdate(const std::vector<double> &residuals_new,
      const std::vector<double> &residuals, const std::vector<double> &hlm,
      Matrix &jacobianT);

  /** Jacobian matrix */
  Matrix jacobian_;
  // MB: I need to have a look to see if it is much slower if we just make the
  // jacobian locally each time we need it, rather than storing it as a member
  // variable (as it is probably cleaner - more pure).
  /** Number of parameters to find */
  std::size_t dimensions_;
  /** Number of experimental observations */
  std::size_t observation_size_;
  /** Variable to store the current cost of the solution*/
  double stored_cost_;
};

}  // namespace Unfit

#endif
