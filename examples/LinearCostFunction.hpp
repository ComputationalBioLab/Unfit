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
#ifndef UNFIT_EXAMPLES_LINEARCOSTFUNCTION_HPP_
#define UNFIT_EXAMPLES_LINEARCOSTFUNCTION_HPP_

#include <vector>
#include "GenericCostFunction.hpp"

namespace Unfit
{
namespace Examples
{
/**
 * A model of a straight line, which has two parameters. The function is
 * given by:
 *
 *   f(x) = c0*x + c1
 *
 * Here we assume we have some observations, y_i at locations x_i. We want to
 * compute the distance (difference) between the model, f(x_i) and the data
 * at y_i. We do this as:
 *
 *   r_i = y_i - f(x_i)
 *
 * r_i is a residual, and we calculate one residual for each data point then
 * return them in a vector.
 */
class LinearCostFunction : public GenericCostFunction
{
 public:
  /**
   * As we are storing data internally, remove the possibility of creating the
   * cost function without any data.
   */
  LinearCostFunction() = delete;

  /**
   * Create the cost function. Here the experimental data is passed in as two
   * vectors, the independent variable (x), and the dependent variable (y). The
   * data is stored within the class as the optimisers will call this thousands
   * of times and we only want to pass the data in once.
   *
   * \param x A vector of independent variable values
   * \param y A vector of experimental observations
   */
  LinearCostFunction(const std::vector<double> &x, const std::vector<double> &y)
    : x_ {x},
      y_ {y}
  {}

  /**
   * Calculate the linear distance (residuals) between our model and the data.
   * Given data y[] and a model m[], we return r[] = y[] - m[]. This method
   * therefore encapsulates the model, and expects the current estimates of
   * the unknown parameters c[] as an input. Here the model is a straight line.
   * For simplicity, we set r[] = y[], then calculate r[] - m[].
   *
   * \param c A vector containing the model parameters
   * \return A vector containing the resulting residuals
   */
  std::vector<double> operator()(const std::vector<double> &c)
  {
    auto residuals = y_;
    for (auto i = 0u; i < residuals.size(); ++i) {
      residuals[i] -= (c[0]*x_[i] + c[1]);
    }
    return residuals;
  }
 private:
  /** A vector to store the experimental x values*/
  const std::vector<double> x_;
  /** A vector to store the experimental y values*/
  const std::vector<double> y_;
};

}  // namespace Examples
}  // namespace Unfit

#endif  // UNFIT_EXAMPLES_LINEARCOSTFUNCTION_HPP_
