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
#ifndef UNFIT_EXAMPLES_PARABOLIC_HPP_
#define UNFIT_EXAMPLES_PARABOLIC_HPP_

#include <vector>
#include "GenericCostFunction.hpp"

namespace Unfit
{
namespace Examples
{
/**
 * \brief Fit a three parameter parabola to some data
 *
 * Here the goal is to fit a parabola with three parameters to some data. The
 * function is given by:
 *
 *   f(x) = A*x*x + B*x + C
 *
 * The goal is to find the values of A, B & C that best fit the data. In terms
 * of the code, A = param[0], B = param[1] and C = param[2].
 */
class Parabolic : public GenericCostFunction
{
 public:
  /**
   * Create the cost function. Here the experimental data is passed in as two
   * vectors, the independent variable (x), and the dependent variable (y). The
   * data is stored within the class so it only has to be initialized once.
   *
   * Intended use :
   *   Parabolic cost_func(data_x, data_y);
   *
   * \param data_x A vector of independent variable values
   * \param data_y A vector of experimental observations
   */
  Parabolic(const std::vector<double> &data_x,
      const std::vector<double> &data_y)
    : data_x_ {data_x},
      data_y_ {data_y}
  {}

  /**
   * Calculate the linear distance (residuals) between our model and the data.
   * Given data y[] and a model m[], we return r[] = y[] - m[]. This method
   * therefore encapsulates the model, and expects the current estimates of
   * the unknown parameters as an input. All cost functions must have this
   * method and must return a vector of residuals. The same cost function can
   * be used for any of the optimization techniques.
   *
   * Intended use :
   *   residuals = cost_func(param)
   *
   * \param param A vector containing the current estimates of the parameters
   *     we are trying to fit
   * \return A vector containing the residuals
   */
  std::vector<double> operator()(const std::vector<double> &param)
  {
    // First set the residuals to be the observed values
    auto residuals = data_y_;
    // Next, subtract the model evaluated at each point to give the residuals
    // as r[i] = y[i] - (A*x[i]*x[i] + B*x[i] + C)
    for (auto i = 0u; i < residuals.size(); ++i) {
      residuals[i] -= (param[0]*data_x_[i]*data_x_[i] + param[1]*data_x_[i]
          + param[2]);
    }
    return residuals;
  }
 private:
  /** A vector to store the experimental x values*/
  const std::vector<double> data_x_;
  /** A vector to store the experimental y values*/
  const std::vector<double> data_y_;
};

}  // namespace Examples
}  // namespace Unfit

#endif
