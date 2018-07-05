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
#ifndef UNFIT_EXAMPLES_OSBORNE_HPP_
#define UNFIT_EXAMPLES_OSBORNE_HPP_

#include <cmath>
#include <vector>
#include "GenericCostFunction.hpp"

namespace Unfit
{
namespace Examples
{
/**
 * \brief Fit a double exponential to experimental data
 *
 * Here the goal is to find a parameter set that best fits the following
 * function to the experimental data:
 *
 *   f(t) = A + B*exp(-D*t) + C*exp(-E*t)
 *
 * The goal is to find the values of A, B, C, D & E that gives a best fit. In
 * terms of the model, A = param[0], B = param[1], C = param[2], D = param[3],
 * and E = param[4].
 */
class Osborne : public GenericCostFunction
{
 public:
  /**
   * Create the cost function. Here the experimental data must be passed in,
   * and cannot be changed (if you want to, just create another cost function
   * object). Here the experimental data is a vector of data.
   *
   * Intended use :
   *   Osborne cost_func(x);
   *
   * \param x A vector of experimental data
   */
  Osborne(const std::vector<double> &x)
  : x_ {x}
  {}

  /**
   * Calculate the linear distance (residuals) between our model and the data.
   * This method encapsulates the model, and expects the current estimates of
   * the unknown parameters as an input. See the class documentation for details
   * about the model.
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
    auto residuals = x_;
    for (auto i = 0u; i < residuals.size(); ++i) {
      double t = 10.0 * static_cast<double>(i);
      residuals[i] -= param[0] + param[1]*exp(-param[3]*t) +
          param[2]*exp(-param[4]*t);
    }
    return residuals;
  }
 private:
  /** A vector to store the experimental data x*/
  const std::vector<double> x_;
};

}  // namespace Examples
}  // namespace Unfit

#endif
