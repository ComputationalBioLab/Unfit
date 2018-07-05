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
#ifndef UNFIT_EXAMPLES_ODE2D_HPP_
#define UNFIT_EXAMPLES_ODE2D_HPP_

#include <vector>
#include "GenericCostFunction.hpp"

namespace Unfit
{
namespace Examples
{
/**
 * \brief Fit a two parameter first order ordinary differential equation (ODE)
 *
 * This is an example data fitting problem with two parameters for a single
 * first order ordinary differential equation (ODE). Here the goal is to fit
 * the following ODE to some experimental data:
 *
 *   dx/dt = Ax + B
 *
 * The goal is to find the values of A & B that best fit the data. In terms
 * of the model, A = param[0] and B = param[1].
 */
class ODE2D : public GenericCostFunction
{
 public:
  /**
   * Create the cost function. Here the experimental data must be passed in,
   * and cannot be changed (if you want to, just create another cost function
   * object). Here the experimental data is a vector, and we also require the
   * time step between observations.
   *
   * Intended use :
   *   ODE2D cost_func(x, dt);
   *
   * \param x A vector of experimental data
   * \param dt The time step between observations
   */
  ODE2D(const std::vector<double> &x, double dt)
    : x_ {x},
      dt_ {dt}
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
    double predicted_x = 0.0;  // Initial condition
    for (auto i = 0u; i < residuals.size(); ++i) {
      residuals[i] -= predicted_x;
      predicted_x += dt_ * (param[0] * predicted_x + param[1]);
    }
    return residuals;
  }
 private:
  /** A vector to store the experimental data x*/
  const std::vector<double> x_;
  /** The time step between observations*/
  const double dt_;
};

}  // namespace Examples
}  // namespace Unfit

#endif
