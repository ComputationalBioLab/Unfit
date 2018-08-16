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
#ifndef UNFIT_EXAMPLES_ODE2DMODEL_HPP_
#define UNFIT_EXAMPLES_ODE2DMODEL_HPP_

#include <vector>
#include "GenericModel.hpp"

namespace Unfit
{
namespace Examples
{
/**
 * This is model that contains a single first order, ordinary differential
 * equation (ODE). The equation has two parameters, c0, and c1, and looks like:
 *
 *   dy/dt = c0*y + c1
 *
 * If we choose to solve this equation numerically with a forward Euler
 * approach, we can discretize the derivative term in the ODE above and
 * rearrange to get an update equation:
 *
 *   y(t+dt) = y(t) + dt*[c0*y(t) + c1]
 *
 * It is this discretized equation that is implemented in the model. Of course
 * you could implement another more sophisticated integration scheme (e.g.
 * Improved Euler or Runge-Kutta) if you wanted to. The result is an initial
 * value problem, so we require an initial value for y (i.e., y(t=0)) to get
 * started.
 */
class Ode2DModel : public GenericModel
{
 public:
  /**
   * Create the model with the default parameters. Here the default time step
   * step is set to 1.0 and the default initial condition is set to 0.0.
   */
  Ode2DModel()
    : dt_ {1.0},
      ic_ {0.0}
  {}

  /**
   * Create the model and provide the time step that we should be using to
   * solve the equation, along with the value of the initial condition for y
   * i.e., y(t=0).
   *
   * \param dt The time step between observations
   * \param ic The value of the initial condition
   */
  Ode2DModel(double dt, double ic)
    : dt_ {dt},
      ic_ {ic}
  {}

  /**
   * This method takes in a vector of model parameters (c), and a vector of
   * times (t). While we could use the times to calculate dt at each iteration,
   * here we assume that the dt passed in during the model construction is
   * constant and correct, so we don't make use of the t vector, apart from
   * using it to get the number of time steps/observations. A three parameter
   * version of this could have the initial condition as an additional unknown.
   *
   * \param c A vector of model parameters/constants
   * \param t A vector of independent variable values (stored in t[0])
   * \return A vector containing the model solution at each time
   */
  std::vector<double> operator()(const std::vector<double> &c,
      const std::vector<std::vector<double>> &t)
  {
    auto model = t[0];  // Create a vector to store the model results
    model[0] = ic_;     // Impose the initial condition
    // Loop over time and solve via forward Euler
    for (auto i = 1u; i < model.size(); ++i) {
      model[i] = model[i-1] + dt_ * (c[0] * model[i-1] + c[1]);
    }
    return model;
  }
 private:
  /** The time step between observations*/
  const double dt_;
  /** The initial condition for the problem*/
  const double ic_;
};

}  // namespace Examples
}  // namespace Unfit

#endif  // UNFIT_EXAMPLES_ODE2DMODEL_HPP_

