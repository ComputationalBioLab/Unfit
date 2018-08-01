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
#ifndef UNFIT_EXAMPLES_ODE3DMODEL_HPP_
#define UNFIT_EXAMPLES_ODE3DMODEL_HPP_

#include <vector>
#include "GenericNDModel.hpp"

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
 * However, we can see this is an initial value problem and we will require an
 * initial value for y to solve this (i.e., y(t=0)). Therefore we add a third
 * parameter, c2, to the problem as follows:
 *
 *   y(t=0) = c2
 *
 * If we choose to solve this equation numerically with a forward Euler
 * approach, we can discretize the derivative term in the ODE above and
 * rearrange to get an update equation:
 *
 *   y(t+dt) = y(t) + dt*[c0*y(t) + c1]
 *
 * It is this discretized equation that is implemented in the model. Of course
 * you could implement another more sophisticated integration scheme (e.g.
 * Improved Euler or Runge-Kutta) if you wanted to.
 */
class Ode3DModel : public GenericNDModel
{
 public:
  /**
   * This method takes in a vector of model parameters (c), and a vector of
   * times (t). The time step used in the integration is calculated from pairs
   * of adjacent entries in the t vector. Note that the third of the three
   * parameters represents the initial condition for this initial value problem.
   *
   * \param c A vector of model parameters/constants
   * \param t A vector of time values (stored in t[0])
   * \return A vector containing the model solution at each time
   */
  std::vector<double> operator()(const std::vector<double> &c,
      const std::vector<std::vector<double>> &t)
  {
    auto model = t[0];  // Create a vector to store the model results
    model[0] = c[2];    // Impose the initial condition
    // Loop over time and solve via forward Euler
    for (auto i = 1u; i < model.size(); ++i) {
      // Get the appropriate time step from the t vector
      const auto dt = t[0][i] - t[0][i-1];
      model[i] = model[i-1] + dt * (c[0] * model[i-1] + c[1]);
    }
    return model;
  }
};

}  // namespace Examples
}  // namespace Unfit

#endif  // UNFIT_EXAMPLES_ODE3DMODEL_HPP_
