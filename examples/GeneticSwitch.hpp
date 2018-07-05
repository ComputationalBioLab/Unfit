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
#ifndef UNFIT_EXAMPLES_GENETICSWITCH_HPP_
#define UNFIT_EXAMPLES_GENETICSWITCH_HPP_

#include <cmath>
#include <vector>
#include "GenericCostFunction.hpp"

namespace Unfit
{
namespace Examples
{
/**
 * \brief Fit the parameters in a model of a genetic switch
 *
 * An example two ordinary differential equations (ODE) data fitting problem.
 * The equations describe genetic switch between two genes, z_a and z_b, for
 * protein regulation. They depend on four constants, u, n, v1 and v2. Here the
 * goal is to fit these constants given experimental data for z_a. The two ODEs
 * are:
 *
 *   dz_a/dt = u/(1 + z_b^n) - z_a - v1
 *   dz_b/dt = u/(1 + z_b^n) - z_b - v2
 *
 * The initial conditions are z_a(0) = 10, and z_b(0) = 10, and the goal is to
 * find the values of u, n, v1 & v2 that best fit the data. In terms
 * of the model, u = param[0], n = param[1], v1 = param[2] and v2 = param[3].
 */
class GeneticSwitch : public GenericCostFunction
{
 public:
  /**
   * Create the cost function. Here the experimental data must be passed in,
   * and cannot be changed (if you want to, just create another cost function
   * object). Here the experimental data is a vector containing the value
   * of z_a, and the time step used in the data collection.
   *
   * Intended use :
   *   GeneticSwitch cost_func(z_a, t);
   *
   * \param z_a A vector of experimental z_a values
   * \param dt The experimental time step
   */
  GeneticSwitch(const std::vector<double> &z_a, double dt)
    : z_a_data_ {z_a},
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
  std::vector<double> operator() (const std::vector<double> &param)
  {
    double z_a = 10.0;          // Initial condition
    double z_b = 10.0;          // Initial condition
    auto residuals = z_a_data_;    // Residuals
    residuals[0] -= z_a;
    for (auto i = 1u; i < residuals.size(); ++i) {
      const auto l_z_a = z_a + dt_ * (param[0] / (1.0 + pow(z_b, param[1])) -
          z_a - param[2]);
      z_b += dt_ * (param[0] / (1.0 + pow(z_a, param[1])) - z_b - param[3]);
      z_a = l_z_a;
      residuals[i] -= z_a;
    }
    return residuals;
  }
 private:
  /** A vector to store the experimental data, Z_A */
  const std::vector<double> z_a_data_;
  /** A vector to store the time step, dt */
  const double dt_;
};

}  // namespace Examples
}  // namespace Unfit

#endif
