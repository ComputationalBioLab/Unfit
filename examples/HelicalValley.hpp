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
#ifndef UNFIT_EXAMPLES_HELICALVALLEY_HPP_
#define UNFIT_EXAMPLES_HELICALVALLEY_HPP_

#include <cmath>
#include <vector>
#include "GenericCostFunction.hpp"

namespace Unfit
{
namespace Examples
{
/**
 * \brief Implements the Helical Valley function
 *
 * Here the goal is to find a parameter set that minimises the Helical Valley
 * function. This function is defined as:
 *
 *   100{[z-10phi(x,y)]^2 + (sqrt(x^2+y^2) -1)^2} + z^2
 *   phi = arctan (y/x)/(2*PI) for x>0
 *
 * The goal is to find the values of x, y & z that give a minimum cost. In terms
 * of the model, x = param[0], y = param[1] and z = param[2].
 */
class HelicalValley : public GenericCostFunction
{
 public:
  /**
   * Calculate the residuals for the Helical Valley function.
   *
   * Intended use :
   *   residuals = cost_func(param)
   *
   * \param param A vector containing the current estimates of the parameters
   * \return A vector containing the residuals
   */
  std::vector<double> operator()(const std::vector<double> &param)
  {
    const double math_pi = acos(-1.0);
    double theta = atan(param[1] / param[0]) / (2.0 * math_pi);
    if (param[0] <= 0.0) theta += 0.5;
    std::vector<double> residuals(3);
    residuals[0] = 10.0*(param[2] - 10.0*theta);
    residuals[1] = sqrt(param[0]*param[0] + param[1]*param[1]) - 1.0;
    residuals[2] = param[2];
    return residuals;
  }
};

}  // namespace Examples
}  // namespace Unfit

#endif
