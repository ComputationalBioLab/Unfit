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
#ifndef UNFIT_EXAMPLES_HELICALVALLEYCOSTFUNCTION_HPP_
#define UNFIT_EXAMPLES_HELICALVALLEYCOSTFUNCTION_HPP_

#include <cmath>
#include <vector>
#include "GenericCostFunction.hpp"

namespace Unfit
{
namespace Examples
{
/**
 * This class implements the three parameter Helical Valley function, which is
 * sometimes known as the Fletcher-Powell Helical Valley function, and is
 * commonly used to test optimisation methods. The function looks like:
 *
 *  f(x) = 100*{[c2 - 10*theta]^2 + [sqrt(c0^2 + c1^2) - 1]^2} + c2^2
 *
 *    where theta = arctan(c2/c1) / (2*pi) for c0 >= 0
 *                = 0.5 + arctan(c2/c1) / (2*pi) for c0 < 0
 *
 * Here the goal is to find the values of the c0, c1, and c2 parameters that
 * minimise this function. The correct answer is c0 = 1, c1 = 0, and c2 = 0.
 */
class HelicalValleyCostFunction : public GenericCostFunction
{
 public:
  /**
   * Evaluates the Helical Valley function, given an estimate of the three
   * parameters, c0, c1, and c2, which are passed in.
   *
   * \param c A vector containing the values of the three parameters
   * \return A vector of length one containing the residual/cost
   */
  std::vector<double> operator()(const std::vector<double> &c)
  {
    const double math_pi = std::acos(-1.0);
    double theta = std::atan(c[1]/c[0]) / (2.0 * math_pi);
    if (c[0] <= 0.0) theta += 0.5;
    const double term_1 = c[2] - 10.0*theta;
    const double term_2 = std::sqrt(c[0]*c[0] + c[1]*c[1]) - 1.0;
    std::vector<double> residuals(1);
    residuals[0] = 100.0*(term_1*term_1 + term_2*term_2) + c[2]*c[2];
    return residuals;
  }
};

}  // namespace Examples
}  // namespace Unfit

#endif  // UNFIT_EXAMPLES_HELICALVALLEYCOSTFUNCTION_HPP_
