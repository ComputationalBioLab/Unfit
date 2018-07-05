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
#ifndef UNFIT_EXAMPLES_WOODS_HPP_
#define UNFIT_EXAMPLES_WOODS_HPP_

#include <cmath>
#include <vector>
#include "GenericCostFunction.hpp"

namespace Unfit
{
namespace Examples
{
/**
 * \brief Implements the Woods test function
 *
 * Here the goal is to find a parameter set that minimises the Woods test
 * function (four dimensions). This function is defined as:
 *
 *   residuals[0] = 10.0*(B - A*A);
 *   residuals[1] = 1.0 - A;
 *   residuals[2] = sqrt(90.0)*(D - C*C);
 *   residuals[3] = 1.0 - C;
 *   residuals[4] = sqrt(10.0)*(B + D - 2.0);
 *   residuals[5] = (B - D) / sqrt(10.0);
 *
 * The goal is to find the values of A, B, C & D that give a minimum cost. In
 * terms of the model, A = param[0], B = param[1], C = param[2], and
 * D = param[3].
 */
class Woods : public GenericCostFunction
{
 public:
  /**
   * Calculate the residuals for the Woods function.
   *
   * Intended use :
   *   residuals = cost_func(param)
   *
   * \param param A vector containing the current estimates of the parameters
   * \return A vector containing the residuals
   */
  std::vector<double> operator()(const std::vector<double> &param)
  {
    std::vector<double> residuals(6);
    residuals[0] = 10.0*(param[1] - param[0]*param[0]);
    residuals[1] = 1.0 - param[0];
    residuals[2] = sqrt(90.0)*(param[3] - param[2]*param[2]);
    residuals[3] = 1.0 - param[2];
    residuals[4] = sqrt(10.0)*(param[1] + param[3] - 2.0);
    residuals[5] = (param[1] - param[3]) / sqrt(10.0);
    return residuals;
  }
};

}  // namespace Examples
}  // namespace Unfit

#endif
