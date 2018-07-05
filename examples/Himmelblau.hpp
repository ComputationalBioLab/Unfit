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
#ifndef UNFIT_EXAMPLES_HIMMELBLAU_HPP_
#define UNFIT_EXAMPLES_HIMMELBLAU_HPP_

#include <cmath>
#include <vector>
#include "GenericCostFunction.hpp"

namespace Unfit
{
namespace Examples
{
/**
 * The Himmelblau Function is defined as (x^2 + y -11)^2 + (x + y^2 - 7)^2
 *
 * Number of dimensions = 2
 *
 * The global minimum is: 0 at (3, 2), (-2.805118, 3.131312),
 * (-3.779310, -3.283186) and (3.584428, -1.848126)T
 * Initial guess: (4, 2), (-4, 2), (-4, -2), (4, -2)
 *
 */
class Himmelblau : public GenericCostFunction
{
 public:
  /**
   * We overload the operator as required in GenericCostFunction to
   * calculate the difference between each experimental
   * value and its corresponding value predicted by the model
   * equation i.e. residual
   *
   * Intended use :
   *   residuals = cost_func(const std::vector<double> &param);
   *
   * \param variable vector containing the current estimates of the variables
   * \return result vector containing the residuals
   */
  std::vector<double> operator()(const std::vector<double> &variable)
  {
    std::vector<double> result(2, 0.0);
    result[0] = (variable[0]*variable[0] + variable[1] - 11)*(variable[0]*
        variable[0] + variable[1] - 11);
    result[1] = (variable[0] + variable[1]*variable[1] - 7)*(variable[0] +
        variable[1]*variable[1] - 7);
    return result;
  }
};

}  // namespace Examples
}  // namespace Unfit

#endif
