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
#ifndef UNFIT_EXAMPLES_SIMPLEPARABOLICCOSTFUNCTION_HPP_
#define UNFIT_EXAMPLES_SIMPLEPARABOLICCOSTFUNCTION_HPP_

#include <cmath>
#include <vector>
#include "GenericCostFunction.hpp"

namespace Unfit
{
namespace Examples
{
/**
 * This class implements a simple one parameter parabolic function. The function
 * looks like:
 *
 *  f(x) = |c0*c0|
 *
 * Here the goal is to find the value of the c0 that minimises this function.
 * The correct answer is, of course, c0 = 0.
 */
class SimpleParabolicCostFunction : public GenericCostFunction
{
 public:
  /**
   * Evaluates a simple one parameter parabolic function, given an estimate of
   * the parameter c0, which is passed in.
   *
   * \param c A vector containing the value of the parameter
   * \return A vector of length one containing the residual/cost
   */
  std::vector<double> operator()(const std::vector<double> &c)
  {
    return {std::fabs(c[0]*c[0])};
    // The one-liner above is equivalent to:
    //   std::vector<double> residuals(1);
    //   residuals[0] = std::fabs(c[0]*c[0]);
    //   return residuals;
    // which may be more use if you are writing your own.
  }
};

}  // namespace Examples
}  // namespace Unfit

#endif  // UNFIT_EXAMPLES_SIMPLEPARABOLICCOSTFUNCTION_HPP_
