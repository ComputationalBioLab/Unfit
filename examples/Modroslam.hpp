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
#ifndef UNFIT_EXAMPLES_MODROSLAM_HPP_
#define UNFIT_EXAMPLES_MODROSLAM_HPP_

#include <cmath>
#include <vector>
#include "GenericCostFunction.hpp"

namespace Unfit
{
namespace Examples
{
/**
 * \brief Implements a Modified Rosenbrock function
 *
 * Here the goal is to find a parameter set that minimises the Modified
 * Rosenbrock function. This function is defined as:
 *
 *   residuals[i] = 10.0*(B - A*A);
 *   residuals[i+1] = 1.0 - A;
 *   residuals[i+2] = 100.0;
 *
 * The goal is to find the values of A & B that give a minimum cost. In terms
 * of the model, A = param[0] and B = param[1].
 */
class Modroslam : public GenericCostFunction
{
 public:
  /**
   * Create the cost function. Here the number of observations must be passed
   * in.
   *
   * Intended use :
   *   Modroslam cost_func(n);
   *
   * \param n The number of observations (must be divisible by 3. If it is not,
   *     the code will truncate back to the nearest multiple of 3).
   */
  Modroslam(unsigned n)
    : n_ {(n < 3) ? 3 : (n - n%3)}
  {}

  /**
   * Calculate the residuals for the Modified Rosenbrock function.
   *
   * Intended use :
   *   residuals = cost_func(param)
   *
   * \param param A vector containing the current estimates of the parameters
   * \return A vector containing the residuals
   */
  std::vector<double> operator()(const std::vector<double> &param)
  {
    std::vector<double> residuals(n_);
    for (auto i = 0u; i < n_; i += 3) {
      residuals[i] = 10.0*(param[1] - param[0]*param[0]);
      residuals[i+1] = 1.0 - param[0];
      residuals[i+2] = 100.0;
    }
    return residuals;
  }
 private:
  /** A variable to store number of observations (should be a multiple of 3)*/
  const unsigned n_;
};

}  // namespace Examples
}  // namespace Unfit

#endif
