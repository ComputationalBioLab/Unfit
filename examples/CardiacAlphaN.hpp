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
#ifndef UNFIT_EXAMPLES_CARDIACALPHAN_HPP_
#define UNFIT_EXAMPLES_CARDIACALPHAN_HPP_

#include <cmath>
#include <vector>
#include "GenericCostFunction.hpp"

namespace Unfit
{
namespace Examples
{
/**
 * \brief Fit the voltage dependence of a cardiac ion channel rate constant
 *
 * Here the goal is to fit the alpha_n rate constant as a function of the
 * cell membrane potential (vm), given a set of experimental data (vm and
 * alpha_n). This channel is a modified slow potassium channel from a cardiac
 * muscle cell. The equation that has been chosen to fit the data is:
 *
 *   alpha_n = A / (1 + exp(B*vm + C))
 *
 * The goal is to find the values of A, B & C that best fit the data. In terms
 * of the model, A = param[0], B = param[1] and C = param[2].
 */
class CardiacAlphaN : public GenericCostFunction
{
 public:
  /**
   * Create the cost function. Here the experimental data must be passed in,
   * and cannot be changed (if you want to, just create another cost function
   * object). Here the experimental data is two vectors, membrane potential (vm)
   * and the rate constant alpha_m.
   *
   * Intended use :
   *   CardiacAlphaN cost_func(vm, alpha_n);
   *
   * \param vm A vector of membrane potentials
   * \param alpha_n A vector of experimental alpha_n data
   */
  CardiacAlphaN(const std::vector<double> &vm,
      const std::vector<double> &alpha_n)
    : vm_ {vm},
      alpha_n_ {alpha_n}
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
    auto residuals = alpha_n_;
    for (unsigned i = 0; i < residuals.size(); ++i) {
      residuals[i] -= param[0] / (1.0 + exp(param[1]*vm_[i] + param[2]));
    }
    return residuals;
  }
 private:
  /** A vector to store the experimental membrane potential, vm*/
  const std::vector<double> vm_;
  /** A vector to store the experimental rate constant, alpha_n*/
  const std::vector<double> alpha_n_;
};

}  // namespace Examples
}  // namespace Unfit

#endif
