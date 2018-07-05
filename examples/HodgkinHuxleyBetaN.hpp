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
#ifndef UNFIT_EXAMPLES_HODGKINHUXLEYBETAN_HPP_
#define UNFIT_EXAMPLES_HODGKINHUXLEYBETAN_HPP_

#include <cmath>
#include <vector>
#include "GenericCostFunction.hpp"

namespace Unfit
{
namespace Examples
{
/**
 * \brief Fit the voltage dependence of a Hodgkin Huxley rate constant
 *
 * Here the goal is to fit the beta_n rate constant (potassium channel
 * inactivation) using data digitized from Hodgkin & Huxley's landmark 1952
 * paper, Figure 4 (PMCID: PMC1392413) as a function of the cell membrane
 * potential (vm). The equation that is used to fit the data is:
 *
 *   beta_n = A*exp(-vm/B)
 *
 * The goal is to find the values of A & B that best fit the data. In terms
 * of the model, A = param[0] and B = param[1].
 */
class HodgkinHuxleyBetaN : public GenericCostFunction
{
 public:
  /**
   * Create the cost function. Here the experimental data must be passed in,
   * and cannot be changed (if you want to, just create another cost function
   * object). Here the experimental data is two vectors, membrane potential (vm)
   * and the rate constant beta_n.
   *
   * Intended use :
   *   HodgkinHuxleyBetaN cost_func(vm, beta_n);
   *
   * \param vm A vector of membrane potentials
   * \param beta_n A vector of experimental beta_n data
   */
  HodgkinHuxleyBetaN(const std::vector<double> &vm,
      const std::vector<double> beta_n)
    : vm_ {vm},
      beta_n_ {beta_n}
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
    auto residuals = beta_n_;
    for (unsigned i = 0; i < residuals.size(); ++i) {
      residuals[i] -= param[0] * exp(-vm_[i] / param[1]);
    }
    return residuals;
  }
 private:
  /**A vector to store the experimental membrane potential, vm*/
  const std::vector<double> vm_;
  /**A vector to store the experimental rate constant, beta_n*/
  const std::vector<double> beta_n_;
};

}  // namespace Examples
}  // namespace Unfit

#endif
