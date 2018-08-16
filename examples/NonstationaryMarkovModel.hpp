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
#ifndef UNFIT_EXAMPLES_NONSTATIONARYMARKOVMODEL_HPP_
#define UNFIT_EXAMPLES_NONSTATIONARYMARKOVMODEL_HPP_

#include <cmath>
#include <vector>
#include "GenericModel.hpp"

namespace Unfit
{
namespace Examples
{
/**
 * A model to describe eta (the opening rate) for a non-stationary Markov ion
 * channel model (sodium channel) in terms of both voltage and time. The
 * equation looks like:
 *
 *   eta = c0 exp(vm/c1) exp(−(t−c2)^2 / 2(c3)^2)
 *
 * In terms of our nonclemature, this is a 3D model as we have eta = f(vm,t),
 * with four parameters (c0, c1, c2, c3).
 */
class NonstationaryMarkovModel : public GenericModel
{
 public:
  /**
   * This method takes in a vector of model parameters (c), and a vector of
   * membrane potentials at which we want the model to be evaluated (vm). It
   * then evaluates the model at each voltage to calculate a vector of the
   * rate constant, alpha_n, as a function of voltage.
   *
   * \param c A vector of model parameters/constants
   * \param vm A vector of membrane potentials (stored in vm[0])
   * \return A vector containing eta at each membrane potential
   */
  std::vector<double> operator()(const std::vector<double> &c,
      const std::vector<std::vector<double>> &data)
  {
    auto model = data[0];
    for (auto i = 0u; i < model.size(); ++i) {
      const auto t = data[0][i];   // Time data is stored in data[0]
      const auto vm = data[1][i];  // Voltage data is stored in data[1]
      model[i] = c[0]*exp(vm/c[1])*exp(-(t-c[2])*(t-c[2])/(2.0*c[3]*c[3]));
    }
    return model;
  }
};

}  // namespace Examples
}  // namespace Unfit

#endif  // UNFIT_EXAMPLES_NONSTATIONARYMARKOVMODEL_HPP_
