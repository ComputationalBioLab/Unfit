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
#ifndef UNFIT_EXAMPLES_CARDIACALPHANMODEL_HPP_
#define UNFIT_EXAMPLES_CARDIACALPHANMODEL_HPP_

#include <cmath>
#include <vector>
#include "GenericNDModel.hpp"

namespace Unfit
{
namespace Examples
{
/**
 * A model to describe the alpha_n rate constant that arises in a Hodgkin-Huxley
 * style description of a cardiac potassium channel. The chosen model is a
 * logistic sigmoid function of the membrane potential (vm):
 *
 *   alpha_n = c0 / [1 + exp(c1*vm + c2)]
 */
class CardiacAlphaNModel : public GenericNDModel
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
   * \return A vector containing alpha_n at each membrane potential
   */
  std::vector<double> operator()(const std::vector<double> &c,
      const std::vector<std::vector<double>> &vm)
  {
    auto model = vm[0];
    for (auto &m : model) {
      m = c[0] / (1.0 + exp(c[1]*m + c[2]));
    }
    return model;
  }
};

}  // namespace Examples
}  // namespace Unfit

#endif  // UNFIT_EXAMPLES_CARDIACALPHAN_MODEL_HPP_
