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
#ifndef UNFIT_EXAMPLES_ODEMODELWITHBAILOUT_HPP_
#define UNFIT_EXAMPLES_ODEMODELWITHBAILOUT_HPP_

#include <vector>
#include "GenericModel.hpp"

namespace Unfit
{
namespace Examples
{
/**
 * This is model is the same as Ode3DModel.hpp, and you should look at that
 * model for more detailed documentation on the ODE solution. This example model
 * showcases Unfit's "bailout" function. If you have a model that takes a long
 * time to evaluate, it may be useful to check that the current parameter set
 * is working as you go. If you find the parameter set generates erroneous
 * results there is a fast way to discard it in Unfit. All you have to do is
 * set the first entry in the model vector to be NaN or Infinite and return the
 * model vector, no matter how incomplete it is. The current parameter set will
 * be discarded. Internally, Unfit checks the first entry that the model returns
 * (or that the cost function returns if you write a cost function directly)
 * before proceeding with the optimisation algorithm. This is particularly
 * useful for e.g. ODE solutions with many steps. If you can detect something
 * is wrong early on and bail out, your optimisation can be much faster. It is
 * also useful to put in a bail out condition if your model seems to produce
 * NaNs or Infs from time to time as performing calculations on these can be
 * very slow.
 *
 */
class OdeModelWithBailout : public GenericModel
{
 public:
  /**
   * This method takes in a vector of model parameters (c), and a vector of
   * times (t). The time step used in the integration is calculated from pairs
   * of adjacent entries in the t vector. Note that the third of the three
   * parameters represents the initial condition for this initial value problem.
   * Here we arbitrarily assume that the model must produce positive results
   * for all time, and thus any negative values are an error and trigger a bail
   * out from the model calculations.
   *
   * \param c A vector of model parameters/constants
   * \param t A vector of time values (stored in t[0])
   * \return A vector containing the model solution at each time
   */
  std::vector<double> operator()(const std::vector<double> &c,
      const std::vector<std::vector<double>> &t)
  {
    auto model = t[0];  // Create a vector to store the model results
    model[0] = c[2];    // Impose the initial condition
    // Loop over time and solve via forward Euler
    for (auto i = 1u; i < model.size(); ++i) {
      // Get the appropriate time step from the t vector
      const auto dt = t[0][i] - t[0][i-1];
      model[i] = model[i-1] + dt * (c[0] * model[i-1] + c[1]);
      // ***BAILOUT ***
      // Let's assume the model must be positive at all times, if we find a
      // negative in the model we want a fast exit
      if (model[i] < 0.0) {
        // Set the first model result to an nan or infinity...
        model[0] = std::numeric_limits<double>::signaling_NaN();
        // ...then return the model. The rest of the model need not be correct.
        // The nan/inf in the first position will cause this parameter set to be
        // discarded.
        return model;
      }
    }
    return model;
  }
};

}  // namespace Examples
}  // namespace Unfit

#endif  // UNFIT_EXAMPLES_ODEMODELWITHBAILOUT_HPP_
