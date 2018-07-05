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
#ifndef UNFIT_EXAMPLES_NONSTATIONARYMARKOV_HPP_
#define UNFIT_EXAMPLES_NONSTATIONARYMARKOV_HPP_

#include <cmath>
#include <vector>
#include "GenericCostFunction.hpp"

namespace Unfit
{
namespace Examples
{
/**
 * \brief An example data fitting problem with four parameters
 *
 * Here the goal is to find a function that describes eta (the opening rate)
 * from a non-stationary Markov ion channel model in terms of both voltage and
 * time (sodium channel). The equation to fit is:
 *
 *   eta = A1 exp(vm/A2) exp(−(t−A3)^2 / 2(A4)^2)
 *
 * The goal is to find the values of A1, A2, A3 & A4 that best fit the data. In
 * terms of the model, A1 = param[0], A2 = param[1], A3 = param[2] and
 * A4 = param[2].
 */
class NonStationaryMarkov : public GenericCostFunction
{
 public:
  /**
   * Create the cost function. Here the experimental data must be passed in,
   * and cannot be changed (if you want to, just create another cost function
   * object). Here the experimental data is a vector of data at eight clamping
   * voltages (a vector of vectors).
   *
   * markov_data[0] stores the time step
   * markov_data[1] stores the eta values at voltage value, -60mV
   * markov_data[2] stores the eta values at voltage value, -50mV
   * markov_data[3] stores the eta values at voltage value, -40mV
   * markov_data[4] stores the eta values at voltage value, -30mV
   * markov_data[5] stores the eta values at voltage value, -20mV
   * markov_data[6] stores the eta values at voltage value, -10mV
   * markov_data[7] stores the eta values at voltage value, 0mV
   *
   * Intended use :
   *   NonStationaryMarkov cost_func(markov_data);
   *
   * \param markov_data A vector of vectors of experimental data (eta)
   */
  NonStationaryMarkov(const std::vector<std::vector<double>> &markov_data)
    : markov_data_ {markov_data}
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
    const std::size_t num_voltages = markov_data_.size();
    const std::size_t num_data_points = markov_data_[0].size();
    std::vector<double> residuals(num_voltages*num_data_points);
    double vm = -60.0;
    // Note that the zero index is time
    for (auto j = 1u; j < num_voltages; ++j) {
      for (unsigned i = 0u; i < num_data_points; ++i) {
        const auto t = markov_data_[0][i];  // time
        residuals[j*num_data_points + i] = markov_data_[j][i]      // experiment
            - param[0] * exp(vm / param[1]) * exp(-(t - param[2]) *
            (t - param[2]) / (2.0 * param[3] * param[3]));         // - model
       }
       vm += 10.0;
    }
    return residuals;
  }
 private:
  /** A vector of vectors to store the experimental data*/
  const std::vector<std::vector<double>> markov_data_;
};

}  // namespace Examples
}  // namespace Unfit

#endif
