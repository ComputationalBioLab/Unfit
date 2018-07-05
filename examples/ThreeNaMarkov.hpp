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
#ifndef UNFIT_EXAMPLES_THREENAMARKOV_HPP_
#define UNFIT_EXAMPLES_THREENAMARKOV_HPP_

#include <vector>
#include "GenericCostFunction.hpp"

namespace Unfit
{
namespace Examples
{
/**
 * \brief Fit the time dependence of a three state hidden Markov model
 *
 * This is an example of fitting a hidden Markov model as can be used to
 * describe an ion channel. There are three states which leads to three ordinary
 * differential equations (ODEs) and four rate constants to fit. The states
 * are labelled C = closed, O = open, and I = inactive. The initial conditions
 * for the states are C = 1.0, O & I = 0.0. The three ODEs are:
 *
 *   dO/dt = kco*C + kio*I - (koc+koi)*O
 *   dC/dt = koc*O - kco*C
 *   dI/dt = koi*O - kio*I
 *
 * The goal is to find the values of the k's (rate constants) that best fit the
 * data, which is in the form of O vs time. In terms of the model,
 * kco = param[0], koc = param[1], koi = param[2], kio = param[3]. These ODEs
 * are integrated with the forward Euler method.
 */
class ThreeNaMarkov : public GenericCostFunction
{
 public:
  /**
   * Create the cost function. Here the experimental data must be passed in,
   * and cannot be changed (if you want to, just create another cost function
   * object). Here the experimental data is a vector containing the value
   * of the open probabilities, and the time step used in the data collection.
   *
   * Intended use :
   *   ThreeNaMarkov cost_func(open_prob, dt);
   *
   * \param open_prob A vector of experimental open probabilities
   * \param dt The experimental time step
   */
  ThreeNaMarkov(const std::vector<double> &open_prob, double dt)
    : open_prob_ {open_prob},
      dt_ {dt}
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
    // Initial conditions
    double O = 0.0;
    double C = 1.0;
    double I = 0.0;
    // Rate constants
    const double kco = param[0];
    const double koc = param[1];
    const double koi = param[2];
    const double kio = param[3];
    // Residual initialization
    auto residuals = open_prob_;
    residuals[0] -= O;

    for (auto i = 1u; i < residuals.size(); ++i) {
      const double deltaC = dt_*(koc*O - kco*C);
      const double deltaI = dt_*(koi*O - kio*I);
      I += deltaI;
      C += deltaC;
      O = 1.0 - C - I;
      // NOTE: equivalent behaviour of "O = 1-C-I" is
      // double deltaO = dt_*(kco*C + kio*I - (koc+koi)*O);
      // O += deltaO;
      residuals[i] -= O;
    }
    return (residuals);
  }
 private:
  /** A vector to store the experimental open probabilities */
  const std::vector<double> open_prob_;
  /** A vector to store the time step, dt */
  const double dt_;
};

}  // namespace Examples
}  // namespace Unfit

#endif
