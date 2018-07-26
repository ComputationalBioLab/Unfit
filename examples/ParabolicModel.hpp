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
#ifndef UNFIT_EXAMPLES_PARABOLICMODEL_HPP_
#define UNFIT_EXAMPLES_PARABOLICMODEL_HPP_

#include <vector>
#include "GenericNDModel.hpp"

namespace Unfit
{
namespace Examples
{
/**
 * A model of a parabola, which has three parameters. The function is given by:
 *
 *   f(x) = c0*x*x + c1*x + c2
 *
 * In terms of terminology, we have a three parameter two-dimensional model. The
 * three parameters come from us having three c's (c0, c1, c2). The
 * two-dimensional part comes from us having a function of two variables, i.e.,
 * we could write this function in the form y = f(x), where x & y are our two
 * variables.
 */
class ParabolicModel : public GenericNDModel
{
 public:
  /**
   * Given estimates of the parameters c0, c1, and c2, plus a vector containing
   * the x-coordinates/values of interest, this method will evaluate the model
   * f(x) = c0*x*x + c1*x + c2 at each x, and return a vector containing the
   * a value for f(x) for each x.
   *
   * In terms of the code, our values for c0, c1 & c2 are stored in c[0], c[1],
   * and c[2], respectively.The x coordinates at which we want to evaluate the
   * model at are stored in a vector in x[0] (x is a vector of vectors), meaning
   * the values can be accessed at x[0][0], x[0][1], x[0][...] etc.
   *
   * \param c A vector of model parameters/constants
   * \param x A vector of x values/coordinates (stored in x[0])
   * \return A vector the model f(x) evaluated at each x
   */
  std::vector<double> operator()(const std::vector<double> &c,
      const std::vector<std::vector<double>> &x)
  {
    // "model" = f(x_i) needs to be the same length as x. Here we making model
    // the same as our x data which guarantees this. We can then use a range
    // based for loop and not have to deal with any indices.
    //
    // FYI: The equivalent code without using a range-based for would look
    // something like this:
    //
    //    std::vector<double> model(x[0].size());
    //    for (auto i = 0u; i < x.size(); ++i) {
    //      model[i] = c[0]*x[0][i]*x[0][i] + c[1]*x[0][i] + c[2];
    //    }
    //
    auto model = x[0];
    for (auto &mx : model) {
      mx = c[0]*mx*mx + c[1]*mx + c[2];
    }
    return model;
  }
};

}  // namespace Examples
}  // namespace Unfit

#endif
