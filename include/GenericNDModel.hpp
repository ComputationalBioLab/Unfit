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
#ifndef GENERICNDMODEL_HPP_
#define GENERICNDMODEL_HPP_

#include <vector>

namespace Unfit
{
/**
 * This class provides an interface to what we call an n-dimensional, or ND
 * model. What is a ND model? The simplest way to think of it is to consider a
 * 2D model which is pretty much anything you could draw on a 2D graph
 * with axes x and y, with a function y = f(x). Of course in your particular
 * application you could have e.g. voltage as a function of time: v = f(t),
 * or stress vs strain, or something completely different, but as long as you
 * have one independent variable (e.g. x) and one independent variable (e.g. y),
 * you have a 2D model. What if I have something more complicated, like stress
 * and strain as a function of time, or stress = f(strain, time)? Well, then you
 * have a 3D model, and we could write it in a general form like y = f(x0, x1).
 * Here we have two independent variables that are needed to calculate our y.
 * You can imagine that we could extend this to 4, or 5 or even N-dimensions.
 * That is where the ND model comes in. In general we can write
 * y = f(x0, x1, x2, ... , xN). Whether you have something like this, or
 * something as simple as y = f(x), your model implementation should derive from
 * this class so we keep a consistent interface to our modelling tools. Note
 * that you can have as many parameters/constants as you like. You should look
 * at our examples for further guidance if needed.
 */
class GenericNDModel
{
 public:
  /**
   * As we are deriving from this class, the destructor should be virtual. The
   * default destructor is fine.
   */
  virtual ~GenericNDModel() = default;

  /**
   * All ND models must be implemented using this interface. For details of what
   * constitutes an ND model please see the GenericNDModel class documentation.
   * Here we want to pass in a vector of parameters (c), and a vector of vectors
   * containing the independent variable values at which we want to evaluate our
   * model (x). x[0] will contain a vector with the values of the first
   * independent variable, x[1] will contain the second, if present, etc. The
   * method should then return a vector containing the results of evaluating
   * the model at each x value, i.e., y = f(x0, x1, ... , xN).
   *
   * \param c A vector of model parameters/constants
   * \param x A vector of vectors containing the independent variable values
   *          at which the model should be evaluated.
   * \return  A vector containing the output from the model when evaluated at
   *          x with parameters c.
   */
  virtual std::vector<double> operator()(const std::vector<double> &c,
      const std::vector<std::vector<double>> &x) = 0;
};

}  // namespace Unfit

#endif  // GENERICNDMODEL_HPP_
