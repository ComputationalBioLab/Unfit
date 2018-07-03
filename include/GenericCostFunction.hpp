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
#ifndef UNFIT_INCLUDE_GENERICCOSTFUNCTION_HPP_
#define UNFIT_INCLUDE_GENERICCOSTFUNCTION_HPP_

#include <vector>

namespace Unfit
{
/**
 * An interface class to ensure all functions coming in to our optimizer are
 * consistent. It is an interface class because it contains only pure virtual
 * methods. This means that any methods defined here must be implemented in
 * a derived class. It also means that this class cannot be instantiated
 * directly.
 */
class GenericCostFunction
{
 public:
  /**
   * As we are deriving from this class, the destructor should be virtual. In
   * this case (at present) we have nothing that will not be deleted when the
   * class goes out of scope so an empty destructor method is fine.
   */
  virtual ~GenericCostFunction() {}

  /**
   * All cost functions must implement this method which overloads operator().
   * Just look at any of the bundled examples to see how to do it. This
   * interface is designed to take in a vector containing the current estimates
   * of the unknown model parameters (those we are fitting) and returns a vector
   * of residuals which is the (signed) linear distance between the model and
   * the data. In other words, it should calculate r[] = data[] - model[].
   *
   * Parameters:
   *   \param x A vector containing the current estimates of the unknown model
   *          parameters
   *   \return A vector of residuals, r[] = data[] - model[].
   */
  virtual std::vector<double> operator()(const std::vector<double> &x) = 0;
};

}  // namespace Unfit

#endif
