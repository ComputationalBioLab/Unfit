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
#ifndef UNFIT_UNITTESTS_TESTFUNCTIONS_HPP_
#define UNFIT_UNITTESTS_TESTFUNCTIONS_HPP_

#include <limits>
#include <vector>
#include "GenericCostFunction.hpp"
#include "GenericModel.hpp"

namespace Unfit
{
namespace UnitTests
{
/**
 * This function simply returns the value that was supplied as the residual for
 * each element. The minimum is always (0,0...). The last entry is not included
 * in what is returned as that is used to store the cost.
 */
class SimpleCostFunction : public GenericCostFunction
{
 public:
  /**
   * Calculate and return the cost of this function given the parameter
   * estimates.
   *
   * \param x A vector containing the parameter estimates
   * \return The residuals are the same as the supplied vector
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    auto residuals = x;
    residuals.pop_back();
    return residuals;
  }
};

/**
 * This function acts like the SimpleCostFunction but inserts a nan in the
 * first element. The last entry is not included in what is returned as that is
 * used to store the cost.
 */
class FirstNanCostFunction : public GenericCostFunction
{
 public:
  /**
   * Calculate and return the cost of this function given the parameter
   * estimates.
   *
   * \param x A vector containing the parameter estimates
   * \return The residuals are the same as the supplied vector
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    auto residuals = x;
    residuals.pop_back();
    residuals[0] = std::numeric_limits<double>::signaling_NaN();
    return residuals;
  }
};

/**
 * This function is a simple model that returns the first set of x data
 * but sets the first value of the model to a nan.
 */
class FirstNanCostModel : public GenericModel
{
 public:
  /**
   * Calculate and return the cost of this function given the parameter
   * estimates.
   *
   * \param c A vector containing the parameter estimates
   * \param x A vector of vectors containing the independent data
   * \return The first residual is always nan
   */
  std::vector<double> operator()(const std::vector<double> &c,
      const std::vector<std::vector<double>> &x)
  {
    auto model = x[0];
    model[0] = c[0];  // to avoid a 'c not used' warning
    model[0] = std::numeric_limits<double>::signaling_NaN();
    return model;
  }
};

/**
 * This function acts like the SimpleCostFunction but inserts a nan in the
 * last element. The last entry is not included in what is returned as that is
 * used to store the cost.
 */
class LastNanCostFunction : public GenericCostFunction
{
 public:
  /**
   * Calculate and return the cost of this function given the parameter
   * estimates.
   *
   * \param x A vector containing the parameter estimates
   * \return The residuals are the same as the supplied vector
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    auto residuals = x;
    residuals.pop_back();
    residuals.back() = std::numeric_limits<double>::signaling_NaN();
    return residuals;
  }
};

/**
 * This function acts like the SimpleCostFunction but inserts a inf in the
 * first element. The last entry is not included in what is returned as that is
 * used to store the cost.
 */
class FirstInfCostFunction : public GenericCostFunction
{
 public:
  /**
   * Calculate and return the cost of this function given the parameter
   * estimates.
   *
   * \param x A vector containing the parameter estimates
   * \return The residuals are the same as the supplied vector
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    auto residuals = x;
    residuals.pop_back();
    residuals[0] = std::numeric_limits<double>::infinity();
    return residuals;
  }
};

/**
 * This function acts like the SimpleCostFunction but inserts a inf in the
 * last element. The last entry is not included in what is returned as that is
 * used to store the cost.
 */
class LastInfCostFunction : public GenericCostFunction
{
 public:
  /**
   * Calculate and return the cost of this function given the parameter
   * estimates.
   *
   * \param x A vector containing the parameter estimates
   * \return The residuals are the same as the supplied vector
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    auto residuals = x;
    residuals.pop_back();
    residuals.back() = std::numeric_limits<double>::infinity();
    return residuals;
  }
};

/**
 * This function returns a residual of 1.0 if the coordinate lies between 0 and
 * 1 (inclusive), and returns infinitiy otherwise. The last entry is not
 * included in what is returned as that is used to store the cost.
 */
class TableTopCostFunction : public GenericCostFunction
{
 public:
  /**
   * Calculate and return the cost of this function given the parameter
   * estimates.
   *
   * \param x A vector containing the parameter estimates
   * \return The residuals are the same as the supplied vector
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    auto residuals = x;
    residuals.pop_back();
    for (auto &residual : residuals) {
      if (residual >= 0.0 && residual <= 1.0) {
        residual = 1.0;
      }
      else {
        residual = std::numeric_limits<double>::infinity();
      }
    }
    return residuals;
  }
};

/**
 * This function returns a residual of 1.0 if the coordinate lies outside 0 and
 * 1 (inclusive), and returns infinitiy otherwise. The last entry is not
 * included in what is returned as that is used to store the cost.
 */
class ReverseTableTopCostFunction : public GenericCostFunction
{
 public:
  /**
   * Calculate and return the cost of this function given the parameter
   * estimates.
   *
   * \param x A vector containing the parameter estimates
   * \return The residuals are the same as the supplied vector
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    auto residuals = x;
    residuals.pop_back();
    for (auto &residual : residuals) {
      if (residual > 0.0 && residual < 1.0) {
        residual = std::numeric_limits<double>::infinity();
      }
      else {
        residual = 1.0;
      }
    }
    return residuals;
  }
};

}  // namespace Unittests
}  // namespace Unfit

#endif

