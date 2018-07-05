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
#ifndef UNFIT_UNITTESTS_LEVENBERGMARQUARDTTESTFUNCTIONS_HPP_
#define UNFIT_UNITTESTS_LEVENBERGMARQUARDTTESTFUNCTIONS_HPP_

#include <limits>
#include <vector>
#include "GenericCostFunction.hpp"

namespace Unfit
{
namespace Examples
{
/**
* \brief Equation1 is defined as x^2 + 2x + 1
*
* Number of dimensions = 1
* The global minimum is: 0 at (-1)T
* Initial guess: (10)
*
*/
class Equation1 : public GenericCostFunction
{
 public:
  /**
   * We overload the operator as is required in GenericCostFunction to
   * calculate the cost of the function.
   *
   * Behaviour:
   *   cost = x^2 + 2x + 1
   *
   * Intended use :
   *   Equation1 Func;
   *   cost = Func(const std::vector<double> variable);
   *
   * Parameters:
   *   \param variable (input) vector containing coordinates of x
   *   \return result i.e. cost
   */
  std::vector<double> operator()(const std::vector<double> &variable)
  {
    std::vector<double> result {0};
    result[0] = variable[0]*variable[0] + 2*variable[0] + 1;
    return result;
  }
};

/**
 * \brief Equation2 is defined as x^4 + 2x^2*y^2 + y^4
 *
 * Number of dimensions = 2
 * The global minimum is: 0 at (0, 0)T
 * Initial guess: (10, 10)
 *
 */
class Equation2 : public GenericCostFunction
{
 public:
  /**
   * We overload the operator as is required in GenericCostFunction to
   * calculate the cost of the function.
   *
   * Behaviour:
   *   cost = x^4 + 2x^2*y^2 + y^4
   *
   * Intended use :
   *   Equation2 Func;
   *   cost = Func(const std::vector<double> variable);
   *
   * Parameters:
   *   \param variable (input) vector containing coordinates of x and y
   *   \return result i.e. cost
   */
  std::vector<double> operator()(const std::vector<double> &variable)
  {
    std::vector<double> result(3, 0.0);
    result[0] = variable[0]*variable[0]*variable[0]*variable[0];
    result[1] = 2*variable[0]*variable[0] * variable[1]*variable[1];
    result[2] = variable[1]*variable[1]*variable[1]*variable[1];
    return result;
  }
};

/**
 * \brief Equation3 is defined as (x - 1)^2 + (y - 1)^2
 *
 * Number of dimensions = 2
 * The global minimum is: 0 at (1, 1)T
 * Initial guess: (10, 10)
 *
 */
class Equation3 : public GenericCostFunction
{
 public:
  /**
   * We overload the operator as is required in GenericCostFunction to calculate
   * the cost of the function.
   *
   * Behaviour:
   *   cost = (x - 1)^2 + (y - 1)^2
   *
   * Intended use :
   *   Equation3 Func;
   *   cost = Func(const std::vector<double> variable);
   *
   * Parameters:
   *   \param variable (input) vector containing coordinates of x and y
   *   \return result i.e. cost
   */
  std::vector<double> operator()(const std::vector<double> &variable)
  {
    std::vector<double> result(2, 0.0);
    result[0] = (variable[0]-1)*(variable[0]-1);
    result[1] = (variable[1]-1)*(variable[1]-1);
    return result;
  }
};

/**
 * \brief Equation4 is defined as 4x^4
 *
 * Number of dimensions = 1
 * The global minimum is: 0 at (0)T
 * Initial guess: (1)
 *
 */
class Equation4 : public GenericCostFunction
{
 public:
  /**
   * We overload the operator as is required in GenericCostFunction to calculate
   * the cost of the function.
   *
   * Behaviour:
   *   cost = 4x^4
   *
   * Intended use :
   *   Equation4 Func;
   *   cost = Func(const std::vector<double> variable);
   *
   * Parameters:
   *   \param variable (input) vector containing coordinates of x
   *   \return result i.e. cost
   */
  std::vector<double> operator()(const std::vector<double> &variable)
  {
    std::vector<double> result {0};
    result[0] = 4*variable[0]*variable[0]*variable[0]*variable[0];
    return result;
  }
};

/**
* This test function is a simple parabola, y=x*x, where we want to find x=0.
* The catch is that this one has a hole filled with NaNs around the solution.
*/
class ParabolaWithHole : public GenericCostFunction
{
 public:
  /**
   * We overload the operator as is required in GenericCostFunction to
   * calculate the cost of the function.
   *
   * \param variable the current estimate of x
   * \return a vector of residuals (length one in this case)
   */
  std::vector<double> operator()(const std::vector<double> &variable)
  {
    std::vector<double> residuals {0};
    if (fabs(variable[0]) < 0.1) {
      residuals[0] = std::numeric_limits<double>::quiet_NaN();
    }
    else {
      residuals[0] = variable[0]*variable[0];
    }
    return residuals;
  }
};

/* The equations below are commented out because LM cannot find minimum that
    is negative
 */

// /**
// * \brief Equation5 is defined as x^2 + 2y + 1
// *
// * Number of dimensions = 2
// * The global minimum is: 0 at (-1)T
// * Initial guess: ???
// *
// */
// class Equation5 : public GenericCostFunction
// {
// public:
//  /**
//   * We overload the operator as is required in GenericCostFunction to
//   * calculate the cost of the function.
//   *
//   * Behaviour:
//   *   cost = x^2 + 2y + 1
//   *
//   * Intended use :
//   *   Equation5 Func;
//   *   cost = Func(const std::vector<double> variable);
//   *
//   * Parameters:
//   *   \param variable (input) vector containing coordinates of x and y
//   *   \return cost
//   */
//  std::vector<double> operator()(const std::vector<double> &variable)
//  {
//    std::vector<double> result {0};
//    result[0] = variable[0]*variable[0] + 2*variable[1] + 1;
//    return result;
//    // doesnt work as the cost here is -ve inf for x or y
//    // since LM will square it, it will find min at (0.0, 0.0) instead
//  }
// };
//
// /**
// * \brief Equation6 is defined as x^3 + y^3
// *
// * Number of dimensions = 2
// * The global minimum is: negative infinity
// * Initial guess: ???
// *
// */
// class Equation6 : public GenericCostFunction
// {
// public:
//  /**
//   * We overload the operator as is required in GenericCostFunction to
//   * calculate the cost of the function.
//   *
//   * Behaviour:
//   *   cost = x^3 + y^3
//   *
//   * Intended use :
//   *   Equation6 Func;
//   *   cost = Func(const std::vector<double> variable);
//   *
//   * Parameters:
//   *   \param variable (input) vector containing coordinates of x and y
//   *   \return cost
//   */
//  std::vector <double> operator()(const std::vector <double> &variable)
//  {
//    std::vector <double> result {0};
//    result[0] = variable[0]*variable[0]*variable[0] + variable[1]*variable[1]*
//        variable[1];
//    return result;
//  }
// };
//
// /**
// * \brief Equation7 is defined as x^3 + y^3 - 2
// *
// * Number of dimensions = 2
// * The global minimum is: negative infinity
// * Initial guess: ???
// *
// */
// class Equation7 : public GenericCostFunction
// {
// public:
//  /**
//   * We overload the operator as is required in GenericCostFunction to
//   * calculate the cost of the function.
//   *
//   * Behaviour:
//   *   cost = x^3 + y^3 - 2
//   *
//   * Intended use :
//   *   Equation7 Func;
//   *   cost = Func(const std::vector<double> variable);
//   *
//   * Parameters:
//   *   \param variable (input) vector containing coordinates of x and y
//   *   \return cost
//   */
//  std::vector <double> operator()(const std::vector <double> &variable)
//  {
//    std::vector <double> result {0};
//    result[0] = variable[0]*variable[0]*variable[0] + variable[1]*variable[1]*
//        variable[1] - 2;
//    return result;
//  }
// };
//
// /**
// * \brief Equation8 is defined as 2*x^2 + 3*y + z
// *
// * Number of dimensions = 3
// * The global minimum is: 0 at (0)T
// * Initial guess: ???
// *
// */
// class Equation8 : public GenericCostFunction
// {
// public:
//  /**
//   * We overload the operator as is required in GenericCostFunction to
//   * calculate the cost of the function.
//   *
//   * Behaviour:
//   *   cost = 2*x^2 + 3*y + z
//   *
//   * Intended use :
//   *   Equation8 Func;
//   *   cost = Func(const std::vector<double> variable);
//   *
//   * Parameters:
//   *   \param variable (input) vector containing coordinates of x, y and z
//   *   \return cost
//   */
//  std::vector<double> operator()(const std::vector<double> &variable)
//  {
//    std::vector<double> result(3, 0.0);
//    result[0] = 2*variable[0]*variable[0];
//    result[1] = 3*variable[1];
//    result[2] = variable[2];
//    return result;
//  }
// };
//
// /**
// * \brief Equation9 is defined as x*y
// *
// * Number of dimensions = 2
// * The global minimum is: ???
// * Initial guess: (10, 10)
// *
// */
// class Equation9 : public GenericCostFunction
// {
// public:
//  /**
//   * We overload the operator as is required in GenericCostFunction to
//   * calculate the cost of the function.
//   *
//   * Behaviour:
//   *   cost = x*y
//   *
//   * Intended use :
//   *   Equation9 Func;
//   *   cost = Func(const std::vector<double> variable);
//   *
//   * Parameters:
//   *   \param variable (input) vector containing coordinates of x and y
//   *   \return cost
//   */
//  std::vector<double> operator()(const std::vector<double> &variable)
//  {
//    std::vector<double> result(2, 0.0);
//    result[0] = variable[0];
//    result[1] = variable[1];
//    return result;
//    // doesnt work as the cost here is -ve inf for x or y
//    // since LM will square it, it will find min at (0.0, 0.0) instead
//  }
// };

}  // namespace Examples
}  // namespace Unfit

#endif
