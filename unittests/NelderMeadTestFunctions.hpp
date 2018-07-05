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
#ifndef UNFIT_UNITTESTS_NELDERMEADTESTFUNCTIONS_HPP_
#define UNFIT_UNITTESTS_NELDERMEADTESTFUNCTIONS_HPP_

#include <cmath>
#include <limits>
#include <vector>
#include "GenericCostFunction.hpp"

namespace Unfit
{
namespace UnitTests
{
/**
This file contains a list of functions that was used to unit-test
the Nelder Mead method. These are normal functions, without any data.

Often, the calculated cost is simply the squre root of a function.
This is because the Nelder Mead method calculates the cost
by doing a sum of squared resisuals.

Here is an index of the functions you may find here:

LineToInfinity
InfinityInBetween
Infinity2D
SampleFunction1
SampleFunction2 (unused and commented out)
SampleFunction3
SampleFunction4 (unused and commented out)
SampleFunction5
SampleFunction6
SampleFunction7
SampleFunction8
SampleFunction9
Asymptote2D
FailFindMin
FailFindMin2
PowellSingular
AlwaysInfinite

*/

/**
 * This function is a simple line for x > 3, infinity otherwise.
 * It is used to unit-test the case of expansion to infinity
 */
class LineToInfinity : public GenericCostFunction
{
 public:
  /**
   * We overload the operator as is required in GenericCostFunction to calculate
   * the cost of the function.
   *
   * Behaviour:
   *   cost = abs(x) if x> 3, infinity otherwise
   *
   * Intended use :
   *   LineToInfinity Func;
   *   cost = Func(const std::vector<double> x);
   *
   * NOTE that the returned cost is the sqare root of the evaluation
   *      due to the fact the Nelder Mead class will square the cost.
   *
   * Parameters:
   *   \param x (input) vector containing coordinates of x (1D)
   *   \return cost as a vector
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    std::vector<double> ret {0};
    if (x[0] > 3.0) {
        ret[0] = sqrt(fabs(x[0]));
    }
    else {
        ret[0] = std::numeric_limits<double>::infinity();
    }
    return ret;
  }
};

/**
 * This function is used to unit-test the case of expansion to infinity
 * It returns infinity when x is in between -1 and 1.
 */
class InfinityInBetween : public GenericCostFunction
{
 public:
  /**
   * We overload the operator as is required in GenericCostFunction to calculate
   * the cost of the function.
   *
   * Behaviour:
   *   cost = x if x > 1 and x < -1, infinity otherwise
   *
   * Intended use :
   *   InfinityInBetween Func;
   *   cost = Func(const std::vector<double> x);
   *
   * NOTE sicne this function actually retrns x, from th Nelder-Mead
   * persepctive it behaves like a parabola, since NM squares the cost
   *
   * Parameters:
   *   \param x (input) vector containing coordinates of x (1D)
   *   \return cost as a vector
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    std::vector<double> ret {0};
    if (x[0] > 1.0) {
        ret[0] = x[0];
    }
    else if (x[0] < -1.0) {
        ret[0] = x[0];
    }
    else {
        ret[0] = std::numeric_limits<double>::infinity();
    }
    return ret;
  }
};

/**
 * This function is used to unit-test the case of invalid outside contraction.
 */
class Infinity2D : public GenericCostFunction
{
 public:
  /**
   * Here we are minimising x+y for 0.9<x<2.1 and 0.9<y<2.1. Outside this the
   * cost function is constant except for a hole in the domain at (1.25,0.5).
   * We can use this geometry to first force an outside contract and then force
   * the outside contract to land in the hole and therefore test we get the
   * correct response to an outside contract failure (a shrink).
   *
   * Parameters:
   *   \param x (input) vector containing coordinates of x (2D)
   *   \return cost as a vector (3D)
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    std::vector<double> f_i(3, 0.0);
    if ((x[0] < 1.26 && x[0] > 1.24) && (x[1] < 0.51 && x[1] > 0.49)) {
      f_i[0] = std::numeric_limits<double>::infinity();
      f_i[1] = std::numeric_limits<double>::infinity();
      f_i[2] = std::numeric_limits<double>::infinity();
    }
    else if (x[0] < 0.9 || x[0] > 2.1 || x[1] < 0.9 || x[1] > 2.1) {
      f_i[0] = 3.5;
      f_i[1] = 3.5;
      f_i[2] = 3.5;
    }
    else {
      f_i[0] = x[0] + x[1];
      f_i[1] = x[0] + x[1];
      f_i[2] = x[0] + x[1];
    }
    return f_i;
  }
};

/**
 * The Sample function 1 is defined as
 *
 *   (x^2 - 4x + y^2 - y - xy).
 *
 * Number of dimensions = 2
 * The global minimum is: (3, 2)
 * Initial guess: (0, 0.8)
 *
 *
 * Reference : https://docs.google.com/viewer?a=v&q=cache:qDWbSzyjw98J:math.
 * fullerton.edu/mathews/n2003/neldermead/NelderMeadProof.pdf+&hl=en&gl=sg&pid=
 * bl&srcid=ADGEESgadIne-2Pom8fA4eyCi-wrcuSOu4Qytf8KctcUBggvGpahFMKcco7UhNUCS
 * -1ie2ChoM3OjO5RSgp05mrK3lBtbsXQFQFrMP426SOH8BGxZO4YcdbRBXsGn4HsyO02eKFIgWIr
 * &sig=AHIEtbR1PU77goer0s2S6zH7HKthD0USmQ
 */
class SampleCostFunction1 : public GenericCostFunction
{
 public:
  /**
   * We overload the operator as is required in GenericCostFunction to calculate
   * the cost of the function.
   *
   * Behaviour:
   *   cost = x^2 - 4x + y^2 - y - xy
   *
   * Intended use :
   *   SampleCostFunction1 Func;
   *   cost = Func(const std::vector<double> x);
   *
   * NOTE that the returned cost is the sqare root of the evaluation
   *      due to the fact the Nelder Mead class will square the cost.
   *
   * Parameters:
   *   \param x (input) vector containing coordinates of x and y
   *   \return cost
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    std::vector<double> ret {0};
    ret[0] = sqrt(fabs(x[0]*x[0] - 4*x[0] + x[1]*x[1] - x[1] - x[0]*x[1]));
    return ret;
  }
};

// /**
// * The Sample function 2 is defined as
// *
// *   (y^2 + 4xy +x^3 - 2x).
// *
// * Number of dimensions = 2
// * The global minimum is: (2.896805253, -5.793610507)
// * Initial guess: (0, 0.8)
// *
// * Reference: http://www.wolframalpha.com/input/
// * ?i=%28y^2+%2B+4xy+%2Bx^3+-+2x%29
// */
// class SampleCostFunction2 : public GenericCostFunction
// {
// public:
//  /**
//   * We overload the operator as is required in GenericCostFunction to
//   * calculate the cost of the function.
//   *
//   * Behaviour:
//   *   cost = y^2 + 4xy +x^3 - 2x
//   *
//   * Intended use :
//   *   SampleCostFunction2 Func;
//   *   cost = Func(const std::vector<double> x);
//   *
//   * NOTE that the returned cost is the sqare root of the evaluation
//   *      due to the fact the Nelder Mead class will square the cost.
//   *
//   * Parameter:
//   *   \param x (input) vector containing coordinates of x and y
//   *   \return cost
//   */
//  std::vector<double> operator()(const std::vector<double> &x)
//  {
//    std::vector<double> ret {0};
//    ret[0] = sqrt(fabs(x[1]*x[1] + 4*x[1]*x[0] + x[0]*x[0]*x[0]- 2*x[0]));
//    return ret;
//  }
// };


/**
 * The Sample function 3 is defined as
 *
 *   (sqrt(y^2 + x^2))
 *
 * Number of dimensions = 2
 * The global minimum is: (0, 0)
 * Initial guess: (5, 0)
 *
 * Reference:
 */
class SampleCostFunction3 : public GenericCostFunction
{
 public:
  /**
   * We overload the operator as is required in GenericCostFunction to calculate
   * the cost of the function.
   *
   * Behaviour:
   *   cost = sqrt(y^2 + x^2)
   *
   * Intended use :
   *   SampleCostFunction3 Func;
   *   cost = Func(const std::vector<double> x);
   *
   * NOTE that the returned cost is the sqare root of the evaluation
   *      due to the fact the Nelder Mead class will square the cost.
   *
   * Parameters:
   *   \param x (input) vector containing coordinates of x and y
   *   \return cost
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    std::vector<double> ret {0};
    ret[0] = sqrt(sqrt(x[1]*x[1] + x[0]*x[0]));
    return ret;
  }
};

// /**
// * The Sample function 4 is defined as
// *
// *   (x^2 + (y - 11)^2 + (y^2 - 7 + x)^2).
// *
// * Number of dimensions = 2
// * The global minimum is: (3, 2)
// * Initial guess: (0, 0)
// *
// * Reference:
// */
// class SampleCostFunction4 : public GenericCostFunction
// {
// public:
//  /**
//   * We overload the operator as is required in GenericCostFunction to
//   * calculate the cost of the function.
//   *
//   * Behaviour:
//   *   cost = x^2 + (y - 11)^2 + (y^2 - 7 + x)^2
//   *
//   * Intended use :
//   *   SampleCostFunction4 Func;
//   *   cost = Func(const std::vector<double> x);
//   *
//   * NOTE that the returned cost is the sqare root of the evaluation
//   *      due to the fact the Nelder Mead class will square the cost.
//   *
//   * Parameters:
//   *   \param x (input) vector containing coordinates of x and y
//   *   \return cost
//   */
//  std::vector<double> operator()(const std::vector<double> &x)
//  {
//    std::vector<double> ret {0};
//    ret[0] = sqrt(fabs(x[0]*x[0] + x[1] - 11)*(x[0]*x[0] + x[1] - 11)
//      + (x[1]*x[1] - 7 + x[0])*(x[1]*x[1] - 7 + x[0]));
//    return ret;
//  }
// };

/**
 * The Sample function 5 is defined as
 *
 *   (x^2 + y^2 + z^2).
 *
 * Number of dimensions = 3
 * The global minimum is: (0, 0, 0)
 * Initial guess:
 *
 * Reference:
 */
class SampleCostFunction5 : public GenericCostFunction
{
 public:
  /**
   * We overload the operator as is required in GenericCostFunction to calculate
   * the cost of the function.
   *
   * Behaviour:
   *   cost = x^2 + y^2 + z^2
   *
   * Intended use :
   *   SampleCostFunction5 Func;
   *   cost = Func(const std::vector<double> x);
   *
   * NOTE that the returned cost is the sqare root of the evaluation
   *      due to the fact the Nelder Mead class will square the cost.
   *
   * Parameters:
   *   \param x (input) vector containing coordinates of x, y and z
   *   \return cost
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    std::vector<double> ret {0};
    ret[0] = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    return ret;
  }
};

/**
 * The Sample function 6 is defined as
 *
 *   (v^2 + w^24 + x^2 + y^2 + z^2).
 *
 * The global minimum is: (0, 0)
 * Initial guess: (0, 0.8)
 *
 * Number of dimensions = 5
 *
 *
 * Reference: http://172.18.53.85/trac/attachment/wiki/
 *     SummerOfCode2012/NelderMead2D.pdf
 */
class SampleCostFunction6 : public GenericCostFunction
{
 public:
  /**
   * We overload the operator as is required in GenericCostFunction to calculate
   * the cost of the function.
   *
   * Behaviour:
   *   cost = v^2 + w^24 + x^2 + y^2 + z^2
   *
   * Intended use :
   *   SampleCostFunction6 Func;
   *   cost = Func(const std::vector<double> x);
   *
   * NOTE that the returned cost is the sqare root of the evaluation
   *      due to the fact the Nelder Mead class will square the cost.
   *
   * Parameters:
   *   \param x (input) vector containing coordinates of v, w, x, y and z
   *   \return cost
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    std::vector<double> ret {0};
    ret[0] = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[4]*x[4]);
    return ret;
  }
};

/**
 * The Sample function 7 is defined as
 *
 *   1/((x-1.05)yz)
 *
 * Number of dimensions = 3
 * Global minimum does not exist
 * Initial guess is thus not required
 *
 * Reference: nil
 */
class SampleCostFunction7 : public GenericCostFunction
{
 public:
  /**
   * We overload the operator as is required in GenericCostFunction to calculate
   * the cost of the Sample function 7 as defined above.
   *
   * Behaviour:
   *   cost = 1/(x-1.05)(yz)
   *
   * Intended use :
   *   SampleCostFunction7 func;
   *   double cost = func(const std::vector<double> x);
   *
   * NOTE that the returned cost is the sqare root of the evaluation
   *      due to the fact the Nelder Mead class will square the cost.
   *
   * Parameters:
   *   \param x (input) vector containing values of the four variables
   *   \return cost
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    std::vector<double> ret {0};
    ret[0] = sqrt(fabs(1/((x[0] - 1.05)*x[1]*x[2])));
    return ret;
  }
};

/**
 * The Sample function 8 is defined as
 *
 *   1/x^2  + 1/y^2
 *
 * Number of dimensions = 2
 * Global minimum does not exist
 * Initial guess is thus not required
 *
 * Reference: nil
 */
class SampleCostFunction8 : public GenericCostFunction
{
 public:
  /**
   * We overload the operator as is required in GenericCostFunction to calculate
   * the cost of the function.
   *
   * Behaviour:
   *   cost = 1/x^2  + 1/y^2
   *
   * Intended use :
   *   SampleCostFunction8 Func;
   *   cost = Func(const std::vector<double> x);
   *
   * NOTE that the returned cost is the sqare root of the evaluation
   *      due to the fact the Nelder Mead class will square the cost.
   *
   * Parameters:
   *   \param x (input) vector containing coordinates of x and y
   *   \return cost
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    std::vector<double> ret {0};
    ret[0] = sqrt(1/(x[0]*x[0]) + 1/(x[1]*x[1]));
    return ret;
  }
};

/**
 * The Sample function 9 is defined as
 *
 *   1/((x-0.95)yz)
 *
 * Number of dimensions = 3
 * Global minimum does not exist
 * Initial guess is thus not required
 *
 * Reference: nil
 */
class SampleCostFunction9 : public GenericCostFunction
{
 public:
  /**
   * We overload the operator as is required in GenericCostFunction to calculate
   * the cost of the Sample function 9 as defined above.
   *
   * Behaviour:
   *   cost = 1/(x-0.95)(yz)
   *
   * Intended use :
   *   SampleCostFunction9 func;
   *   double cost = func(const std::vector<double> x);
   *
   * NOTE that the returned cost is the sqare root of the evaluation
   *      due to the fact the Nelder Mead class will square the cost.
   *
   * Parameters:
   *   \param x (input) vector containing values of the four variables
   *   \return cost
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    std::vector<double> ret {0};
    ret[0] = sqrt(fabs(1/((x[0] - 0.95)*x[1]*x[2])));
    return ret;
  }
};

/**
 * The 2 Dimensional Asymptote function is defined as
 *
 *  -1/(xy+x(x-3))
 *
 * The function has an asymptote at any (x, y) combination whereby
 * (xy + x(x - 3)) = 0.
 *
 * Number of dimensions = 2
 * Example of an asymptote at (0, 0)
 */
class Asymptote2D : public GenericCostFunction
{
 public:
  /**
   * We overload the operator as is required in GenericCostFunction to calculate
   * the cost of the 2 Dimensional Asymptote function as defined above.
   *
   * Behaviour:
   *   cost = -1/(xy+x(x-3))
   *
   * Intended use :
   *   Asymptote2D Func;
   *   cost = Func(const std::vector<double> x);
   *
   * Parameters:
   *   \param x (input) vector containing values of the 3 variables
   *   \return cost
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    std::vector<double> ret {0};
    ret[0] = (fabs(-1/(x[0]*x[1] + x[0]*(x[0] - 3))));
    return ret;
  }
};

/**
 * This Piecewise 2D function is defined as
 *
 *  INF for x<1 and y>1.025
 *  else
 *  5x^2+y^2
 *
 * Number of dimensions = 2
 *
 * This function gives a inf value at (x,y) for x<1 and y>1.025
 */
class FailFindMin : public GenericCostFunction
{
 public:
  /**
   * We overload the operator as is required in GenericCostFunction to calculate
   * the cost of the 2D piecewise function as defined above.
   *
   * Behaviour:
   *   cost = INF for x<1 and y>1.025
   *       5x^2+y^2 for all values
   *
   * Intended use :
   *   FailFindMin Func;
   *   cost = Func(const std::vector<double> x);
   *
   * NOTE that the returned cost is the sqare root of the evaluation
   *      due to the fact the Nelder Mead class will square the cost.
   *
   * Parameters:
   *   \param x (input) vector containing values of the 3 variables
   *   \return cost
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    std::vector<double> ret {0};
    if (x[0] < 1 && x[1]>1.025) {
      ret[0] = 1/(x[0]-x[0]);
    }
    else {
      ret[0] = 5*x[0]*x[0] + x[1]*x[1];
    }
    ret[0] = sqrt(ret[0]);
    return ret;
  }
};

/**
 * This Piecewise 2D function is defined as
 *
 *  (x+5)^2+(y+5)^2+1/(x-102) for all values
 *  and
 *  1/(1-(0.5*x-y-151.5)) for x <105 and y < -100
 *
 * Number of dimensions = 2
 *
 * This function gives a inf value at (x,y) for x<1 and y>1.025
 */
class FailFindMin2 : public GenericCostFunction
{
 public:
  /**
   * We overload the operator as is required in GenericCostFunction to calculate
   * the cost of the 2D piecewise function as defined above.
   *
   * Behaviour:
   *   cost = (x+5)^2+(y+5)^2+1/(x-102) for all values
   *       1/(1-(0.5*x-y-151.5)) for x <105 and y < -100
   *
   * Intended use :
   *   FailFindMin2 Func;
   *   cost = Func(const std::vector<double> x);
   *
   * NOTE that the returned cost is the sqare root of the evaluation
   *      due to the fact the Nelder Mead class will square the cost.
   *
   * Parameters:
   *   \param x (input) vector containing values of the 3 variables
   *   \return cost
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    std::vector<double> ret {0};
    ret[0] = (x[0]-90)*(x[0]-90) + (x[1]+50)*(x[1]+50) + 1/(x[0]-102);
    if (x[0] < 105) {
      ret[0] += 1/(1 - (0.5*x[0]+x[1]+48.5));
    }

    ret[0] = sqrt(fabs(ret[0]));
    return ret;
  }
};

/**
 * The PowellSingularFunction is defined as
 *
 *   (x[0]+10x[1])^2 + 5(x[2]-x[3])^2 + (x[1]-2x[2])^4 + 10(x[0]-x[3])^4
 *
 * Number of dimensions = 4
 * The global minimum is: (0, 0, 0, 0)T
 * Initial guess: (3, -1, 0, 1)T
 *
 * Reference: folk.uib.no/ssu029/Pdf_file/Fletcher77.pdf
 */
class PowellSingular : public GenericCostFunction
{
 public:
  /**
   * We overload the operator as is required in GenericCostFunction to calculate
   * the cost of the function.
   *
   * Behaviour:
   *   cost = (x[0]+10x[1])^2 + 5(x[2]-x[3])^2 +
   *       (x[1]-2x[2])^4 + 10(x[0]-x[3])^4
   *
   * Intended use :
   *   PowellSingular Func;
   *   cost = Func(const std::vector<double> x);
   *
   * Parameters:
   *   \param x (input) vector containing coordinates of x and y
   *   \return cost
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    std::vector<double> ret {0};
    ret[0] = ((x[0] + 10*x[1])*(x[0] + 10*x[1]) + 5*(x[2] -
     x[3])*(x[2] - x[3]) + (x[1] - 2*x[2])*(x[1] - 2*x[2])*(x[1] -
     2*x[2])*(x[1] - 2*x[2]) + 10*(x[0] - x[3])*(x[0] -
     x[3])*(x[0] - x[3])*(x[0] - x[3]));

    ret[0] = sqrt(fabs(ret[0]));
    return ret;
  }
};

/**
 * This function gives an infinite cost for all values of x.
 */
class AlwaysInfinite : public GenericCostFunction
{
 public:
  /**
   * We overload the operator as is required in GenericCostFunction to calculate
   * the cost of the AlwaysInfinite function as defined above.
   *
   * Behaviour:
   *   cost = y = 1000*double_max
   *
   * Intended use :
   *   AlwaysInfinite Func;
   *   cost = Func(const std::vector<double> x);
   *
   * Parameters:
   *   \param x (input) vector containing values of the variable
   *   \return cost
   */
  std::vector<double> operator()(const std::vector<double> &x)
  {
    const double temp = (x[0] >= 1.0) ? x[0] : 1.0;
    std::vector<double> ret {0};
    ret[0] = (std::numeric_limits<double>::max()*1000.0*temp);
    return ret;
  }
};
}  // UnitTests namespace
}  // Unfit namespace
#endif

