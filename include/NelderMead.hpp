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
#ifndef UNFIT_INCLUDE_NELDERMEAD_HPP_
#define UNFIT_INCLUDE_NELDERMEAD_HPP_

#include <vector>
#include "GenericCostFunction.hpp"
#include "GenericOptimizer.hpp"

namespace Unfit
{
/**
 * \brief A class to implement the NelderMead optimization method
 *
 * This class implements the NelderMead optimization method. It requires a
 * cost function (written by the user) and an initial guess of the unknown
 * parameters. It then proceeds to try to find a set of parameters that minimize
 * the cost function.
 *
 * The inspiration for this implementation was derived from the paper:
 *   Implementing the Nelder-Mead simplex algorithm with adaptive parameters,
 *   Journal of Computational Optimization and Applications,
 *   Fuchang Gao & Lixing Han, 51(1):259-277, 2012
 */
class NelderMead : public GenericOptimizer
{
  friend class TestNelderMead;
 public:
  /**
   * The constructor sets all of the parameters to their default values.
   */
  NelderMead();

  /**
   * As we are deriving from this class, the destructor should be virtual. In
   * this case an empty destructor is just fine as we are putting everything
   * into std::vectors which take care of themselves.
   */
  virtual ~NelderMead() {}

  /**
   * \brief Implements the Nelder-Mead optimization method.
   *
   * This method takes in a cost function (user supplied) and an initial guess
   * for the unknown parameters and uses the Nelder-Mead simplex algorithm to
   * try to find the set of parameters that provides a minimum cost. This is
   * implemented as a nonlinear least squares problem, so the cost will always
   * be >= 0. This is known as a local optimizer, so the results you get will
   * depend on how good the initial guess is.
   *
   * Intended use:
   *   auto rc = object.FindMin(CostFunction, coordinates)
   *
   * \param CostFunction A user supplied cost function that returns a vector of
   *        residuals.
   * \param coordinates A vector that contains the initial guess (coordinates)
   *        of the unknown parameters.
   * \return rc return code; zero upon success, non-zero upon failure. \n
   *       0 = Success \n
   *       1 = Max function evaluations reached \n
   *       2 = Max iteration reached \n
   *       3 = Optimisation was unable to find a valid simplex \n
   *       4 = Convergence was achieved but the tolerance(s) were not met \n
   *       GENERATED SIMPLEX: \n
   *      -1 = Invalid initial guess (check cost and bounds) \n
   *      -2 = Unable to generate initial simplex within the bounds \n
   *      -3 = Initial simplex had has invalid costs \n
   *      -4 = Number of coordinates is zero \n
   *       USER SUPPLIED SIMPLEX: \n
   *      -1 = Empty population \n
   *      -2 = Empty population member(s) \n
   *      -3 = Population is the wrong size (vertices != dimensions + 1) \n
   *      -4 = Vertex lengths are inconsistent \n
   *      -5 = One or more vertices has invalid costs
   */
  int FindMin(GenericCostFunction &CostFunction,
      std::vector<double> &coordinates) override;

  /**
   * Resets the variables to default values of NelderMead. Equivalent to
   * destroying the object and creating a new one.
   */
  void Reset() override;

 private:
  /**
   * Operation contains the various operations that are possible at each
   * iteration of the NelderMead algorithm.
   */
  enum Operation {
    Reflected,
    Expanded,
    ContractedIn,
    ContractedOut,
    Shrunk,
    Restarted
  };

  /**
   * \brief Generate the initial simplex by accepting a vector which contains
   * the initial guess and the CostFunction to compute costs.
   *
   * The initial simplex is generated following L.Pfeffer at Stanford. The
   * initial simplex is scaled with respect to the characteristic lengths of
   * the problem. The scaling is defined with respect to the initial point.
   * The method proceeds by defining the scale.
   *   Scale :
   *     non-zero scale : 0.05
   *     zero scale : 0.00025
   *
   * The first vertex is always the initial guess.
   *
   * For the subsequent vertices, we examine ith coordinate and
   * if nth coordinate is zero, scale it up by zero scale while keeping the
   * other coordinates constant. The new set of coordinates will be the
   * nth vertex.
   *
   * If nth coordinate is non-zero, scale the value up by non-zero scale while
   * keeping the other coordinates constant. The new set of coordinates will
   * be the nth vertex.
   *
   * If the coordinates cross a bound, the algorithm tries to find a valid
   * simplex if it can. If the cost of any of the vertices in a geometrically
   * valid simplex are NAN or INF then the algorithm is considered to have
   * failed.
   *
   * NOTES:
   * Vertices of the initial simplex are not sorted under this function.
   * Generated simplex is stored in the 2D vector vertices_
   *
   * \param CostFunction (input) the function to be implemented
   * \param initial_point (input) the initial guess point for the FindMin
   *        function
   * \return  0 = Success. \n
   *          1 = Failure. Invalid initial point \n
   *          2 = Failure. Bounds are too close together to form a simplex \n
   *          3 = Failure. Invalid vertex costs in the initial simplex
   */
  int GeneratePopulation(GenericCostFunction &CostFunction,
      const std::vector<double> &initial_point);

  /**
   * Compute the centroid of n good vertices without computing the cost
   * function by averaging the positions of all of the vertices except the
   * vertex with the worst cost.
   */
  void ComputeCentroid();

  /**
   * ContractInside is called by ProcessFindMin. As the name suggests, it
   * performs an inside contraction, which creates a new candidate point for
   * the simplex between the worst point and the centroid. The distance is
   * controlled by the gamma parameter.
   *
   *   contractin = centroid - gamma*(reflect - centroid)
   *
   * \param CostFunction (input) the function to be implemented
   * \return true = Success \n
   *         false = Failure. Cannot contract and maintain a valid simplex
   */
  bool ContractInside(GenericCostFunction &CostFunction);

  /**
   * ContractOutside is called by ProcessFindMin. As the name suggests, it
   * performs an outside contraction, which creates a new candidate point for
   * the simplex between the reflect point and the centroid. The distance is
   * controlled by the gamma parameter.
   *
   *   contractout = centroid + gamma*(reflect - centroid_[i])
   *
   * \param CostFunction (input) the function to be implemented
   * \return true = Success \n
   *         false = Failure. Cannot contract and maintain a valid simplex
   */
  bool ContractOutside(GenericCostFunction &CostFunction);

  /**
   * Compute the centroid and the reflection point and evaluate the value of R
   * based on alpha (alpha>0). This computes the centroid of best and good
   * vertices, and get the reflected point with respect to the worst vertex and
   * alpha.
   *
   *   reflect = centroid_ + alpha*(centroid_ - worst)
   *
   * \param CostFunction (input) the function to be implemented
   * \return true = Success \n
   *         false = Failure. Cannot reflect and maintain a valid simplex
   */
  bool Reflect(GenericCostFunction &CostFunction);

  /**
   * Computes the expanded coordinates and their cost function based on beta
   *
   *   expand = centroid + beta*(reflect - centroid)
   *
   * \param CostFunction (input) the function specified by the user
   * \return true = Success \n
   *         false = Failure. Cannot expand and maintain a valid simplex
   */
  bool Expand(GenericCostFunction &CostFunction);

  /**
   * Shrink shifts Good and Worst vertices towards Best vertex
   * It shifts all the vertices except the Best vertex, based on delta
   *
   * \param CostFunction (input) the function to be implemented
   * \return true = Success \n
   *         false = Failure. Cannot shrink and maintain a valid simplex
   */
  bool Shrink(GenericCostFunction &CostFunction);

  /**
   * Check if the vertices in the simplex have become colinear. If they have
   * the simplex is considered degenerate.
   *
   * \return true Simplex is degenerate \n
   *         false Simplex is not degenerate
   */
  bool IsDegenerate();

  /**
   * Initialise the working vectors vertices_, contract_, centroid_,
   * reflect_, expand_ to size depending on the dimension of the problem.
   */
  void InitialiseVectors();

  /**
   * Set the parameters to be used based on the number of dimensions
   */
  void UseAdaptiveParameters();

  /**
   * Prints the column headers relating to iteration by iteration output to the
   * screen if the output level has been set to be >= 1.
   *
   * \param best_cost The best cost after the initial setup
   */
  void PrintInitialOutput(double best_cost) const override;

  /**
   * At each iteration, prints nothing if the output level is zero, the
   * iteration number if the output level is one, and the iteration number plus
   * the number of function evaluations and the best cost if the output level is
   * greater than one, along with the operation from the last iteration.
   *
   * \param best_cost The best cost after the latest iteration
   */
  void PrintIterationOutput(double best_cost) const override;

  /**
   * Here is the main implementation of the Nelder Mead simplex algorithm to
   * which the FindMin methods are a front-end.
   *
   * \param CostFunction (input) a user function
   * \return rc return code; zero upon success, non-zero upon failure. \n
   *       0 = Success \n
   *       1 = Max function evaluations reached \n
   *       2 = Max iteration reached \n
   *       3 = Optimisation was unable to find a valid simplex \n
   *       4 = Convergence was achieved but the tolerance(s) were not met
   */
  int ProcessFindMin(GenericCostFunction &CostFunction);

  /**
   * Creates a new simplex around the best vertex by making a call to
   * GeneratePopulation internally.
   *
   * \param CostFunction (input) the function to be implemented
   * \return 0 = Success \n
   *         1 = Unable to generate valid simplex
   */
  int RegeneratePopulation(GenericCostFunction &CostFunction);

  /** Vector to store the contracted point*/
  std::vector<double> contract_;
  /** Vector to store the centroid*/
  std::vector<double> centroid_;
  /** Vector to store the reflected point*/
  std::vector<double> reflect_;
  /** Vector to store the expanded point*/
  std::vector<double> expand_;
  /** Variable to store the number of dimensions*/
  std::size_t dimensions_;
  /** The index of the best vertex*/
  std::size_t best_vertex_;
  /** The index of the worst vertex*/
  std::size_t worst_vertex_;
  /** The index of the next-to-worst vertex*/
  std::size_t next_worst_vertex_;
  /** The index of the cost of each vertex*/
  std::size_t cost_;
  /** Variable to store the best cost before a restart*/
  double restart_best_cost_;
  /** Enum data type to store the operation executed*/
  enum Operation process_;
};

}  // namespace Unfit

#endif
