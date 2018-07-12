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
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>
#include "GenericCostFunction.hpp"
#include "NelderMead.hpp"
#include "NelderMeadTestFunctions.hpp"
#include "UnitTest++.h"

static const double epsilon(1e-12);

namespace Unfit
{
/**
 * \brief Access to the private member functions in the neldermead class
 *
 * The main goal of this class is to access the private member functions in the
 * neldermead class in tests. See the documentation of each function for more
 * details.
 */
class TestNelderMead : public NelderMead
{
 public:
  TestNelderMead()=default;

  /**
   * create a function to access the private method(GeneratePopulation) in
   * class NelderMead.
   *
   * Intended use:
   *   AccessGeneratePopulation(CostFunction, min_point)
   *
   * Parameters:
   *   \param CostFunction (input) the function to be implemented
   *   \param initial_point (input) vector which contain the initial guess
   *
   * Return Codes
   *   \return GeneratePopulation(CostFunction, initial_point)
   */
  int AccessGeneratePopulation(GenericCostFunction &CostFunction,
      std::vector<double> &initial_point)
  {
    return GeneratePopulation(CostFunction, initial_point);
  }

/**
 * create a function to access the private method(ComputeCentroid) in
 * class NelderMead.
 *
 * Intended use:
 *   rc = AccessComputeCentroid()
 */
  void AccessComputeCentroid()
  {
    ComputeCentroid();
  }

  /**
   * create a function to access the private method (ContractInside) in
   * class NelderMead.
   *
   * Intended use:
   *   bool rc = AccessContractInside(CostFunction);
   *
   * Parameters:
   *   \param CostFunction (input) the function to be implemented
   *
   * Return Codes:
   *  \return true = Success; false = Cannot compute
   */
  bool AccessContractInside(GenericCostFunction &CostFunction)
  {
    return ContractInside(CostFunction);
  }

  /**
   * create a function to access the private method (ContractOutside) in
   * class NelderMead.
   *
   * Intended use:
   *   bool rc = AccessContractOutside(CostFunction);
   *
   * Parameters:
   *   \param CostFunction (input) the function to be implemented
   *
   * Return Codes:
   *  \return true = Success; false = Cannot compute
   */
  bool AccessContractOutside(GenericCostFunction &CostFunction)
  {
    return ContractOutside(CostFunction);
  }

  /**
   * create a function to access the private method(Reflect) in
   * class NelderMead.
   *
   * Intended use:
   *   bool rc = AccessReflect(CostFunction);
   *
   * Parameters:
   *   \param CostFunction (input) the function to be implemented
   *
   * Return Codes:
   *  \return true = Success; false = Cannot compute
   */
  bool AccessReflect(GenericCostFunction &CostFunction)
  {
    return Reflect(CostFunction);
  }

  /**
   * create a function to access the private method(Expand) in
   * class NelderMead.
   *
   * Intended use:
   *   bool rc = AccessExpand(CostFunction);
   *
   * Parameters:
   *   \param CostFunction (input) the function specified by the user
   *
   * Return Codes:
   *   \return true = Success; false = Cannot compute
   */
  bool AccessExpand(GenericCostFunction &CostFunction)
  {
    return Expand(CostFunction);
  }

/**
 * create a function to access the private method(Shrink) in class NelderMead.
 * Intended use:
 *   bool rc = AccessShrink(CostFunction);
 *
 * Parameters:
 *   \param CostFunction (input) the function to be implemented
 *
 * Return Codes:
 *   \return true = Success; false = When the points are the same
 */
  bool AccessShrink(GenericCostFunction &CostFunction)
  {
    return Shrink(CostFunction);
  }

  /**
   * Initialise the working vectors population_, contract_, centroid_,
   * reflect_, expand_ to size depending on the dimension of the problem.
   * Elements are assigned 0 value.
   *
   * Behaviour:
   *  For population_:
   *      population_[number_of_dimensions+1][number_of_dimensions+1]
   *  For the rest:
   *      workingvector_[number_of_dimensions+1]
   *
   * Intended use:
   *   AccessInitialiseVectors();
   */
  void AccessInitialiseVectors()
  {
    InitialiseVectors();
  }

  /**
   * Calls the PrintIterationOutput method of NelderMead. Requires that the
   * number of iterations, number of function evaluations, cost of the best
   * vertex, and the process be set.
   *
   * \param best_cost The best cost at the current iteration
   */
  void AccessPrintIterationOutput(double best_cost)
  {
    PrintIterationOutput(best_cost);
  }

  /**
   * create a function to access the private method(RegeneratePopulation) in
   * class NelderMead.
   *
   * Intended use:
   *   int rc = AccessRegeneratePopulation(CostFunction);
   *
   * Parameters:
   *   \param CostFunction (input) the function specified by the user
   *   \return rc return code; zero upon success, non-zero upon failure
   *
   * Return Codes:
   *   0 = Success
   *   1 = invalid simplex
   */
  int AccessRegeneratePopulation(GenericCostFunction &CostFunction)
  {
    return RegeneratePopulation(CostFunction);
  }

  /**
   * Set the values stored in the vector centroid_.
   * NOTE: the size of the vector that is passed in must be 3.
   *
   * Intended use:
   *   SetCentroid(centroid_values);
   *
   * Parameters:
   *   \param centroid_values (input) a vector of values to be put inside the
   *          centroid_ vector
   */
  void SetCentroid(std::vector<double> &centroid_values)
  {
     centroid_ = centroid_values;
  }

  /**
   * Set the values stored in the vector reflect_.
   * NOTE: the size of the vector that is passed in must be 3.
   *
   * Intended use:
   *   SetReflect(Reflect_values);
   *
   * Parameters:
   *  \param reflect_values (input) a vector of values to be put inside the
   *         reflect_ vector
   */
  void SetReflect(std::vector<double> &reflect_values)
  {
    reflect_.assign(reflect_values.size(), 0);
    reflect_ = reflect_values;
  }

 /**
   * Set the values stored in the vector contract_.
   * NOTE: the size of the vector that is passed in must be 3.
   *
   * Intended use:
   *   SetContract(Contract_values);
   *
   * Parameters:
   *  \param contract_values (input) a vector of values to be put inside the
   *         contract_ vector
   */
  void SetContract(std::vector<double> &contract_values)
  {
    contract_.assign(contract_values.size(), 0);
    contract_ = contract_values;
  }

  /**
   * Set the number of dimensions stored in the variable dimensions_.
   * NOTE: the size of the variable must be 2 and above.
   *
   * Intended use:
   *   SetNumberOfDimensions_(number of dimensions);
   *
   * Parameters:
   *  \param number_of_dimensions (input) an integer to be the number of
   *         dimensions
   */
  void SetNumberOfDimensions_(std::size_t number_of_dimensions)
  {
    dimensions_ = number_of_dimensions;
    cost_ = dimensions_;
  }

  /**
   * Get the number of dimensions stored in the variable dimensions_.
   * NOTE: the size of the variable must be 2 and above.
   *
   * Intended use:
   *   unsigned rc = GetNumberOfDimensions_();
   *
   * Return Codes:
   *  \return rc = number of dimensions
   */
  std::size_t GetNumberOfDimensions_()
  {
    return dimensions_;
  }

  /**
   * Gets the expanded coordinates and f(x, y) after the expand function. The
   * vector is stored in a temporary vector store_expand. It returns 0 when
   * it is successful.
   *
   * Intended use:
   *   unsigned rc = GetExpand(store_expand);
   *
   * Parameters:
   *   \param store_expand (input) an empty vector to store the expand_ vector
   *          that contains the expand point
   *
   * Return Codes:
   *  \return 0 = Success; 1 = Cannot get the vector
   */
  unsigned GetExpand(std::vector<double> &store_expand) const
  {
    store_expand = expand_;
    return 0;
  }

  /**
   * Gets the centroid coordinates and f(x, y) of the mid point. The mid point
   * vector is stored in a temporary vector store_centroid_. Returns 0 when
   * it is successful.
   *
   * Intended use:
   *   unsigned rc = GetCentroid(store_mid_point);
   *
   * Parameters:
   *   \param store_centroid (input) a vector of values to be put inside the
   *          centroid_ vector
   *
   * Return Codes:
   *  \return 0 = Success; 1 = Cannot compute
   */
  unsigned GetCentroid(std::vector<double> &store_centroid) const
  {
    store_centroid = centroid_;
    return 0;
  }

  /**
   * Gets the reflect point coordinates and f(x, y) of the reflect point. The
   * reflect vector is stored in a temporary vector store_reflect_
   *
   * Intended use:
   *   unsigned rc = GetCentroid(store_reflect);
   *
   * Parameters:
   *   \param store_reflect (input) a vector of values to be put inside the
   *          reflect_ vector
   *
   * Return Codes:
   *  \return 0 = Success; 1 = Cannot get the vector
   */
  unsigned GetReflect(std::vector<double> &store_reflect) const
  {
    store_reflect = reflect_;
    return 0;
  }

  /**
   * Gets the expanded coordinates and f(x, y) after the contract function. The
   * vector is stored in a temporary vector store_contract. It returns 0 when
   * it is successful.
   *
   * Intended use:
   *   unsigned rc = GetContract(store_contract);
   *
   * Parameters:
   *   \param store_contract (input) vector to be filled in with contract_
   *          values obtained from the Contract function
   *
   * Return Codes:
   *  \return 0 = Success; 1 = Cannot get the vector
   */
  unsigned GetContract(std::vector<double> &store_contract) const
  {
    store_contract = contract_;
    return 0;
  }

  /**
   * Check if the termination criteria are met and stops the NelderMead
   * iteration.
   * Two termination criteria are checked against a user-defined tolerance, and
   * both have to be met.
   * Criteria 1:
   * Check how different the minima calculated are - difference in costs
   * comparing the best vertex to the other vertices
   * Criteria 2:
   * Check how different the coordinates are relative to one another, comparing
   * the best vertex to the other vertices
   *
   * Intended use:
   *   AccessIsDegenerate()
   *
   * Return Codes:
   *  \return false = Termination criteria not met, NOT okay to stop iterating;
   *          true = Termination criteria met, okay to stop iterating
   */
  bool AccessIsDegenerate()
  {
    return IsDegenerate();
  }

  /**
   * Set the values stored in the simplex as a vector of vectors
   * NOTE: the size of the vector that is passed in must be at least 2 and
   * include the cost of the vertex. Size of the vector vector should be a
   * matrix of (n+1) x (n+1) where n is the number of dimensions.
   * This function also helps to set the number of dimension based on the size
   * of the vectors given.
   *
   * Intended use:
   *   SetVertices(vertex);
   *
   * Parameters:
   *  \param vertex (input) a vector of vectors to be put inside the
   *         population_ vector
   */
  void SetVertices(std::vector< std::vector<double> > &vertex)
  {
    if (vertex.size() >= 2 && vertex[0].size() <= vertex.size()) {
      population_.assign(vertex.size(),
          std::vector<double>(vertex[0].size(), 0));
      population_ = vertex;
      SetNumberOfDimensions_(vertex[0].size()-1);
    }
  }

  /**
   * Get the value of the parameter current_function_evaluations
   *
   * Intended use:
   *   unsigned rc = GetCurrentFunctionEvaluation();
   *
   * Return Code:
   *   \return number of function evaluations
   */
  std::size_t GetCurrentFunctionEvaluation()
  {
    return (function_evaluations_);
  }

  /**
   * Get the value of the parameter current_iterations_
   *
   * Intended use:
   *   unsigned rc = GetCurrentIterations();
   *
   * Return Code:
   *   \return number of iterations;
   */
  std::size_t  GetCurrentIterations()
  {
    return (iterations_);
  }

  /**
   * Set the value of the Operation enum (process_ variable). See the
   * documentation for NelderMead.hpp to see what entries are valid.
   *
   * Parameters:
   *   \param process Must be selected from the Operation enum
   */
  void SetProcess(int process)
  {
    // Use static case just to make it clear we are using an int to set the
    // value of the enum. This is because Operation is private and therefore
    // we can't see the names associated with Operation/process_.
    process_ = static_cast<Operation>(process);
  }

  /**
   * Set the vertex that has the worst cost. If it is not clear which one it
   * should be, just call AccessSortPopulation, and then the index of the vertex
   * with the worst cost will be the same as the number of dimensions.
   *
   * Parameters:
   *   \param vertex_number The index of the vertex with the worst cost
   *     (index starts from zero and ends at dimensions_)
   */
  void SetWorstVertex(unsigned vertex_number)
  {
    if (vertex_number <= dimensions_) worst_vertex_ = vertex_number;
  }
};

namespace UnitTests
{
SUITE(UnitTestNelderMead)
{
// Test whether the constructor correctly inputs the given vector to population_
// by checking each elements in the population_
TEST_FIXTURE(TestNelderMead, SettingVertices)
{
  // Initialise the vector of vector to set the private member population_
  std::vector<std::vector<double>> v1 = {{2, 4, 4},
      {56, 234234, -4}, {45, 12412, 2}};
  SetVertices(v1);
  // Create three empty vector to be used to get the vertex out from the private
  // member
  std::vector<double> best = GetSolution(0);
  best.push_back(GetCost(0));
  std::vector<double> good = GetSolution(1);
  good.push_back(GetCost(1));
  std::vector<double> worst = GetSolution(2);
  worst.push_back(GetCost(2));
  // Check the vector that is extracted out with what was being set.
  CHECK_CLOSE(2, best[0], epsilon);
  CHECK_CLOSE(4, best[1], epsilon);
  CHECK_CLOSE(4, best[2], epsilon);
  CHECK_CLOSE(56, good[0], epsilon);
  CHECK_CLOSE(234234, good[1], epsilon);
  CHECK_CLOSE(-4, good[2], epsilon);
  CHECK_CLOSE(45, worst[0], epsilon);
  CHECK_CLOSE(12412, worst[1], epsilon);
  CHECK_CLOSE(2, worst[2], epsilon);
}

// Test whether SortPopulation correctly sorts vertices in the correct order
// according to the costs by checking each element in the population_ after sort
TEST_FIXTURE(TestNelderMead, SortPopulationNegativeIntegers)
{
  // Initialise vertices. Set them and sort them
  std::vector<std::vector<double>> v1 = {{2, 4, 4},
      {56, 234234, -4}, {45, 12412, 2}};
  SetVertices(v1);
  SortPopulation();

  // Empty vectors to get sorted vertices
  std::vector<double> v4 = GetSolution(0);
  v4.push_back(GetCost(0));
  std::vector<double> v5 = GetSolution(1);
  v5.push_back(GetCost(1));
  std::vector<double> v6 = GetSolution(2);
  v6.push_back(GetCost(2));

  // Check if vertices have been sorted properly
  CHECK_CLOSE(56, v4[0], epsilon);
  CHECK_CLOSE(234234, v4[1], epsilon);
  CHECK_CLOSE(-4, v4[2], epsilon);
  CHECK_CLOSE(45, v5[0], epsilon);
  CHECK_CLOSE(12412, v5[1], epsilon);
  CHECK_CLOSE(2, v5[2], epsilon);
  CHECK_CLOSE(2, v6[0], epsilon);
  CHECK_CLOSE(4, v6[1], epsilon);
  CHECK_CLOSE(4, v6[2], epsilon);
}

// Test whether SortPopulation correctly sorts vertices in the correct order
// this test uses a vector contains vectices in a different order
TEST_FIXTURE(TestNelderMead, SortPopulationPositiveIntegers)
{
  // Initialise vertices. Set them and sort them
  std::vector<std::vector<double>> v1 = {{2, 4, 3},
      {56, 234234, 2}, {45, 12412, 1}};
  SetVertices(v1);
  SortPopulation();

  // Empty vectors to get sorted vertices
  std::vector<double> v4 = GetSolution(0);
  v4.push_back(GetCost(0));
  std::vector<double> v5 = GetSolution(1);
  v5.push_back(GetCost(1));
  std::vector<double> v6 = GetSolution(2);
  v6.push_back(GetCost(2));

  // Check if vertices have been sorted properly
  CHECK_CLOSE(45, v4[0], epsilon);
  CHECK_CLOSE(12412, v4[1], epsilon);
  CHECK_CLOSE(1, v4[2], epsilon);
  CHECK_CLOSE(56, v5[0], epsilon);
  CHECK_CLOSE(234234, v5[1], epsilon);
  CHECK_CLOSE(2, v5[2], epsilon);
  CHECK_CLOSE(2, v6[0], epsilon);
  CHECK_CLOSE(4, v6[1], epsilon);
  CHECK_CLOSE(3, v6[2], epsilon);
}

// this test involves cost function with negative values
TEST_FIXTURE(TestNelderMead, SortPopulationNegativeDecimals)
{
  // Initialise vertices. Set them and sort them
  std::vector<std::vector<double>> v1 = {{2, 4, 2.1},
      {2, 5, -1.2}, {-2, -3, 3.5}};
  SetVertices(v1);
  SortPopulation();

  // Empty vectors to get sorted vertices
  std::vector<double> v4 = GetSolution(0);
  v4.push_back(GetCost(0));
  std::vector<double> v5 = GetSolution(1);
  v5.push_back(GetCost(1));
  std::vector<double> v6 = GetSolution(2);
  v6.push_back(GetCost(2));

  // Check if vertices have been sorted properly
  CHECK_CLOSE(2, v4[0], epsilon);
  CHECK_CLOSE(5, v4[1], epsilon);
  CHECK_CLOSE(-1.2, v4[2], epsilon);
  CHECK_CLOSE(2, v5[0], epsilon);
  CHECK_CLOSE(4, v5[1], epsilon);
  CHECK_CLOSE(2.1, v5[2], epsilon);
  CHECK_CLOSE(-2, v6[0], epsilon);
  CHECK_CLOSE(-3, v6[1], epsilon);
  CHECK_CLOSE(3.5, v6[2], epsilon);
}

//
TEST_FIXTURE(TestNelderMead, SortPopulationNegativeCoordinates)
{
  // Initialise vertices. Set them and sort them
  std::vector<std::vector<double>> v1 = {{2, 4, -2.5},
      {2, 5, 2}, {-2, -3, -3}};
  SetVertices(v1);
  SortPopulation();

  // Empty vectors to get sorted vertices
  std::vector<double> v4 = GetSolution(0);
  v4.push_back(GetCost(0));
  std::vector<double> v5 = GetSolution(1);
  v5.push_back(GetCost(1));
  std::vector<double> v6 = GetSolution(2);
  v6.push_back(GetCost(2));

  // Check if vertices have been sorted properly
  CHECK_CLOSE(-2, v4[0], epsilon);
  CHECK_CLOSE(-3, v4[1], epsilon);
  CHECK_CLOSE(-3, v4[2], epsilon);
  CHECK_CLOSE(2, v5[0], epsilon);
  CHECK_CLOSE(4, v5[1], epsilon);
  CHECK_CLOSE(-2.5, v5[2], epsilon);
  CHECK_CLOSE(2, v6[0], epsilon);
  CHECK_CLOSE(5, v6[1], epsilon);
  CHECK_CLOSE(2, v6[2], epsilon);
}

TEST_FIXTURE(TestNelderMead, SortPopulationNegativeCoordinatesAndCost)
{
  // Initialise vertices. Set them and sort them
  std::vector<std::vector<double>> v1 = {{2, 4, -100},
      {2, 5, -99}, {-2, -3, -98}};
  SetVertices(v1);
  SortPopulation();

  // Empty vectors to get sorted vertices
  std::vector<double> v4 = GetSolution(0);
  v4.push_back(GetCost(0));
  std::vector<double> v5 = GetSolution(1);
  v5.push_back(GetCost(1));
  std::vector<double> v6 = GetSolution(2);
  v6.push_back(GetCost(2));

  // Check if vertices have been sorted properly
  CHECK_CLOSE(2, v4[0], epsilon);
  CHECK_CLOSE(4, v4[1], epsilon);
  CHECK_CLOSE(-100, v4[2], epsilon);
  CHECK_CLOSE(2, v5[0], epsilon);
  CHECK_CLOSE(5, v5[1], epsilon);
  CHECK_CLOSE(-99, v5[2], epsilon);
  CHECK_CLOSE(-2, v6[0], epsilon);
  CHECK_CLOSE(-3, v6[1], epsilon);
  CHECK_CLOSE(-98, v6[2], epsilon);
}

TEST_FIXTURE(TestNelderMead, SortPopulationDecimalCoordinateAndCost)
{
  // Initialise vertices. Set them and sort them
  std::vector<std::vector<double>> v1 = {{2.1, 4.2, 3.3},
      {2.8, 5.1, 3.33}, {-2.3, -3.3, 3.31}};
  SetVertices(v1);
  SortPopulation();

  // Empty vectors to get sorted vertices
  std::vector<double> v4 = GetSolution(0);
  v4.push_back(GetCost(0));
  std::vector<double> v5 = GetSolution(1);
  v5.push_back(GetCost(1));
  std::vector<double> v6 = GetSolution(2);
  v6.push_back(GetCost(2));

  // Check if vertices have been sorted properly
  CHECK_CLOSE(2.1, v4[0], epsilon);
  CHECK_CLOSE(4.2, v4[1], epsilon);
  CHECK_CLOSE(3.3, v4[2], epsilon);
  CHECK_CLOSE(-2.3, v5[0], epsilon);
  CHECK_CLOSE(-3.3, v5[1], epsilon);
  CHECK_CLOSE(3.31, v5[2], epsilon);
  CHECK_CLOSE(2.8, v6[0], epsilon);
  CHECK_CLOSE(5.1, v6[1], epsilon);
  CHECK_CLOSE(3.33, v6[2], epsilon);
}

// Tests SortPopulation of 10 dimensions.
TEST_FIXTURE(TestNelderMead, SortPopulation10D)
{
  // Initialise vertices. Set them and sort them
  std::vector<std::vector<double>> v1 = {
      {2.1, 4.2, 3.3, 2, 34, 1.2, 4.3, 12, 42, 2, 2},
      {2.1, 4.2, 3.3, 2, 34, 1.2, 4.3, 12, 42, 1, 1},
      {2.1, 4.2, 3.3, 2, 34, 1.2, 4.3, 12, 42, 1, 3},
      {2.1, 4.2, 3.3, 2, 34, 1.2, 4.3, 12, 42, 1, 4},
      {2.1, 4.2, 3.3, 2, 34, 1.2, 4.3, 12, 42, 1, 5},
      {2.1, 4.2, 3.3, 2, 34, 1.2, 4.3, 12, 42, 1, 6},
      {2.1, 4.2, 3.3, 2, 34, 1.2, 4.3, 12, 42, 1, 7},
      {2.1, 4.2, 3.3, 2, 34, 1.2, 4.3, 12, 42, 1, 8},
      {2.1, 4.2, 3.3, 2, 34, 1.2, 4.3, 12, 42, 1, 10},
      {2.1, 4.2, 3.3, 2, 34, 1.2, 4.3, 12, 42, 1, 9},
      {2.1, 4.2, 3.3, 2, 34, 1.2, 4.3, 12, 42, 1, 11},
      {2.1, 4.2, 3.3, 2, 34, 1.2, 4.3, 12, 42, 1, 55}};
  SetVertices(v1);
  SortPopulation();

  // Empty vectors to get sorted vertices
  std::vector<double> v2 = GetSolution(0);
  v2.push_back(GetCost(0));
  std::vector<double> v3 = GetSolution(1);
  v3.push_back(GetCost(1));
  std::vector<double> v4 = GetSolution(2);
  v4.push_back(GetCost(2));
  std::vector<double> v5 = GetSolution(3);
  v5.push_back(GetCost(3));
  std::vector<double> v6 = GetSolution(4);
  v6.push_back(GetCost(4));
  std::vector<double> v7 = GetSolution(5);
  v7.push_back(GetCost(5));
  std::vector<double> v8 = GetSolution(6);
  v8.push_back(GetCost(6));
  std::vector<double> v9 = GetSolution(7);
  v9.push_back(GetCost(7));
  std::vector<double> v10 = GetSolution(8);
  v10.push_back(GetCost(8));
  std::vector<double> v11 = GetSolution(9);
  v11.push_back(GetCost(9));
  std::vector<double> v12 = GetSolution(10);
  v12.push_back(GetCost(10));

  // Check if vertices have been sorted properly
  CHECK_CLOSE(1, v2[10], epsilon);
  CHECK_CLOSE(2, v3[10], epsilon);
  CHECK_CLOSE(3, v4[10], epsilon);
  CHECK_CLOSE(4, v5[10], epsilon);
  CHECK_CLOSE(5, v6[10], epsilon);
  CHECK_CLOSE(6, v7[10], epsilon);
  CHECK_CLOSE(7, v8[10], epsilon);
  CHECK_CLOSE(8, v9[10], epsilon);
  CHECK_CLOSE(9, v10[10], epsilon);
  CHECK_CLOSE(10, v11[10], epsilon);
  CHECK_CLOSE(11, v12[10], epsilon);
}

// Test whether GetCentroid get the correct value stored in centroid_
TEST_FIXTURE(TestNelderMead, GetCentroid)
{
  // Initialise and set centroid values
  std::vector<double> centroid_values = {0.6, 0.4, 0};
  SetCentroid(centroid_values);

  // Get centroid values and check
  std::vector<double> getcentroid;
  unsigned rc = GetCentroid(getcentroid);
  CHECK_EQUAL(0u, rc);
  CHECK_CLOSE(0.6, getcentroid[0], epsilon);
  CHECK_CLOSE(0.4, getcentroid[1], epsilon);
}

// Test whether GetReflect get the correct value stored in reflect_
TEST_FIXTURE(TestNelderMead, GetReflect)
{
  // Initialise and set reflect values
  std::vector<double> reflect_values = {1.2, 0.8, 0};
  SetReflect(reflect_values);

  // Get reflect values and check
  std::vector<double> getreflect;
  unsigned rc = GetReflect(getreflect);
  CHECK_CLOSE(1.2, getreflect[0], epsilon);
  CHECK_CLOSE(0.8, getreflect[1], epsilon);
  CHECK_EQUAL(0u, rc);
}

// Test whether ComputeCentroid correctly computes the centroid (the midpoint
// of a line in this case) by checking the value stored in centroid_ in
// 2 dimensions
TEST_FIXTURE(TestNelderMead, ComputeCentroid2D)
{
  // Create our object with defined vertices and cost
  std::vector<std::vector<double>> v1 = {{1.2, 0, -3.36},
      {0, 0.8, -0.16}, {0, 0.9, -0.16}};
  SetVertices(v1);
  SetWorstVertex(2);

  // Compute centroid and test centroid coordinates
  AccessComputeCentroid();
  std::vector<double> centroid;
  GetCentroid(centroid);
  // Check x-coordinate of centroid
  CHECK_CLOSE(0.6, centroid[0], epsilon);
  // Check y-coordinate of centroid
  CHECK_CLOSE(0.4, centroid[1], epsilon);
}

// Test whether ComputeCentroid works in 4 dimensions
TEST_FIXTURE(TestNelderMead, ComputeCentroid4D)
{
  // Create our object with defined vertices and cost
  std::vector<std::vector<double>> v1 = {{1, 1, 0, 0, 5},
      {0, 1, 0, 0, 4}, {0, 0, 1, 0, 3}, {0, 0, 0, 1, 2}, {0, 0, 0, 0, 6}};
  SetVertices(v1);
  // Set Number of Dimensions to allow computation of centroid
  SetNumberOfDimensions_(4);
  SetWorstVertex(4);

  // Sort and then ComputeCentroid
  SortPopulation();
  AccessComputeCentroid();

  // Initialise an empty vector and get the best vertice and test
  std::vector<double> v2 = GetSolution(0);
  v2.push_back(GetCost(0));
  CHECK_CLOSE(0, v2[0], epsilon);
  CHECK_CLOSE(0, v2[1], epsilon);
  CHECK_CLOSE(0, v2[2], epsilon);
  CHECK_CLOSE(1, v2[3], epsilon);
  CHECK_CLOSE(2, v2[4], epsilon);

  //  Initialise an empty vector to get the computed centroid
  std::vector<double> centroid;
  GetCentroid(centroid);
  CHECK_CLOSE(0.25, centroid[0], epsilon);
  CHECK_CLOSE(0.5, centroid[1], epsilon);
  CHECK_CLOSE(0.25, centroid[2], epsilon);
  CHECK_CLOSE(0.25, centroid[3], epsilon);
}

// Test whether Contract_out correctly computes outward contracted point by
// checking the value stored in contract_ for 3D function
TEST_FIXTURE(TestNelderMead, Contract_out1)
{
  // Initialise the population_ with a generated simplex
  std::vector<std::vector<double>> v1 = { {0, 0, 0, 0}, {0, 0, 0, 0},
      {1, 0, 0, 1}, {1, 1, 1, 3} };

  // Set number of dimensions and initialise all the working vectors
  SetNumberOfDimensions_(3);
  AccessInitialiseVectors();
  SetVertices(v1);
  SetWorstVertex(3);

  // Compute the Centroid
  AccessComputeCentroid();
  std::vector<double> centroid;
  GetCentroid(centroid);
  CHECK_CLOSE((1.0/3.0), centroid[0], epsilon);
  CHECK_CLOSE(0.0, centroid[1], epsilon);
  CHECK_CLOSE(0.0, centroid[2], epsilon);

  // Compute the contract point using the outward contract function
  SampleCostFunction5 Func;
  CHECK(AccessContractOutside(Func));

  // Initialise an empty vector and get the contracted point and test the values
  std::vector<double> store_contract;
  GetContract(store_contract);
  CHECK_CLOSE(0.0, store_contract[0], 1e-6);
  CHECK_CLOSE(-0.5, store_contract[1], 1e-6);
  CHECK_CLOSE(-0.5, store_contract[2], 1e-6);
}

// Test whether Contract_out correctly computes outward contracted point by
// checking the value stored in contract_ for 5D function
TEST_FIXTURE(TestNelderMead, Contract_out2)
{
  // Initialise the population_ with a generated simplex
  std::vector<std::vector<double>> v1 = { {0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0},
      {1, 1, 1, 0, 0, 3}, {1, 1, 1, 1, 1, 5} };

  // Set number of dimensions and initialise all the working vectors
  SetNumberOfDimensions_(5);
  AccessInitialiseVectors();
  SetVertices(v1);
  SetWorstVertex(5);
  AccessComputeCentroid();

  // Compute the contract point using the outward contract function
  SampleCostFunction6 Func;
  CHECK(AccessContractOutside(Func));

  // Initialise an empty vector and get the contracted point and test the values
  std::vector<double> store_contract;
  GetContract(store_contract);
  CHECK_CLOSE(-0.2, store_contract[0], epsilon);
  CHECK_CLOSE(-0.2, store_contract[1], epsilon);
  CHECK_CLOSE(-0.2, store_contract[2], epsilon);
  CHECK_CLOSE(-0.5, store_contract[3], epsilon);
  CHECK_CLOSE(-0.5, store_contract[4], epsilon);
  CHECK_CLOSE(0.62, store_contract[5], epsilon);
}


TEST_FIXTURE(TestNelderMead, ContractOutsideInvalidDomain)
{
  // Initialise the population_ with a generated simplex
  std::vector<std::vector<double>> v1 = { {0, 0, 0}, {1, 0, 2}, {1, 1, 3} };

  // Set number of dimensions and initialise all the working vectors
  SetNumberOfDimensions_(2);
  AccessInitialiseVectors();
  SetVertices(v1);
  SetWorstVertex(2);
  bounds.SetBounds(1, -0.1, 10);
  AccessComputeCentroid();
  SampleCostFunction5 Func;

  // Compute the contract point using the outward contract function
  CHECK(!AccessContractOutside(Func));
}

// Test whether Contract_in correctly computes inward contracted point by
// checking the value stored in contract_ for 3D function
TEST_FIXTURE(TestNelderMead, Contract_in1)
{
  // Initialise the population_ with a generated simplex
  std::vector<std::vector<double>> v1 = { {0, 0, 0, 0}, {0, 0, 0, 0},
      {1, 0, 0, 1}, {1, 1, 1, 2} };

  // Set number of dimensions and initialise all the working vectors
  SetNumberOfDimensions_(3);
  AccessInitialiseVectors();
  SetVertices(v1);
  SetWorstVertex(3);
  AccessComputeCentroid();

  // Compute the contract point using the inward contract function
  SampleCostFunction5 Func;
  CHECK(AccessContractInside(Func));

  // Initialise an empty vector and get the contracted point and test the values
  std::vector<double> store_contract;
  GetContract(store_contract);
  CHECK_CLOSE((2.0/3.0), store_contract[0], epsilon);
  CHECK_CLOSE(0.5, store_contract[1], epsilon);
  CHECK_CLOSE(0.5, store_contract[2], epsilon);
}

// Test whether Contract_in correctly computes inward contracted point by
// checking the value stored in contract_ for 5D function
TEST_FIXTURE(TestNelderMead, Contract_in2)
{
  // Initialise the population_ with a generated simplex
  std::vector<std::vector<double>> v1 = { {0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0},
      {1, 1, 1, 0, 0, 3}, {1, 1, 1, 1, 0, 4} };

  // Set number of dimensions and initialise all the working vectors
  SetNumberOfDimensions_(5);
  AccessInitialiseVectors();
  SetVertices(v1);
  SetWorstVertex(5);
  AccessComputeCentroid();

  // Compute the contract point using the inward contract function
  SampleCostFunction6 Func;
  CHECK(AccessContractInside(Func));

  // Initialise an empty vector and get the contracted point and test the values
  std::vector<double> store_contract;
  GetContract(store_contract);
  CHECK_CLOSE(0.6, store_contract[0], epsilon);
  CHECK_CLOSE(0.6, store_contract[1], epsilon);
  CHECK_CLOSE(0.6, store_contract[2], epsilon);
  CHECK_CLOSE(0.5, store_contract[3], epsilon);
  CHECK_CLOSE(0, store_contract[4], epsilon);
  CHECK_CLOSE(1.33, store_contract[5], epsilon);
}

TEST_FIXTURE(TestNelderMead, ContractInvalidCost)
{
  // Initialise the population_ with a generated simplex
  std::vector<std::vector<double>> v1 = { {0, 0, 0, 0}, {0, 0, 0, 0},
      {1, 0, 0, 1}, {1, 1, 0, 2} };

  // Set number of dimensions and initialise all the working vectors
  SetNumberOfDimensions_(3);
  AccessInitialiseVectors();
  SetVertices(v1);

  AlwaysInfinite Func;
  // Compute the contract point using the inward contract function
  CHECK(!AccessContractInside(Func));
  // Compute the contract point using the outward contract function
  CHECK(!AccessContractOutside(Func));
}

// Test whether Shrink correctly shrinks the simplex towards the best point and
// updates the vertex by newly computed shrink value by checking the value
// stored in the worst vertex
TEST_FIXTURE(TestNelderMead, Shrink1)
{
  // Initialise the population_ with a generated simplex
  std::vector<std::vector<double>> v1 = {{2, 0, -3.36}, {0, 6, -0.16},
      {0, 4, 0}};
  // Set number of dimensions and initialise all the working vectors
  SetNumberOfDimensions_(2);
  AccessInitialiseVectors();
  SetVertices(v1);

  // Compute all the shrunk points using the shrink function
  SampleCostFunction1 Func;
  CHECK(AccessShrink(Func));

  // Initialise an empty vector and get the shrunk points and test the values
  std::vector<double> v2 = GetSolution(1);
  v2.push_back(GetCost(1));
  CHECK_CLOSE(1, v2[0], epsilon);
  CHECK_CLOSE(3, v2[1], epsilon);
  CHECK_CLOSE(0, v2[2], epsilon);
  v2 = GetSolution(2);
  v2.push_back(GetCost(2));
  CHECK_CLOSE(1, v2[0], epsilon);
  CHECK_CLOSE(2, v2[1], epsilon);
  CHECK_CLOSE(3, v2[2], epsilon);
}

// Test Shrink in 4D by checking the values stored in all the vertices except
// for the best vertex
TEST_FIXTURE(TestNelderMead, Shrink4D)
{
  // Initialise the population_ with a generated simplex
  std::vector<std::vector<double>> v1 = {{1, 1, 0, 0, 5}, {0, 1, 0, 0, 4},
      {0, 0, 1, 0, 3}, {0, 0, 0, 1, 2}, {0, 0, 0, 0, 6}};

  // Set number of dimensions and initialise all the working vectors
  SetNumberOfDimensions_(4);
  AccessInitialiseVectors();
  SetVertices(v1);
  SortPopulation();

  // Compute all the shrunk points using the shrink function
  SampleCostFunction1 Func;
  CHECK(AccessShrink(Func));

  // Initialise an empty vector and get the shrunk points and test the values
  std::vector<double> v2 = GetSolution(1);
  CHECK_CLOSE(0, v2[0], epsilon);
  CHECK_CLOSE(0, v2[1], epsilon);
  CHECK_CLOSE(0.5, v2[2], epsilon);
  CHECK_CLOSE(0.5, v2[3], epsilon);
  v2 = GetSolution(2);
  CHECK_CLOSE(0, v2[0], epsilon);
  CHECK_CLOSE(0.5, v2[1], epsilon);
  CHECK_CLOSE(0, v2[2], epsilon);
  CHECK_CLOSE(0.5, v2[3], epsilon);
  v2 = GetSolution(3);
  CHECK_CLOSE(0.5, v2[0], epsilon);
  CHECK_CLOSE(0.5, v2[1], epsilon);
  CHECK_CLOSE(0, v2[2], epsilon);
  CHECK_CLOSE(0.5, v2[3], epsilon);
  v2 = GetSolution(4);
  CHECK_CLOSE(0, v2[0], epsilon);
  CHECK_CLOSE(0, v2[1], epsilon);
  CHECK_CLOSE(0, v2[2], epsilon);
  CHECK_CLOSE(0.5, v2[3], epsilon);
}

TEST_FIXTURE(TestNelderMead, ShrinkInvalidCost)
{
  // Initialise the population_ with a generated simplex
  std::vector<std::vector<double>> v1 = { {0, 0, 0, 0}, {0, 0, 0, 0},
      {1, 0, 0, 1}, {1, 1, 0, 2} };

  // Set number of dimensions and initialise all the working vectors
  SetNumberOfDimensions_(3);
  AccessInitialiseVectors();
  SetVertices(v1);
  AlwaysInfinite Func;
  CHECK(!AccessShrink(Func));
}

// Test whether Expand correctly expands the simplex by checking the value
// stored in expand_
TEST_FIXTURE(TestNelderMead, Expand)
{
  // Initialise the population_ with a generated simplex
  std::vector<std::vector<double>> v1 = {{2, 4, 2}, {56, 234234, -4},
      {45, 12412, 2}};

  // Set number of dimensions and initialise all the working vectors
  SetNumberOfDimensions_(2);
  AccessInitialiseVectors();
  SetVertices(v1);

  // Compute and set reflect and centroid values
  std::vector<double> reflect_values = {1.2, 0.8, 0};
  std::vector<double> centroid_values = {0.6, 0.4, 0};
  SetReflect(reflect_values);
  SetCentroid(centroid_values);

  // Compute the expanded point using the expand function
  SampleCostFunction1 Func;
  CHECK(AccessExpand(Func));

  // Initialise an empty vector and get the expanded points and test the values
  std::vector<double> expand_values;
  GetExpand(expand_values);
  CHECK_CLOSE(1.8, expand_values[0], epsilon);
  CHECK_CLOSE(1.2, expand_values[1], epsilon);
  CHECK_CLOSE(5.88, expand_values[2], epsilon);
}

// Test with different predefined centroid and reflect values
TEST_FIXTURE(TestNelderMead, Expand2)
{
  // Initialise the population_ with a generated simplex
  std::vector<std::vector<double>> v1 = {{2, 4, 2}, {56, 234234, -4},
      {45, 12412, 2}};

  // Set number of dimensions and initialise all the working vectors
  SetNumberOfDimensions_(2);
  AccessInitialiseVectors();
  SetVertices(v1);

  // Compute and set reflect and centroid values
  std::vector<double> reflect_values = {3.0, 0.4, 0};
  std::vector<double> centroid_values = {1.5, 0.6, 0};
  SetReflect(reflect_values);
  SetCentroid(centroid_values);

  // Compute the expanded point using the expand function
  SampleCostFunction1 Func;
  CHECK(AccessExpand(Func));

  // Initialise an empty vector and get the expanded points and test the values
  std::vector<double> expand_values;
  GetExpand(expand_values);
  CHECK_CLOSE(4.5, expand_values[0], epsilon);
  CHECK_CLOSE(0.2, expand_values[1], epsilon);
  CHECK_CLOSE(1.19, expand_values[2], epsilon);
}

// Test with different predefined centroid and reflect values
TEST_FIXTURE(TestNelderMead, Expand3)
{
  // Initialise the population_ with a generated simplex
  std::vector<std::vector<double>> v1 = {{2, 4, 2}, {56, 234234, -4},
      {45, 12412, 2}};

  // Set number of dimensions and initialise all the working vectors
  SetNumberOfDimensions_(2);
  AccessInitialiseVectors();
  SetVertices(v1);

  // Compute and set reflect and centroid values
  std::vector<double> reflect_values = {3.6, 1.6, 0};
  std::vector<double> centroid_values = {2.4, 0.8, 0};
  SetReflect(reflect_values);
  SetCentroid(centroid_values);

  // Compute the expanded point using the expand function
  SampleCostFunction1 Func;
  CHECK(AccessExpand(Func));

  // Initialise an empty vector and get the expanded points and test the values
  std::vector<double> expand_values;
  GetExpand(expand_values);
  CHECK_CLOSE(4.8, expand_values[0], epsilon);
  CHECK_CLOSE(2.4, expand_values[1], epsilon);
  CHECK_CLOSE(4.32, expand_values[2], epsilon);
}

TEST_FIXTURE(TestNelderMead, ExpandInfinity)
{
  // define a simplex
  std::vector<double> v1 = {7};
  std::vector<double> v2 = {5.5};
  LineToInfinity cost_func;
  std::vector<std::vector<double> > simplex = {v1, v2};
  for (auto &vertex : simplex) {
    auto residuals = cost_func(vertex);
    vertex.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  Unfit::NelderMead nm;
  nm.options.SetMaxIterations(2);
  nm.SetPopulation(simplex);
//  auto rc = nm.FindMin(cost_func, v1);

//  CHECK_EQUAL(2, nm.FindMin(cost_func, simplex));
  CHECK_EQUAL(2, nm.FindMin(cost_func, v1));
  CHECK_CLOSE(4.0, nm.GetCost(0), 1e-4);  // one reflect
}

TEST_FIXTURE(TestNelderMead, ContractInfinity)
{
  std::vector<double> v1 = {3};
  std::vector<double> v2 = {-3};
  std::vector<std::vector<double> > simplex_contarct_out = {v1, v2};
  InfinityInBetween cost_func;
  for (auto &vertex : simplex_contarct_out) {
    auto residuals = cost_func(vertex);
    vertex.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  Unfit::NelderMead nm;
  nm.options.SetMaxIterations(2);
  nm.SetPopulation(simplex_contarct_out);

  CHECK_EQUAL(2, nm.FindMin(cost_func, v1));
  CHECK_CLOSE(9.0, nm.GetCost(0), 1e-4);  // 3^2 no improvement...
}

TEST_FIXTURE(TestNelderMead, ContractOutsideInfiniteCost)
{
  Unfit::NelderMead nm;
  nm.options.SetMaxIterations(2);
  // Initial simplex
  std::vector<std::vector<double>> simplex {{1.0, 1.0}, {2.0, 1.0}, {2.0, 2.0}};
  Infinity2D cost_func;
  for (auto &vertex : simplex) {
    auto residuals = cost_func(vertex);
    vertex.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  // Costs of the initial simplex are (12, 27 & 48)
  // Cost of the reflect point in 36.75 forcing an outside contract
  // Cost of the outside contract is infinite forcing a shrink
  nm.SetPopulation(simplex);
  std::vector<double> vertex {1.0, 1.0};
  CHECK_EQUAL(2, nm.FindMin(cost_func, vertex));
}

// Test Expand in 4D
TEST_FIXTURE(TestNelderMead, MultipleExpand1)
{
  // Initialise the population_ with a generated simplex
  std::vector<std::vector<double>> v1 = {{1, 2, 3, 4, 1}, {5, 6, 7, 8, 4},
      {9, 10, 11, 12, 2}, {13, 14, 15, 16, 2}, {17, 18, 19, 20, 4}};

  // Set number of dimensions and initialise all the working vectors
  SetNumberOfDimensions_(4);
  AccessInitialiseVectors();
  SetVertices(v1);

  // Check whether the vertex has been set properly
  for (unsigned i = 0; i < 5; ++i) {
    std::vector<double> vertex = GetSolution(i);
    for (unsigned j = 0; j < 4; ++j) {
      CHECK_CLOSE(vertex[j], v1[i][j], 1e-16);
    }
  }
  // Compute and set reflect and centroid values
  std::vector<double> v2 = {1, 2, 3, 4, 5};
  SetReflect(v2);
  v2 = {11.25, 12.5, 13.75, 15, 0};
  SetCentroid(v2);

  // Compute the expanded point using the expand function
  PowellSingular func;
  CHECK(AccessExpand(func));

  // Initialise an empty vector and get the expanded points and test the values
  GetExpand(v2);
  CHECK_CLOSE(-9.25 , v2[0], 1e-16);
  CHECK_CLOSE(-8.5, v2[1], 1e-16);
  CHECK_CLOSE(-7.75, v2[2], 1e-16);
  CHECK_CLOSE(-7, v2[3], 1e-16);
}

// Test whether Reflect correctly reflects the simplex by checking the value
// stored in reflect_
TEST_FIXTURE(TestNelderMead, Reflect1)
{
  // Initialise the population_ with a generated simplex
  std::vector<std::vector<double>> v1 = {{1.2, 0, -3.36}, {0, 0.8, -0.16},
      {0, 0, -0.16}};

  // Set number of dimensions and initialise all the working vectors
  SetNumberOfDimensions_(2);
  AccessInitialiseVectors();
  SetVertices(v1);
  SetWorstVertex(2);
  AccessComputeCentroid();

  // Compute reflect with centroid values
  SampleCostFunction1 Func;
  CHECK(AccessReflect(Func));

  // Check computed centroid coordinates
  std::vector<double> test_mid_point;
  GetCentroid(test_mid_point);
  CHECK_CLOSE(0.6, test_mid_point[0], epsilon);
  CHECK_CLOSE(0.4, test_mid_point[1], epsilon);

  // Initialise an empty vector and get the reflected points and test the values
  std::vector<double> test_reflect_point;
  GetReflect(test_reflect_point);
  CHECK_CLOSE(1.2, test_reflect_point[0], epsilon);
  CHECK_CLOSE(0.8, test_reflect_point[1], epsilon);
}

// Test whether Reflect works for 4D simplex
TEST_FIXTURE(TestNelderMead, Reflect4D)
{
  // Initialise the population_ with a generated simplex
  std::vector<std::vector<double>> v1 = {{1, 2, 3, 4, 5}, {6, 7, 8, 9, 10},
      {11, 12, 13, 14, 15}, {16, 17, 18, 19, 20}, {21, 22, 23, 24, 25}};

  // Set number of dimensions and initialise all the working vectors
  SetNumberOfDimensions_(4);
  AccessInitialiseVectors();
  SetVertices(v1);
  SetWorstVertex(4);
  SortPopulation();
  AccessComputeCentroid();

  // Compute reflect with centroid values
  SampleCostFunction1 Func;
  CHECK(AccessReflect(Func));

  // Initialise an empty vector and get the reflected points and test the values
  std::vector<double> test_reflect_point;
  GetReflect(test_reflect_point);
  CHECK_CLOSE(-4, test_reflect_point[0], epsilon);
  CHECK_CLOSE(-3, test_reflect_point[1], epsilon);
  CHECK_CLOSE(-2, test_reflect_point[2], epsilon);
  CHECK_CLOSE(-1, test_reflect_point[3], epsilon);
}

// Test with a new initial 4D simplex
TEST_FIXTURE(TestNelderMead, Reflect4D1)
{
  // Initialise the population_ with a generated simplex
  std::vector<std::vector<double>> v1 = {{1, 1, 0, 0, 5}, {0, 1, 0, 0, 4},
      {0, 0, 1, 0, 3}, {0, 0, 0, 1, 2}, {0, 0, 0, 1, 6}};
  // Set number of dimensions and initialise all the working vectors
  SetNumberOfDimensions_(4);
  AccessInitialiseVectors();
  SetVertices(v1);
  SetWorstVertex(4);
  AccessComputeCentroid();

  // Compute reflect with centroid values
  SampleCostFunction1 Func;
  CHECK(AccessReflect(Func));

  // Initialise an empty vector and get the reflected points and test the values
  std::vector<double> test_reflect_point;
  GetReflect(test_reflect_point);
  CHECK_CLOSE(0.5, test_reflect_point[0], epsilon);
  CHECK_CLOSE(1, test_reflect_point[1], epsilon);
  CHECK_CLOSE(0.5, test_reflect_point[2], epsilon);
  CHECK_CLOSE(-0.5, test_reflect_point[3], epsilon);
}

// Test the overload operator of the generic functions check if it returns
// the cost
TEST_FIXTURE(TestNelderMead, Overload1)
{
  // Initialise initial point to be passed into the function
  std::vector<double> initialpoint = {-2, -2, 1};

  // Check that the overloaded operator computes the cost correctly
  SampleCostFunction5 Func;
  CHECK_CLOSE(3, Func(initialpoint)[0], 1e-6);
}

// Check if IsAreaConverged returns the right return codes.
// If vertices are far apart (i.e. area of triangle is big), return 0
TEST_FIXTURE(TestNelderMead, IsConverged1)
{
  // Initialise and set a 2D vector that has not converged
  std::vector<std::vector<double>> v1 = {{2, 1, 9.5}, {3, 1, 111},
      {10, 10, 10}};
  SetVertices(v1);

  // Check that IsConverged checks the vector correctly
  CHECK(!IsConverged(v1[0]));
}

// Check if IsAreaConverged returns the right return codes.
// If vertices are far apart (i.e. area of triangle is big), return 0
TEST_FIXTURE(TestNelderMead, IsConverged2)
{
  // Initialise and set a 2D vector that has not converged
  std::vector<std::vector<double>> v1 = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
  SetVertices(v1);

  // Check that if vertices are far apart, return 0
  CHECK(IsConverged(v1[0]));
}

// Test whether the default constructor contains the correct values
TEST_FIXTURE(TestNelderMead, DefaultConstructor)
{
  CHECK_EQUAL(10000u, options.GetMaxIterations());
  CHECK_EQUAL(100000u, options.GetMaxFunctionEvaluations());
  CHECK_CLOSE(1e-12, options.GetCostTolerance(), 1e-16);
  CHECK_CLOSE(1e-4, options.GetGeometricTolerance(), 1e-16);
  CHECK_CLOSE(1e-8, options.GetDegenerateTolerance(), 1e-16);
  CHECK_EQUAL(0u, options.GetOutputLevel());
  CHECK_EQUAL(0u, GetNumberOfDimensions_());
  CHECK_EQUAL(0u, GetCurrentFunctionEvaluation());
  CHECK_EQUAL(0u, GetCurrentIterations());
  CHECK(!options.GetUseAdaptiveParameters());
  double alpha, beta, delta, gamma;
  options.GetNelderMeadStepSizes(alpha, beta, delta, gamma);
  CHECK_CLOSE(1, alpha, 1e-16);
  CHECK_CLOSE(2, beta, 1e-16);
  CHECK_CLOSE(0.5, delta, 1e-16);
  CHECK_CLOSE(0.5, gamma, 1e-16);
}

// Test whether InitialiseWorkingVector initialise vectors in the NelderMead
// class to the correct size
TEST_FIXTURE(TestNelderMead, TestInitialiseVectors)
{
  // Set number of dimensions and initialise working vector
  SetNumberOfDimensions_(3);
  AccessInitialiseVectors();

  // Initialise empty vectors and get from the private members and check size
  std::vector<double> contract;
  GetContract(contract);
  CHECK_EQUAL(4u, contract.size());
  std::vector<double> reflect;
  GetContract(reflect);
  CHECK_EQUAL(4u, reflect.size());
  std::vector<double> expand;
  GetContract(expand);
  CHECK_EQUAL(4u, expand.size());
  std::vector<double> centroid;
  GetContract(centroid);
  CHECK_EQUAL(4u, centroid.size());
}

// Check return code of findmin when the NelderMead encounter insufficient
// iterations
TEST_FIXTURE(TestNelderMead, MaxIterationsReachedOutput3)
{
  // Initialise initial point
  std::vector<double> v2 = {5.0, 4.0, 3.0};

  // Set number of dimensions and initialise working vector
  SetNumberOfDimensions_(2);
  AccessInitialiseVectors();

  // Setting max iteration and output level
  options.SetMaxIterations(3);
  options.SetOutputLevel(0);
  SampleCostFunction5 Func;

  // Check that max iteration has been set correctly
  CHECK_EQUAL(3u, options.GetMaxIterations());

  // Minimise and check return code
  CHECK_EQUAL(2, FindMin(Func, v2));
}

// Check return code of findmin when the NelderMead encounter insufficent
// function evaluations
TEST_FIXTURE(TestNelderMead, MaxFunctionEvaluationsReachedOutput3)
{
  // Initialise initial point
  std::vector<double> v2 = {0, 0, 0};

  // Set number of dimensions and initialise working vector
  SetNumberOfDimensions_(2);
  AccessInitialiseVectors();

  // Setting max function evaluation and output level
  options.SetMaxFunctionEvaluations(5);
  options.SetOutputLevel(0);
  SampleCostFunction5 Func;

  // Check that max function evaluations has been set correctly
  CHECK_EQUAL(5u, options.GetMaxFunctionEvaluations());

  // Minimize and check return code
  CHECK_EQUAL(1, FindMin(Func, v2));
}

// Test for invalid initial point due to boundary
TEST_FIXTURE(TestNelderMead, GeneratePopulation0)
{
  // Initialise initial point
  std::vector<double> vertex0= {1, 0, 0};

  // Set number of dimensions and initialise working vectors
  SetNumberOfDimensions_(3);
  AccessInitialiseVectors();
  SampleCostFunction5 Func;

  // Set boundary so that initial point is invalid
  bounds.SetBounds(0, 10, 13);

  // Calls the function GeneratePopulation and check the return code
  CHECK_EQUAL(1, AccessGeneratePopulation(Func, vertex0));
}

// Test for invalid initial point due to invalid cost
TEST_FIXTURE(TestNelderMead, GeneratePopulation1)
{
  // Initialise initial point
  std::vector<double> vertex0= {1.05, 1, 1};

  // Set number of dimensions and initialise working vectors
  SetNumberOfDimensions_(3);
  AccessInitialiseVectors();
  SampleCostFunction7 Func;

  // Set boundary so that initial point is invalid
  bounds.SetBounds(0, -1, 2);

  // Calls the function GeneratePopulation and check the return code
  CHECK_EQUAL(1, AccessGeneratePopulation(Func, vertex0));
}

// case 2 & positive coordinate & valid costs
TEST_FIXTURE(TestNelderMead, GeneratePopulation2)
{
  // Initialise initial point
  std::vector<double> vertex0= {0, 1, 0};

  // Set number of dimensions and initialise working vectors
  SetNumberOfDimensions_(3);
  AccessInitialiseVectors();
  SampleCostFunction5 Func;

  // Set boundaries
  bounds.SetBounds(0, -1, 10);
  bounds.SetBounds(1, -1, 10);

  // Calls the function GeneratePopulation
  AccessGeneratePopulation(Func, vertex0);

  // Initialise vectors and get the vertices from the NelderMead class
  vertex0 = GetSolution(0);
  vertex0.push_back(GetCost(0));
  std::vector<double> vertex1 = GetSolution(1);
  vertex1.push_back(GetCost(1));
  std::vector<double> vertex2 = GetSolution(2);
  vertex2.push_back(GetCost(2));
  std::vector<double> vertex3 = GetSolution(3);
  vertex3.push_back(GetCost(3));
  // vertex0
  CHECK_CLOSE(0, vertex0[0], epsilon);
  CHECK_CLOSE(1, vertex0[1], epsilon);
  CHECK_CLOSE(0, vertex0[2], epsilon);
  CHECK_CLOSE(1, vertex0[3], epsilon);
  // vertex1
  CHECK_CLOSE(0.00025, vertex1[0], epsilon);
  CHECK_CLOSE(1, vertex1[1], epsilon);
  CHECK_CLOSE(0, vertex1[2], epsilon);
  CHECK_CLOSE(1, vertex1[3], 1.0e-6);
  // vertex2
  CHECK_CLOSE(0, vertex2[0], epsilon);
  CHECK_CLOSE(1.05, vertex2[1], epsilon);
  CHECK_CLOSE(0, vertex2[2], epsilon);
  CHECK_CLOSE(1.1025, vertex2[3], epsilon);
  // vertex3
  CHECK_CLOSE(0, vertex3[0], epsilon);
  CHECK_CLOSE(1, vertex3[1], epsilon);
  CHECK_CLOSE(0.00025, vertex3[2], epsilon);
  CHECK_CLOSE(1, vertex3[3], 1.0e-6);
}

// case 2 & negative coordinate & valid costs
TEST_FIXTURE(TestNelderMead, GeneratePopulation3)
{
  // Initialise initial point
  std::vector<double> vertex0= {0, -1, 0};

  // Set number of dimensions and initialise working vectors
  SetNumberOfDimensions_(3);
  AccessInitialiseVectors();
  SampleCostFunction5 Func;

  // Set boundaries
  bounds.SetBounds(1, -5, 10);

  // Calls the function GeneratePopulation
  AccessGeneratePopulation(Func, vertex0);

  // Initialise vectors and get the vertices from the NelderMead class
  vertex0 = GetSolution(2);
  vertex0.push_back(GetCost(2));
  // vertex0
  CHECK_CLOSE(0, vertex0[0], epsilon);
  CHECK_CLOSE(-1.05, vertex0[1], epsilon);
  CHECK_CLOSE(0, vertex0[2], epsilon);
  CHECK_CLOSE(1.1025, vertex0[3], epsilon);
}

// Case1 & valid cost
TEST_FIXTURE(TestNelderMead, GeneratePopulation4)
{
  // Initialise initial point
  std::vector<double> vertex0= {1, 0, 0};

  // Set number of dimensions and initialise working vectors
  SetNumberOfDimensions_(3);
  AccessInitialiseVectors();
  SampleCostFunction5 Func;

  // Set boundary and create the simplex
  bounds.SetBounds(0, -10, 1);
  AccessGeneratePopulation(Func, vertex0);

  // Get the vertices out from the NelderMead class
  vertex0 = GetSolution(1);
  vertex0.push_back(GetCost(1));
  // vertex1
  CHECK_CLOSE(0.95, vertex0[0], epsilon);
  CHECK_CLOSE(0, vertex0[1], epsilon);
  CHECK_CLOSE(0, vertex0[2], epsilon);
  CHECK_CLOSE(0.9025, vertex0[3], epsilon);
}

// negative coordinate & Case1 & valid cost
TEST_FIXTURE(TestNelderMead, GeneratePopulation5)
{
  // Initial guess
  std::vector<double> vertex0= {-1, 0, 0};

  SampleCostFunction5 Func;
  SetNumberOfDimensions_(3);
  AccessInitialiseVectors();
  bounds.SetBounds(0, -1, 1);
  AccessGeneratePopulation(Func, vertex0);
  vertex0 = GetSolution(1);
  vertex0.push_back(GetCost(1));

  // vertex1
  // Check the result matches what we expect
  CHECK_CLOSE(-0.95, vertex0[0], epsilon);
  CHECK_CLOSE(0, vertex0[1], epsilon);
  CHECK_CLOSE(0, vertex0[2], epsilon);
  CHECK_CLOSE(0.9025, vertex0[3], epsilon);
}

// Case2 & zero coordinate & valid cost
TEST_FIXTURE(TestNelderMead, GeneratePopulation6)
{
  // Initialise initial point
  std::vector<double> vertex0= {0, 0, 0};

  // Set number of dimensions and initialise working vectors
  SetNumberOfDimensions_(3);
  AccessInitialiseVectors();
  SampleCostFunction5 Func;

  // Set boundary and create the simplex
  bounds.SetBounds(0, -10, 3);
  AccessGeneratePopulation(Func, vertex0);

  // Get the vertices out from the NelderMead class
  vertex0 = GetSolution(1);
  vertex0.push_back(GetCost(1));
  // vertex1
  CHECK_CLOSE(0.00025, vertex0[0], epsilon);
  CHECK_CLOSE(0, vertex0[1], epsilon);
  CHECK_CLOSE(0, vertex0[2], epsilon);
  CHECK_CLOSE(6.25e-8, vertex0[3], epsilon);
}

// Case1 & zero coordinate & valid cost
TEST_FIXTURE(TestNelderMead, GeneratePopulation7)
{
  // Initialise initial point
  std::vector<double> vertex0= {0, 0, 0};

  // Set number of dimensions and initialise working vectors
  SetNumberOfDimensions_(3);
  AccessInitialiseVectors();
  SampleCostFunction5 Func;

  // Set boundary and create the simplex
  bounds.SetBounds(0, -10, 0);
  AccessGeneratePopulation(Func, vertex0);

  // Get the vertices out from the NelderMead class
  vertex0 = GetSolution(1);
  vertex0.push_back(GetCost(1));
  // vertex1
  CHECK_CLOSE(-0.00025, vertex0[0], epsilon);
  CHECK_CLOSE(0, vertex0[1], epsilon);
  CHECK_CLOSE(0, vertex0[2], epsilon);
  CHECK_CLOSE(6.25e-8, vertex0[3], epsilon);
}

// Case2 & invalid cost
TEST_FIXTURE(TestNelderMead, GeneratePopulation8)
{
  // Initialise initial point
  std::vector<double> vertex0= {1, 1, 1};

  // Set number of dimensions and initialise working vectors
  SetNumberOfDimensions_(3);
  AccessInitialiseVectors();
  SampleCostFunction7 Func;

  // Set boundary and create the simplex
  bounds.SetBounds(0, -10, 3);
  int rc = AccessGeneratePopulation(Func, vertex0);

  // The resulting simplex has invalid costs, so the return code is 3
  CHECK_EQUAL(3, rc);
}

// Case1  & invalid cost
TEST_FIXTURE(TestNelderMead, GeneratePopulation9)
{
  // Initialise initial point
  std::vector<double> vertex0= {1, 1, 1};

  // Set number of dimensions and initialise working vectors
  SetNumberOfDimensions_(3);
  AccessInitialiseVectors();
  SampleCostFunction9 Func;

  // Set boundary and create the simplex
  bounds.SetBounds(0, -10, 1);
  int rc = AccessGeneratePopulation(Func, vertex0);

  // The resulting simplex has invalid costs, so the return code is 3
  CHECK_EQUAL(3, rc);
}

TEST_FIXTURE(TestNelderMead, GeneratePopulationCloseToBounds)
{
  // Initialise initial point
  std::vector<double> vertex0 = {1, 1};
  std::vector<double> vertex1 = {1, 1};

  // Set number of dimensions and initialise working vectors
  SetNumberOfDimensions_(2);
  AccessInitialiseVectors();
  SampleCostFunction1 cost_func;

  // Set boundary between nonzero_scale and zero_scale above
  bool bounds_set = bounds.SetBounds(0, -10, 1.01);
  CHECK(bounds_set);
  int rc = AccessGeneratePopulation(cost_func, vertex0);
  CHECK_EQUAL(0, rc);
  vertex1 = GetSolution(1);
  CHECK_CLOSE(0.95, vertex1[0], epsilon);
  CHECK_CLOSE(1, vertex1[1], epsilon);

  // Set boundary between nonzero_scale and zero_scale below
  bounds_set = bounds.SetBounds(0, 0.99, 10);
  CHECK(bounds_set);
  rc = AccessGeneratePopulation(cost_func, vertex0);
  CHECK_EQUAL(0, rc);
  vertex1 = GetSolution(1);
  CHECK_CLOSE(1.05, vertex1[0], epsilon);
  CHECK_CLOSE(1, vertex1[1], epsilon);

  // Set boundary between nonzero_scale and zero_scale above and below
  bounds_set = bounds.SetBounds(0, 0.99, 1.01);
  CHECK(bounds_set);
  rc = AccessGeneratePopulation(cost_func, vertex0);
  CHECK_EQUAL(0, rc);
  vertex1 = GetSolution(1);
  CHECK_CLOSE(1.00025, vertex1[0], epsilon);
  CHECK_CLOSE(1, vertex1[1], epsilon);

  // Set boundary closer than zero_scale above
  bounds_set = bounds.SetBounds(0, -10, 1.0001);
  CHECK(bounds_set);
  rc = AccessGeneratePopulation(cost_func, vertex0);
  CHECK_EQUAL(0, rc);
  vertex1 = GetSolution(1);
  CHECK_CLOSE(0.95, vertex1[0], epsilon);
  CHECK_CLOSE(1, vertex1[1], epsilon);

  // Set boundary closer than zero_scale below
  bounds_set = bounds.SetBounds(0, 0.9999, 10);
  CHECK(bounds_set);
  rc = AccessGeneratePopulation(cost_func, vertex0);
  CHECK_EQUAL(0, rc);
  vertex1 = GetSolution(1);
  CHECK_CLOSE(1.05, vertex1[0], epsilon);
  CHECK_CLOSE(1, vertex1[1], epsilon);

  // Set boundary closer than zero_scale above and below
  bounds_set = bounds.SetBounds(0, 0.9999, 1.0001);
  CHECK(bounds_set);
  rc = AccessGeneratePopulation(cost_func, vertex0);
  CHECK_EQUAL(2, rc);  // bounds are too close together

  // Set boundary closer than zero_scale above and nonzero_scale below
  bounds_set = bounds.SetBounds(0, 0.99, 1.0001);
  CHECK(bounds_set);
  rc = AccessGeneratePopulation(cost_func, vertex0);
  CHECK_EQUAL(0, rc);
  vertex1 = GetSolution(1);
  CHECK_CLOSE(0.99975, vertex1[0], epsilon);
  CHECK_CLOSE(1, vertex1[1], epsilon);

  // Set boundary closer than nonzero_scale above and zero_scale below
  bounds_set = bounds.SetBounds(0, 0.9999, 1.01);
  CHECK(bounds_set);
  rc = AccessGeneratePopulation(cost_func, vertex0);
  CHECK_EQUAL(0, rc);
  vertex1 = GetSolution(1);
  CHECK_CLOSE(1.00025, vertex1[0], epsilon);
  CHECK_CLOSE(1, vertex1[1], epsilon);

  // Now choose a negative starting point
  vertex0 = {-1, -1};

  // Set boundary between nonzero_scale and zero_scale above
  bounds_set = bounds.SetBounds(0, -10, -0.99);
  CHECK(bounds_set);
  rc = AccessGeneratePopulation(cost_func, vertex0);
  CHECK_EQUAL(0, rc);
  vertex1 = GetSolution(1);
  CHECK_CLOSE(-1.05, vertex1[0], epsilon);
  CHECK_CLOSE(-1, vertex1[1], epsilon);

  // Set boundary between nonzero_scale and zero_scale below
  bounds_set = bounds.SetBounds(0, -1.01, 10);
  CHECK(bounds_set);
  rc = AccessGeneratePopulation(cost_func, vertex0);
  CHECK_EQUAL(0, rc);
  vertex1 = GetSolution(1);
  CHECK_CLOSE(-0.95, vertex1[0], epsilon);
  CHECK_CLOSE(-1, vertex1[1], epsilon);

  // Set boundary between nonzero_scale and zero_scale above and below
  bounds_set = bounds.SetBounds(0, -1.01, -0.99);
  CHECK(bounds_set);
  rc = AccessGeneratePopulation(cost_func, vertex0);
  CHECK_EQUAL(0, rc);
  vertex1 = GetSolution(1);
  CHECK_CLOSE(-0.99975, vertex1[0], epsilon);
  CHECK_CLOSE(-1, vertex1[1], epsilon);

  // Set boundary closer than zero_scale above
  bounds_set = bounds.SetBounds(0, -10, -0.9999);
  CHECK(bounds_set);
  rc = AccessGeneratePopulation(cost_func, vertex0);
  CHECK_EQUAL(0, rc);
  vertex1 = GetSolution(1);
  CHECK_CLOSE(-1.05, vertex1[0], epsilon);
  CHECK_CLOSE(-1, vertex1[1], epsilon);

  // Set boundary closer than zero_scale below
  bounds_set = bounds.SetBounds(0, -1.0001, 10);
  CHECK(bounds_set);
  rc = AccessGeneratePopulation(cost_func, vertex0);
  CHECK_EQUAL(0, rc);
  vertex1 = GetSolution(1);
  CHECK_CLOSE(-0.95, vertex1[0], epsilon);
  CHECK_CLOSE(-1, vertex1[1], epsilon);

  // Set boundary closer than zero_scale above and below
  bounds_set = bounds.SetBounds(0, -1.0001, -0.9999);
  CHECK(bounds_set);
  rc = AccessGeneratePopulation(cost_func, vertex0);
  CHECK_EQUAL(2, rc);  // bounds are too close together

  // Set boundary closer than zero_scale above and nonzero_scale below
  bounds_set = bounds.SetBounds(0, -1.01, -0.9999);
  CHECK(bounds_set);
  rc = AccessGeneratePopulation(cost_func, vertex0);
  CHECK_EQUAL(0, rc);
  vertex1 = GetSolution(1);
  CHECK_CLOSE(-1.00025, vertex1[0], epsilon);
  CHECK_CLOSE(-1, vertex1[1], epsilon);

  // Set boundary closer than nonzero_scale above and zero_scale below
  bounds_set = bounds.SetBounds(0, -1.0001, -0.99);
  CHECK(bounds_set);
  rc = AccessGeneratePopulation(cost_func, vertex0);
  CHECK_EQUAL(0, rc);
  vertex1 = GetSolution(1);
  CHECK_CLOSE(-0.99975, vertex1[0], epsilon);
  CHECK_CLOSE(-1, vertex1[1], epsilon);

  // Now choose a zero starting point
  vertex0 = {0, 0};

  // Set boundary closer than zero_scale above
  bounds_set = bounds.SetBounds(0, -10, 0.0001);
  CHECK(bounds_set);
  rc = AccessGeneratePopulation(cost_func, vertex0);
  CHECK_EQUAL(0, rc);
  vertex1 = GetSolution(1);
  CHECK_CLOSE(-0.00025, vertex1[0], epsilon);
  CHECK_CLOSE(0, vertex1[1], epsilon);

  // Set boundary closer than zero_scale below
  bounds_set = bounds.SetBounds(0, -0.0001, 10);
  CHECK(bounds_set);
  rc = AccessGeneratePopulation(cost_func, vertex0);
  CHECK_EQUAL(0, rc);
  vertex1 = GetSolution(1);
  CHECK_CLOSE(0.00025, vertex1[0], epsilon);
  CHECK_CLOSE(0, vertex1[1], epsilon);

  // Set boundary closer than zero_scale above and below
  bounds_set = bounds.SetBounds(0, -0.0001, 0.0001);
  CHECK(bounds_set);
  rc = AccessGeneratePopulation(cost_func, vertex0);
  CHECK_EQUAL(2, rc);  // bounds are too close together
}

// Check that 2D asymptote will not find a minimum point
TEST_FIXTURE(TestNelderMead, Asymptote2Dtest)
{
  // Initialise initial point and asymptote function
  Asymptote2D Func;
  std::vector<double> min_point = {1, 1};

  // Check that operator returns correct cost function
  CHECK_CLOSE(1, Func(min_point)[0], 1e-4);

  options.SetMaxIterations(10);

  // Check that max iterations have been reached
  CHECK_EQUAL(2, FindMin(Func, min_point));
}

TEST_FIXTURE(TestNelderMead, CheckVariousOutPutLevel)
{
  SampleCostFunction1 Func;
  // To cover for the printoutput
  std::vector<double> min_point = {-1.2, -10};
  Reset();
  options.SetMaxIterations(3);
  // Cover for outputlevel 1
  options.SetOutputLevel(1);
  FindMin(Func, min_point);
  // Reset
  Reset();
  options.SetMaxIterations(3);
  // Cover for outputlevel 2
  options.SetOutputLevel(2);
  FindMin(Func, min_point);
  // Reset
  Reset();
  options.SetMaxIterations(100);
  // Cover for outputlevel 3
  options.SetOutputLevel(3);
  FindMin(Func, min_point);
}

// This test write information to the screen, which means that the output is not
// tested by unittest++. However, it will fail if anything unexpected happens
// (e.g. crash).
TEST_FIXTURE(TestNelderMead, PrintIterationOutputProcesses)
{
  SampleCostFunction1 func;
  std::vector<double> min_point = {-1.2, -10};
  FindMin(func, min_point);
  options.SetOutputLevel(2);
  const double dummy_cost = 0.0;

  // Reflect
  SetProcess(0);
  AccessPrintIterationOutput(dummy_cost);

  // Expand
  SetProcess(1);
  AccessPrintIterationOutput(dummy_cost);

  // Contract Inside
  SetProcess(2);
  AccessPrintIterationOutput(dummy_cost);

  // Contract Outside
  SetProcess(3);
  AccessPrintIterationOutput(dummy_cost);

  // Shrink
  SetProcess(4);
  AccessPrintIterationOutput(dummy_cost);

  // Restart
  SetProcess(5);
  AccessPrintIterationOutput(dummy_cost);

  // Undefined
  SetProcess(12);
  AccessPrintIterationOutput(dummy_cost);
}

// Check that RegeneratePopulation produces the correct simplex
TEST_FIXTURE(TestNelderMead, RegeneratePopulation1)
{
  // Create vectors containing the vertices and their cost
  std::vector<double> vertex0= {0, 0, 0};
  std::vector<double> vertex1 = {1, 0, 0};
  std::vector<double> vertex2 = {0, 1, 0};
  std::vector<double> vertex3 = {0, 0, 1};
  std::vector<std::vector<double>> vertices = {vertex0, vertex1, vertex2,
      vertex3};

  // Initialise cost function and InitialiseVectors
  SampleCostFunction5 Func;
  SetNumberOfDimensions_(3);
  AccessInitialiseVectors();

  // Assign the simplex
  for (auto &vertex : vertices) {
    auto residuals = Func(vertex);
    vertex.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  SetPopulation(vertices);

  // Regenerate the simplex using the best point
  AccessRegeneratePopulation(Func);

  // Get the vertices of the regenerated simplex and check
  vertex0 = GetSolution(0);
  vertex0.push_back(GetCost(0));
  vertex1 = GetSolution(1);
  vertex1.push_back(GetCost(1));
  vertex2 = GetSolution(2);
  vertex2.push_back(GetCost(2));
  vertex3 = GetSolution(3);
  vertex3.push_back(GetCost(3));
  // vertex0
  CHECK_CLOSE(0, vertex0[0], epsilon);
  CHECK_CLOSE(0, vertex0[1], epsilon);
  CHECK_CLOSE(0, vertex0[2], epsilon);
  CHECK_CLOSE(0, vertex0[3], epsilon);
  // vertex1
  CHECK_CLOSE(0.00025, vertex1[0], epsilon);
  CHECK_CLOSE(0, vertex1[1], epsilon);
  CHECK_CLOSE(0, vertex1[2], epsilon);
  CHECK_CLOSE(6.25e-8, vertex1[3], epsilon);
  // vertex2
  CHECK_CLOSE(0, vertex2[0], epsilon);
  CHECK_CLOSE(0.00025, vertex2[1], epsilon);
  CHECK_CLOSE(0, vertex2[2], epsilon);
  CHECK_CLOSE(6.25e-8, vertex2[3], epsilon);
  // vertex3
  CHECK_CLOSE(0, vertex3[0], epsilon);
  CHECK_CLOSE(0, vertex3[1], epsilon);
  CHECK_CLOSE(0.00025, vertex3[2], epsilon);
  CHECK_CLOSE(6.25e-8, vertex3[3], epsilon);
}

// Check if the function re-produces simplex which previously has invalid cost
TEST_FIXTURE(TestNelderMead, RegeneratePopulation2)
{
  // Create vectors containing the vertices and their cost
  std::vector<double> vertex0= {0, 1, 1};
  std::vector<double> vertex1 = {4.2, 2, 1};
  std::vector<double> vertex2 = {1, 1, 3};
  std::vector<double> vertex3 = {1, 2, 2};
  std::vector<std::vector<double>> vertices = {vertex0, vertex1, vertex2,
      vertex3};

  // Initialise cost function and InitialiseVectors
  SampleCostFunction7 Func;
  SetNumberOfDimensions_(3);
  AccessInitialiseVectors();

  // Assign the simplex
  for (auto &vertex : vertices) {
    auto residuals = Func(vertex);
    vertex.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  SetPopulation(vertices);

  // Regenerate the simplex using the best point
  CHECK_EQUAL(0, AccessRegeneratePopulation(Func));

  // Get the vertices of the regenerated simplex and check
  vertex0 = GetSolution(1);
  vertex0.push_back(GetCost(1));

  // vertex0
  CHECK_CLOSE(0.00025, vertex0[0], epsilon);
  CHECK_CLOSE(1, vertex0[1], epsilon);
  CHECK_CLOSE(1, vertex0[2], epsilon);
  CHECK_CLOSE(0.952608, vertex0[3], 1e-4);
}

// Check whether Degenerate detects the degenerate case correctly
TEST_FIXTURE(TestNelderMead, IsDegenerateCoordinateCheck)
{
  // Generate a simplex with cost and set it into the NelderMead Class
  std::vector<std::vector<double>> v = {{2, 3, 1}, {2, 1, 2}, {2, 2, 3}};
  SetVertices(v);
  SetWorstVertex(2);

  // Compute centroid
  AccessComputeCentroid();

  // Check whether the case has been degenerated
  CHECK_EQUAL(true, AccessIsDegenerate());
}

// Test the function GetPopulation
TEST_FIXTURE(TestNelderMead, GetPopulation)
{
  // Create a simplex and set it to the population_ in NelderMead class
  std::vector<std::vector<double>> v = {{2, 3, 1}, {4, 5, 6}, {7, 8, 9}};
  SetVertices(v);

  // Create a empty vector and calls GetFinalSimplex
  auto coordinates = GetPopulation();

  // Check the values that are retrieved against the values that is set
  // previously
  CHECK_CLOSE(2, coordinates[0][0], 1e-16);
  CHECK_CLOSE(3, coordinates[0][1], 1e-16);
  CHECK_CLOSE(1, coordinates[0][2], 1e-16);
  CHECK_CLOSE(4, coordinates[1][0], 1e-16);
  CHECK_CLOSE(5, coordinates[1][1], 1e-16);
  CHECK_CLOSE(6, coordinates[1][2], 1e-16);
  CHECK_CLOSE(7, coordinates[2][0], 1e-16);
  CHECK_CLOSE(8, coordinates[2][1], 1e-16);
  CHECK_CLOSE(9, coordinates[2][2], 1e-16);
}

// Test the various behaviours of FindMin which accepts a vector of vector
TEST_FIXTURE(TestNelderMead, FindMinWithAssignSimplexBehaviourTest)
{
  std::vector<std::vector<double>> simplex = {{1, 1, 1}, {2, 2, 2}, {5, 6, 7},
    {11, 18, 21}};
  SampleCostFunction5 Func;
  for (auto &vertex : simplex) {
    auto residuals = Func(vertex);
    vertex.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  // Valid simplex
  options.SetMaxIterations(0);  // Check that simplex has been assigned properly
  SetPopulation(simplex);
  std::vector<double> vertex1 {1, 1, 1};
  CHECK_EQUAL(2, FindMin(Func, vertex1));
  CHECK_EQUAL(3u, GetNumberOfDimensions_());
  auto finalsimplex = GetPopulation();
  CHECK_CLOSE(simplex[0][3], finalsimplex[0][3], 1e-16);
  CHECK_CLOSE(simplex[1][3], finalsimplex[1][3], 1e-16);
  CHECK_CLOSE(simplex[2][3], finalsimplex[2][3], 1e-16);
  CHECK_CLOSE(simplex[3][3], finalsimplex[3][3], 1e-16);
  Reset();

  // No of columns greater than number of rows
  simplex.erase(simplex.begin() + 2);
  options.SetMaxIterations(0);  // Check that simplex has been assigned properly
  SetPopulation(simplex);
  CHECK_EQUAL(-3, FindMin(Func, vertex1));
  Reset();

  // No of rows is greater than the number of columns by 2
  simplex = {{1, 1, 1}, {2, 2, 2}, {5, 6, 7}, {11, 18, 19}, {3, 4, 5}};
  options.SetMaxIterations(0);  // Check that simplex has been assigned properly
  for (auto &vertex : simplex) {
    auto residuals = Func(vertex);
    vertex.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  SetPopulation(simplex);
  CHECK_EQUAL(-3, FindMin(Func, vertex1));
  Reset();

  // No of columns is inconsistent for each row
  simplex = {{1, 1, 1}, {2, 2, 2}, {5, 6, 7}, {11, 18, 21}};
  for (auto &vertex : simplex) {
    auto residuals = Func(vertex);
    vertex.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  simplex[2].erase(simplex[2].begin() + 2);
  CHECK_EQUAL(3u, simplex[2].size());
  options.SetMaxIterations(0);  // Check that simplex has been assigned properly
  SetPopulation(simplex);
  CHECK_EQUAL(-4, FindMin(Func, vertex1));
  Reset();

  // Empty simplex
  simplex.clear();
  options.SetMaxIterations(0);  // Check that simplex has been assigned properly
  SetPopulation(simplex);
  CHECK_EQUAL(-1, FindMin(Func, vertex1));
}

TEST_FIXTURE(TestNelderMead, EdgeCasesMaxFunctionAndIteration)
{
  // Check that the consturctor has constructed properly
  CHECK_EQUAL(10000u, options.GetMaxIterations());
  CHECK_EQUAL(100000u, options.GetMaxFunctionEvaluations());

  // Check behaviour if user keyed in a number more than INT_MAX
  options.SetMaxIterations(std::numeric_limits<unsigned>::max()+3);
  options.SetMaxFunctionEvaluations(std::numeric_limits<unsigned>::max()+1);
  CHECK_EQUAL(2u, options.GetMaxIterations());
  CHECK_EQUAL(0u, options.GetMaxFunctionEvaluations());

// Actually this test is not useful as you cannot legally negate an unsigned
//  // Check behaviour if user keyed in a negative number
//  options.SetMaxIterations(-std::numeric_limits<unsigned>::max());
//  options.SetMaxFunctionEvaluations(-std::numeric_limits<unsigned>::max()+2);
//  CHECK_EQUAL(1u, options.GetMaxIterations());
//  CHECK_EQUAL(3u, options.GetMaxFunctionEvaluations());

// Actually this test is not useful/needed as the truncation rules are clear
//  // Check behaviour if user keyed in not an integer
//  options.SetMaxIterations(59.5);
//  options.SetMaxFunctionEvaluations(
//      -std::numeric_limits<unsigned>::max()+0.1);
//  CHECK_EQUAL(59u, options.GetMaxIterations());
//  CHECK_EQUAL(1u, options.GetMaxFunctionEvaluations());

  /* We cannot predict what the number will become when a large number beyond
   the max of an unsigned is entered or when a negative number is keyed in. */
}

TEST(NelderMead_FailFindMin)
{
  Unfit::NelderMead object;

  FailFindMin cost_func;

  // Initial guess
  std::vector<double> min_point = {1.0, 1.0};

  // Check the function calculates the correct cost at the initial point
  CHECK_EQUAL(0, object.FindMin(cost_func, min_point));
  CHECK_CLOSE(0, min_point[0], object.options.GetGeometricTolerance());
  CHECK_CLOSE(0, min_point[1], object.options.GetGeometricTolerance());
}

TEST(NelderMead_FailFindMin2)
{
  Unfit::NelderMead object;
  FailFindMin2 cost_func;
  object.options.SetCostTolerance(1.0e-24);

  // Initial guess
  std::vector<double> min_point = {100.0, -100.0};

  // Minimise
  int rc = object.FindMin(cost_func, min_point);

  // Check the result matches what we expect
  CHECK_EQUAL(4, rc);
  CHECK_CLOSE(90.1, min_point[0], 1e-1);
  CHECK_CLOSE(-50, min_point[1], 1);
}

// TEST(FailFindMin3)
// {
//  Unfit::NelderMead object;
//
//  FailFindMin3 cost_func;
//  object.SetOutputLevel(3);
//
//  // Initial guess
//  std::vector<double> min_point = {1.0, 1.0};
//
//  // Check the function calculates the correct cost at the initial point
//  CHECK_EQUAL(0, object.FindMin(cost_func, min_point));
//  CHECK_CLOSE(-0.26580459, min_point[0], object.GetGeometricTolerance());
//  CHECK_CLOSE(0.39693543, min_point[1], object.GetGeometricTolerance());
// }

// tests collinearity of points with worst point far from the centroid
TEST(NelderMead_FindMinDegenAssignAlmostCollinearBGW)
{
  Unfit::NelderMead object;
  std::vector<std::vector<double>> coordinates = {{0, 12}, {1, 13},
      {2, 14.0001}};
  double tolerance = 1e-4;
  SampleCostFunction3 cost_func;
  for (auto &vertex : coordinates) {
    auto residuals = cost_func(vertex);
    vertex.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  object.SetPopulation(coordinates);
  std::vector<double> coordinate {0, 12};
  int rc = object.FindMin(cost_func, coordinate);
  CHECK_EQUAL(0, rc);
  coordinates = object.GetPopulation();
  CHECK_CLOSE(0, coordinates[0][0], tolerance);
  CHECK_CLOSE(0, coordinates[0][1], tolerance);
  CHECK_CLOSE(0, coordinates[0][2], tolerance);
}

// tests collinearity of points with worst point near the centroid
TEST(NelderMead_FindMinDegenAssignAlmostCollinearBWG)
{
  Unfit::NelderMead object;
  std::vector<std::vector<double>> coordinates = {{3, 12}, {4, 13},
     {3.5, 12.50001}};
  double tolerance = 1e-4;
  SampleCostFunction3 cost_func;
  for (auto &vertex : coordinates) {
    auto residuals = cost_func(vertex);
    vertex.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  object.SetPopulation(coordinates);
  std::vector<double> coordinate {3, 12};
  int rc = object.FindMin(cost_func, coordinate);
  CHECK_EQUAL(0, rc);
  coordinates = object.GetPopulation();
  CHECK_CLOSE(0, coordinates[0][0], tolerance);
  CHECK_CLOSE(0, coordinates[0][1], tolerance);
  CHECK_CLOSE(0, coordinates[0][2], tolerance);
}


//
// ****************************************************************************
// Unit tests for FindMin (2 versions): 100% line, function & branch coverage
// ****************************************************************************
//
TEST(NelderMead_FindMinGivenVertex)
{
  Unfit::NelderMead object;
  SampleCostFunction3 cost_func;

  std::vector<double> vertex = {-1, -1};
  int rc = object.FindMin(cost_func, vertex);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0, vertex[0], 1e-4);
  CHECK_CLOSE(0, vertex[1], 1e-4);
}

TEST(NelderMead_FindMinGivenVertexWithRestart)
{
  Unfit::NelderMead object;
  SampleCostFunction3 cost_func;

  std::vector<double> vertex = {-1, -1};
  int rc = object.FindMin(cost_func, vertex);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0, vertex[0], 1e-4);
  CHECK_CLOSE(0, vertex[1], 1e-4);

  vertex = {-1, -1};
  rc = object.FindMin(cost_func, vertex);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0, vertex[0], 1e-4);
  CHECK_CLOSE(0, vertex[1], 1e-4);
}

TEST(NelderMead_FindMinGivenVertexWithAdaptiveParameters)
{
  Unfit::NelderMead object;
  object.options.SetUseAdaptiveParameters(true);
  SampleCostFunction3 cost_func;

  std::vector<double> vertex = {-1, -1};
  int rc = object.FindMin(cost_func, vertex);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0, vertex[0], 1e-4);
  CHECK_CLOSE(0, vertex[1], 1e-4);
  CHECK_CLOSE(cost_func(vertex)[0], sqrt(object.GetCost(0)), 1e-4);
}

TEST(NelderMead_FindMinGivenVertexWithUnworkableCostFunction)
{
  Unfit::NelderMead object;
  object.options.SetUseAdaptiveParameters(true);
  AlwaysInfinite cost_func;

  std::vector<double> vertex = {-1, -1};
  int rc = object.FindMin(cost_func, vertex);
  CHECK_EQUAL(-1, rc);
}

TEST(NelderMead_FindMinGivenVertexWithZeroLengthVertex)
{
  Unfit::NelderMead object;
  object.options.SetUseAdaptiveParameters(true);
  AlwaysInfinite cost_func;

  std::vector<double> vertex;
  int rc = object.FindMin(cost_func, vertex);
  CHECK_EQUAL(-4, rc);
}

TEST(NelderMead_FindMinGivenSimplex)
{
  Unfit::NelderMead object;
  SampleCostFunction3 cost_func;

  std::vector<std::vector<double>> simplex = {{-1, -1}, {1, -1}, {-1, 1}};
  for (auto &vertex : simplex) {
    auto residuals = cost_func(vertex);
    vertex.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  object.SetPopulation(simplex);
  std::vector<double> vertex {-1, -1};
  int rc = object.FindMin(cost_func, vertex);
  CHECK_EQUAL(0, rc);
  simplex = object.GetPopulation();
  CHECK_CLOSE(0, simplex[0][0], 1e-4);
  CHECK_CLOSE(0, simplex[0][1], 1e-4);
}

TEST(NelderMead_FindMinGivenSimplexWithRestart)
{
  Unfit::NelderMead object;
  SampleCostFunction3 cost_func;

  std::vector<std::vector<double>> simplex = {{-1, -1}, {1, -1}, {-1, 1}};
  for (auto &vertex : simplex) {
    auto residuals = cost_func(vertex);
    vertex.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  object.SetPopulation(simplex);
  std::vector<double> vertex1 {-1, -1};
  int rc = object.FindMin(cost_func, vertex1);
  CHECK_EQUAL(0, rc);
  simplex = object.GetPopulation();
  CHECK_CLOSE(0, simplex[0][0], 1e-4);
  CHECK_CLOSE(0, simplex[0][1], 1e-4);

  simplex = {{-1, -1}, {1, -1}, {-1, 1}};
  for (auto &vertex : simplex) {
    auto residuals = cost_func(vertex);
    vertex.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  object.SetPopulation(simplex);
  rc = object.FindMin(cost_func, vertex1);
  CHECK_EQUAL(0, rc);
  simplex = object.GetPopulation();
  CHECK_CLOSE(0, simplex[0][0], 1e-4);
  CHECK_CLOSE(0, simplex[0][1], 1e-4);
}

TEST(NelderMead_FindMinGivenSimplexWithAdaptiveParameters)
{
  Unfit::NelderMead object;
  object.options.SetUseAdaptiveParameters(true);
  SampleCostFunction3 cost_func;

  std::vector<std::vector<double>> simplex = {{-1, -1}, {1, -1}, {-1, 1}};
  for (auto &vertex : simplex) {
    auto residuals = cost_func(vertex);
    vertex.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  object.SetPopulation(simplex);
  std::vector<double> vertex {-1, -1};
  int rc = object.FindMin(cost_func, vertex);
  CHECK_EQUAL(0, rc);
  simplex = object.GetPopulation();
  CHECK_CLOSE(0, simplex[0][0], 1e-4);
  CHECK_CLOSE(0, simplex[0][1], 1e-4);
}

TEST(NelderMead_FindMinGivenSimplexOfZeroSize)
{
  Unfit::NelderMead object;
  SampleCostFunction3 cost_func;

  std::vector<std::vector<double>> simplex;
  object.SetPopulation(simplex);
  std::vector<double> vertex {-1, -1};
  int rc = object.FindMin(cost_func, vertex);
  CHECK_EQUAL(-1, rc);
}

TEST(NelderMead_FindMinGivenSimplexWithZeroLengthCoordinates)
{
  Unfit::NelderMead object;
  SampleCostFunction3 cost_func;

  std::vector<std::vector<double>> simplex = {{-1, -1}, {1, -1}, {-1, 1}};
  for (auto &vertex : simplex) {
    auto residuals = cost_func(vertex);
    vertex.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  simplex[0].clear();
  object.SetPopulation(simplex);
  std::vector<double> vertex {-1, -1};
  int rc = object.FindMin(cost_func, vertex);
  CHECK_EQUAL(-2, rc);
}

TEST(NelderMead_FindMinGivenSimplexWithWrongNumberOfVertices)
{
  Unfit::NelderMead object;
  SampleCostFunction3 cost_func;

  std::vector<std::vector<double>> simplex = {{-1, -1}, {1, -1}};
  for (auto &vertex : simplex) {
    auto residuals = cost_func(vertex);
    vertex.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  object.SetPopulation(simplex);
  std::vector<double> vertex {-1, -1};
  int rc = object.FindMin(cost_func, vertex);
  CHECK_EQUAL(-3, rc);
}

TEST(NelderMead_FindMinGivenSimplexWithInconsistentCoordinates)
{
  Unfit::NelderMead object;
  SampleCostFunction3 cost_func;

  std::vector<std::vector<double>> simplex = {{-1, -1}, {1}, {-1, 1}};
  for (auto &vertex : simplex) vertex.emplace_back(1.0);  // Dummy cost
  object.SetPopulation(simplex);
  std::vector<double> vertex {-1, -1};
  int rc = object.FindMin(cost_func, vertex);
  CHECK_EQUAL(-4, rc);
}

TEST(NelderMead_FindMinGivenSimplexWithUnworkableCostFunction)
{
  Unfit::NelderMead object;
  AlwaysInfinite cost_func;

  std::vector<std::vector<double>> simplex = {{-1, -1}, {1, -1}, {-1, 1}};
  for (auto &vertex : simplex) {
    auto residuals = cost_func(vertex);
    vertex.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  object.SetPopulation(simplex);
  std::vector<double> vertex {-1, -1};
  int rc = object.FindMin(cost_func, vertex);
  CHECK_EQUAL(-5, rc);
}
//
// ****************************************************************************
// End of unit tests for FindMin (2 versions)
// ****************************************************************************
//

TEST_FIXTURE(TestNelderMead, toUseAdaptiveParameters)
{
  // Initialise the initial guess and the object for the cost function
  std::vector<std::vector<double>> v1 = {{2, 4, 4},
    {56, 234234, -4}, {45, 12412, 2}};
  SampleCostFunction1 Func;
  // Change the flag so that findmin uses adaptive parameters which is
  // documented in the pdf
  options.SetUseAdaptiveParameters(true);

  // Minimise so that the parameters would be changed
  FindMin(Func, v1[0]);

  // Check whether adaptive values have been used
  double alpha, beta, delta, gamma;
  options.GetNelderMeadStepSizes(alpha, beta, delta, gamma);
  CHECK_CLOSE(1, alpha, 1e-3);
  CHECK_CLOSE(1.666667, beta, 1e-3);
  CHECK_CLOSE(0.666667, delta, 1e-3);
  CHECK_CLOSE(0.583333, gamma, 1e-3);
}

TEST(NelderMead_GetSolutionGetCostIndexOutOfBounds)
{
  Unfit::NelderMead object;
  SampleCostFunction3 cost_func;

  std::vector<double> vertex = {-1, -1};
  int rc = object.FindMin(cost_func, vertex);
  CHECK_EQUAL(0, rc);
  // Index out of bounds returns an infinite cost
  auto cost = object.GetCost(10);
  CHECK(!std::isfinite(cost));
  // Index out of bounds returns an empty vertex
  auto solution = object.GetSolution(10);
  CHECK(solution.empty());
}

//// Test Expand with a different beta value
// TEST_FIXTURE(TestNelderMead, MultipleExpandwithdifferentbeta)
// {
//  // Set Beta
//  double alpha, beta, delta, gamma;
//  options.GetNelderMeadStepSizes(alpha, beta, delta, gamma);
//  options.SetNelderMeadStepSizes(alpha, 3.0, delta, gamma);
//
//  // Initialise the population_ with a generated simplex
//  std::vector<std::vector<double>> v1 = {{1, 2, 3, 4, 1}, {5, 6, 7, 8, 4},
//      {9, 10, 11, 12, 2}, {13, 14, 15, 16, 2}, {17, 18, 19, 20, 4}};
//
//  // Set number of dimensions and initialise all the working vectors
//  SetNumberOfDimensions_(4);
//  AccessInitialiseVectors();
//  SetVertices(v1);
//
//  // Compute and set reflect and centroid values
//  std::vector<double> v2 = {1, 2, 3, 4, 5};
//  SetReflect(v2);
//  v2 = {11.25, 12.5, 13.75, 15, 0};
//  SetCentroid(v2);
//
//  // Compute the expanded point using the expand function
//  PowellSingular func;
//  CHECK(AccessExpand(func));
//
//  // Initialise an empty vector, get the expanded points and test the values
//  GetExpand(v2);
//  CHECK_CLOSE(-19.5 , v2[0], 1e-16);
//  CHECK_CLOSE(-19, v2[1], 1e-16);
//  CHECK_CLOSE(-18.5, v2[2], 1e-16);
//  CHECK_CLOSE(-18, v2[3], 1e-16);
// }
}  // suite UnitTestNelderMead
}  // namespace UnitTests
}  // namespace Unfit
