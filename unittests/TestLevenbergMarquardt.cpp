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
#include <vector>
#include "DataFileReader.hpp"
#include "GaussianEquation.hpp"
#include "GenericCostFunction.hpp"
#include "Exponential.hpp"
#include "Matrix.hpp"
#include "LevenbergMarquardt.hpp"
#include "LevenbergMarquardtTestFunctions.hpp"
#include "Parabolic.hpp"
#include "UnitTest++.h"
#include "Woods.hpp"

namespace Unfit
{
/**
 * This class is designed solely to provide the functionality to access the
 * private methods and member variables of the LevenbergMarquardt class.
 */
class TestLevenbergMarquardt : public LevenbergMarquardt
{
 public:
  TestLevenbergMarquardt() : LevenbergMarquardt() {}

  /**
   * A function to access the private method (FindJacobian) in the
   * LevenbergMarquardt class.
   *
   * \param cost_function sample functions to be tested
   * \param residuals a vector of residuals at the current guess. Only used for
   *        one-sided differencing; can be empty for two-sided differencing.
   * \param variables vector of the coefficients of the function
   * \param one_sided_difference select between one sided (forward) and two
   *        sided finite different approximations for the Jacobian calculation
   */
  void AccessFindJacobian(GenericCostFunction &cost_function,
      const std::vector<double> &residuals,
      const std::vector<double> &variables, bool one_sided_difference = true)
  {
    FindJacobian(cost_function, residuals, variables, one_sided_difference);
  }

  /**
   * A function to access the private method (BroydenUpdate) in the
   * LevenbergMarquardt class.
   *
   * \param residuals_new residual vector from the candidate solution
   * \param residuals residual vector from the current solution
   * \param hlm vector containing the step between the current and new solutions
   * \param jacobianT transpose of the Jacobian matrix
   */
  void AccessBroydenUpdate(const std::vector<double> &residuals_new,
    const std::vector<double> &residuals, const std::vector<double> &hlm,
    Matrix &jacobianT)
  {
    BroydenUpdate(residuals_new, residuals, hlm, jacobianT);
  }

  /**
   * A function to access the private member variable (jacobian_) in the
   * LevenbergMarquardt class.
   *
   * \return The jacobian matrix (transposed)
   */
  Matrix AccessJacobian()
  {
    return jacobian_;
  }

  /**
   * A function to access the private method (GetSolution) in the
   * LevenbergMarquardt class.
   *
   * \return An empty vector, as the solution is not stored internally
   */
  std::vector<double> AccessGetSolution()
  {
    return GetSolution();
  }

  /**
   * A function to set the private member variable (dimensions_) in the
   * LevenbergMarquardt class. This is the number of variables to be fitted.
   *
   * \param dimensions the number of variables in the problem
   */
  void SetDimensions(std::size_t dimensions)
  {
    dimensions_ = dimensions;
  }

  /**
   * A function to get the private member variable (dimensions_) in the
   * LevenbergMarquardt class. This is the number of variables to be fitted.
   *
   * \return the number of variables in the problem
   */
  std::size_t GetDimensions() const
  {
    return dimensions_;
  }

  /**
   * A function to set the private member variable (observation_size_) in the
   * LevenbergMarquardt class. This is the number of data points in the fit.
   *
   * \param obs_size the number of data points to be fitted
   */
  void SetObservationSize(std::size_t obs_size)
  {
    observation_size_ = obs_size;
  }

  /**
   * A function to get the private member variable (observation_size_) in the
   * LevenbergMarquardt class. This is the number of data points in the fit.
   *
   * \return the number of data points to be fitted
   */
  std::size_t GetObservationSize() const
  {
    return observation_size_;
  }
};

namespace UnitTests
{
SUITE(UnitTestLevenbergMarquardt)
{
TEST_FIXTURE(TestLevenbergMarquardt, FindJacobian_NullVariables)
{
  std::vector<double> x {1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> y {6.0, 11.0, 18.0, 27.0, 38.0};
  Unfit::Examples::Parabolic parabolic_func(x, y);
  std::vector<double> variables = {};
  std::vector<double> residuals = {};
  AccessFindJacobian(parabolic_func, residuals, variables, false);
}

TEST_FIXTURE(TestLevenbergMarquardt, FindJacobian_Parabolic)
{
  SetObservationSize(5);
  SetDimensions(3);
  std::vector<double> x {1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> y {6.0, 11.0, 18.0, 27.0, 38.0};
  Unfit::Examples::Parabolic parabolic_func(x, y);
  std::vector<double> variables = {1.0, 1.0, 1.0};
  std::vector<double> residuals(y);

  AccessFindJacobian(parabolic_func, residuals, variables, false);
  Matrix expected {3, 5, { -1, -4, -9, -16, -25,
                           -1, -2, -3, -4, -5,
                           -1, -1, -1, -1, -1}};

  auto num_col = x.size();
  auto num_row = variables.size();

  CHECK_EQUAL(15u, AccessJacobian().Size());
  for (auto i = 0u; i < num_row*num_col; ++i) {
     CHECK_CLOSE(expected.values_[i], AccessJacobian().values_[i], 1e-4);
  }
}

TEST_FIXTURE(TestLevenbergMarquardt, FindJacobian_Exponential)
{
  SetObservationSize(4);
  SetDimensions(3);
  std::vector<double> variables {0.0, 0.35, -1.0};
  std::vector<double> x {-1.0, -0.5, 0.0, 0.5};
  std::vector<double> y {-0.259240, -0.308818, -0.367879, -0.438234};
  std::vector<double> residuals(y);
  Unfit::Examples::Exponential exponential_func(x, y);
  AccessFindJacobian(exponential_func, residuals, variables, false);
  Matrix expected {3, 4, {-1, -1, -1, -1,
                          -0.259240, -0.154409, 0, 0.219117,
                           0.259240,  0.308818, 0.367879, 0.438234}};
  auto num_col = x.size();
  auto num_row = variables.size();
  for (auto i = 0u; i < num_row*num_col; ++i) {
     CHECK_CLOSE(expected.values_[i], AccessJacobian().values_[i], 1e-4);
  }
}

TEST_FIXTURE(TestLevenbergMarquardt, FindJacobian_Exponential2)
{
  SetObservationSize(4);
  SetDimensions(3);
  std::vector<double> variables {2.0, -1.2, 0.5};
  std::vector<double> x {-1.0, -0.5, 0.0, 0.5};
  std::vector<double> y {7.4739, 5.0041, 3.6487, 2.9048};
  std::vector<double> residuals(y);
  Unfit::Examples::Exponential exponential_func(x, y);
  AccessFindJacobian(exponential_func, residuals, variables, false);
  Matrix expected {3, 4, {-1, -1, -1, -1,
                      -5.4739473917, -1.502083012, 0, 0.452418709,
                      5.4739473917, 3.0041660239, 1.6487212707, 0.904837418}};
  auto num_col = x.size();
  auto num_row = variables.size();
  for (auto i = 0u; i < num_row*num_col; ++i) {
    CHECK_CLOSE(expected.values_[i], AccessJacobian().values_[i], 1e-4);
  }
}

TEST_FIXTURE(TestLevenbergMarquardt, FindJacobian_Gaussian)
{
  SetObservationSize(6);
  SetDimensions(4);
  options.SetCostTolerance(1e-8);
  std::vector<double> variables {1.0, 2.0, 3.0, 0.0};
  std::vector<double> x {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> y {0.800737, 0.945959, 1.0, 0.945959, 0.80073, 0.60653};
  std::vector<double> residuals(y);
  Unfit::Examples::GaussianEquation gaussian_func(x, y);
  AccessFindJacobian(gaussian_func, residuals, variables, false);
  Matrix expected {6, 4, {-0.80073740291,  0.17794164500,  -0.11862776357, -1,
                          -0.94595946890,  0.10510660759, -0.035035535957, -1,
                          -0.99999999999,              0,               0, -1,
                          -0.94595946890, -0.10510660760, -0.035035535957, -1,
                          -0.80073740291, -0.17794164500,  -0.11862776357, -1,
                          -0.60653065971, -0.20217688649,  -0.20217688672, -1}};
  Matrix expected_trans = Transpose(expected);
  auto num_col = x.size();
  auto num_row = variables.size();
  for (auto i = 0u; i < num_row*num_col; ++i) {
    CHECK_CLOSE(expected_trans.values_[i], AccessJacobian().values_[i], 1e-6);
  }
}

TEST_FIXTURE(TestLevenbergMarquardt, FindJacobian_Gaussian1)
{
  SetObservationSize(1);
  SetDimensions(4);
  std::vector<double> variables {2.0, 2.0, 3.0, 0.0};
  std::vector<double> x {1.0};
  std::vector<double> y {2.0};
  std::vector<double> residuals(y);
  Unfit::Examples::GaussianEquation gaussian_func(x, y);
  AccessFindJacobian(gaussian_func, residuals, variables, false);
  Matrix expected {1, 4, {-0.94595, 0.2102132, -0.07007107, -1}};

  auto num_col = x.size();
  auto num_row = variables.size();
  for (auto i = 0u; i < num_row*num_col; ++i) {
    CHECK_CLOSE(expected.values_[i], AccessJacobian().values_[i], 1e-3);
  }
}

TEST_FIXTURE(TestLevenbergMarquardt, FindJacobian_Gaussian2)
{
  SetObservationSize(2);
  SetDimensions(4);
  std::vector<double> variables {2.0, 2.0, 3.0, 0.0};
  std::vector<double> x {1.0, 2.0};
  std::vector<double> y {2.0, 1.60653};
  std::vector<double> residuals(y);
  Unfit::Examples::GaussianEquation gaussian_func(x, y);
  AccessFindJacobian(gaussian_func, residuals, variables, false);
  Matrix expected {2, 4, {-0.94595, 0.21021, -0.07007, -1,
                          -1,            0,        0,  -1}};
  Matrix expected_trans = Transpose(expected);
  auto num_col = x.size();
  auto num_row = variables.size();
  for (auto i = 0u; i < num_row*num_col; ++i) {
    CHECK_CLOSE(expected_trans.values_[i], AccessJacobian().values_[i], 1e-3);
  }
}

TEST_FIXTURE(TestLevenbergMarquardt, FindJacobian_Gaussian3)
{
  SetObservationSize(4);
  SetDimensions(4);
  std::vector<double> variables {2.0, 2.0, 3.0, 0.0};
  std::vector<double> x {1.0, 2.0, 3.0, 4.0};
  std::vector<double> y {1.89191, 2.0, 1.89191, 1.60147};
  std::vector<double> residuals(y);
  Unfit::Examples::GaussianEquation gaussian_func(x, y);
  AccessFindJacobian(gaussian_func, residuals, variables, false);
  Matrix expected {4, 4, {-0.94595, 0.21021,  -0.07007, -1,
                                -1,       0,         0, -1,
                          -0.94595, -0.21021, -0.07007, -1,
                          -0.80074, -0.35588, -0.23726, -1}};
  Matrix expected_trans = Transpose(expected);
  auto num_col = x.size();
  auto num_row = variables.size();
  for (auto i = 0u; i < num_row*num_col; ++i) {
    CHECK_CLOSE(expected_trans.values_[i], AccessJacobian().values_[i], 1e-3);
  }
}

TEST(LevenbergMarquardt_CheckMaxIterations)
{
  LevenbergMarquardt test;
  test.options.SetMaxIterations(123);
  CHECK_EQUAL(123u, test.options.GetMaxIterations());
}

TEST(LevenbergMarquardt_CheckMaxFunctionEvaluation)
{
  LevenbergMarquardt test;
  test.options.SetMaxFunctionEvaluations(124);
  CHECK_EQUAL(124u, test.options.GetMaxFunctionEvaluations());
}

TEST_FIXTURE(TestLevenbergMarquardt, CheckDimension)
{
  SetDimensions(4);
  CHECK_EQUAL(4u, GetDimensions());
}

TEST_FIXTURE(TestLevenbergMarquardt, CheckObservationSize)
{
  SetObservationSize(4);
  CHECK_EQUAL(4u, GetObservationSize());
}

TEST_FIXTURE(TestLevenbergMarquardt, CheckObservationSize2)
{
  LevenbergMarquardt test;
  SetObservationSize(123);
  CHECK_EQUAL(123u, GetObservationSize());
}

TEST(LevenbergMarquardt_CheckOutputLevel)
{
  LevenbergMarquardt test;
  test.options.SetOutputLevel(1);
  CHECK_EQUAL(1u, test.options.GetOutputLevel());
}

TEST_FIXTURE(TestLevenbergMarquardt, GetSolution)
{
  const auto solution = AccessGetSolution();
  CHECK(solution.empty());
}

TEST(LevenbergMarquardt_CheckMaxFunctionEvaluations)
{
  LevenbergMarquardt test;
  test.options.SetMaxFunctionEvaluations(1);
  CHECK_EQUAL(1u, test.options.GetMaxFunctionEvaluations());
}

TEST(LevenbergMarquardt_InitialGuessIsCorrect)
{
  LevenbergMarquardt test;
  Unfit::Examples::Woods cost_func;
  // Initial guess (correct answer)
  std::vector<double> min_point = {1.0, 1.0, 1.0, 1.0};
  // Minimise
  auto rc = test.FindMin(cost_func, min_point);
  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1, min_point[0], 1e-3);
  CHECK_CLOSE(1, min_point[1], 1e-3);
  CHECK_CLOSE(1, min_point[2], 1e-3);
  CHECK_CLOSE(1, min_point[3], 1e-3);

  // Repeat with output on
  test.options.SetOutputLevel(3);
  min_point = {1.0, 1.0, 1.0, 1.0};
  rc = test.FindMin(cost_func, min_point);
  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1, min_point[0], 1e-3);
  CHECK_CLOSE(1, min_point[1], 1e-3);
  CHECK_CLOSE(1, min_point[2], 1e-3);
  CHECK_CLOSE(1, min_point[3], 1e-3);
}

TEST(LevenbergMarquardt_DifferentOutputLevels)
{
  Unfit::Examples::Woods cost_func;

  // No output
  LevenbergMarquardt test1;
  std::vector<double> min_point = {1.01, 1.0, 1.0, 1.0};
  auto rc = test1.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1, min_point[0], 1e-3);
  CHECK_CLOSE(1, min_point[1], 1e-3);
  CHECK_CLOSE(1, min_point[2], 1e-3);
  CHECK_CLOSE(1, min_point[3], 1e-3);

  // Set output to level 1 (iteration count)
  LevenbergMarquardt test2;
  test2.options.SetOutputLevel(1);
  min_point = {1.01, 1.0, 1.0, 1.0};
  rc = test2.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1, min_point[0], 1e-3);
  CHECK_CLOSE(1, min_point[1], 1e-3);
  CHECK_CLOSE(1, min_point[2], 1e-3);
  CHECK_CLOSE(1, min_point[3], 1e-3);

  // Set output to level 2 (iteration count + costs)
  LevenbergMarquardt test3;
  test3.options.SetOutputLevel(2);
  min_point = {1.01, 1.0, 1.0, 1.0};
  rc = test3.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1, min_point[0], 1e-3);
  CHECK_CLOSE(1, min_point[1], 1e-3);
  CHECK_CLOSE(1, min_point[2], 1e-3);
  CHECK_CLOSE(1, min_point[3], 1e-3);

  // Set output to level 3 (iteration count + costs + final result)
  LevenbergMarquardt test4;
  test4.options.SetOutputLevel(3);
  min_point = {1.01, 1.0, 1.0, 1.0};
  rc = test4.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1, min_point[0], 1e-3);
  CHECK_CLOSE(1, min_point[1], 1e-3);
  CHECK_CLOSE(1, min_point[2], 1e-3);
  CHECK_CLOSE(1, min_point[3], 1e-3);
}

// TEST_FIXTURE(TestLevenbergMarquardt, Constructor1)
// {
//  Matrix matrix {3, 3, {0, 1, 2,
//                        3, 4, 5,
//                        6, 7, 8}};
//  std::vector<double> jacobian = AccessJacobian();
//  jacobian = {0.0, 1.0, 2.0};
//  for (unsigned i = 0; i < 3; ++i) {
//    CHECK_EQUAL(i, jacobian[i]);
//  }
//  for (unsigned i = 0; i < 3; ++i) {
//    for (unsigned j = 0; j < 3; ++j) {
//      CHECK_EQUAL(i*3u+j, matrix.GetValue(i, j));
//    }
//  }
// }

TEST(LevenbergMarquardt_MaxFunctionEvaluation)
{
  Unfit::LevenbergMarquardt test;
  test.options.SetMaxFunctionEvaluations(1);
  std::vector<double> x;
  std::vector<double> y;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/parabolic_data.txt"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, x));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, y));
  Unfit::Examples::Parabolic parabolic_func(x, y);
  std::vector<double> init_guess {1.0, 1.0, 1.0};
  CHECK_EQUAL(4, test.FindMin(parabolic_func, init_guess));
}

TEST(LevenbergMarquardt_MaxIterations)
{
  Unfit::LevenbergMarquardt test;
  test.options.SetMaxIterations(1);
  std::vector<double> x {-5.0, -4.0, -3.0, -2.0, -1.0, 0.0,
                         1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> y {18.0, 11.0, 6.0, 3.0, 2.0, 3.0,
                         6.0, 11.0, 18.0, 27.0, 38.0};
  Unfit::Examples::Parabolic para_eq(x, y);
  std::vector<double> init_guess {1.0, 1.0, 1.0};
  CHECK_EQUAL(5, test.FindMin(para_eq, init_guess));
}

TEST(LevenbergMarquardt_InvalidInitialGuessEmpty)
{
  Unfit::LevenbergMarquardt test;
  test.options.SetMaxIterations(1);
  std::vector<double> x {-5.0, -4.0, -3.0, -2.0, -1.0, 0.0,
                         1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> y {18.0, 11.0, 6.0, 3.0, 2.0, 3.0,
                         6.0, 11.0, 18.0, 27.0, 38.0};
  Unfit::Examples::Parabolic para_eq(x, y);
  std::vector<double> init_guess {};
  CHECK_EQUAL(1, test.FindMin(para_eq, init_guess));
}

TEST(LevenbergMarquardt_InvalidInitialGuessExceedsDataLength)
{
  Unfit::LevenbergMarquardt test;
  std::vector<double> x {-5.0, -4.0, -3.0, -2.0, -1.0, 0.0};
  std::vector<double> y {18.0, 11.0, 6.0, 3.0, 2.0, 3.0};
  Unfit::Examples::Parabolic para_eq(x, y);
  std::vector<double> init_guess {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  CHECK_EQUAL(3, test.FindMin(para_eq, init_guess));
}

TEST(LevenbergMarquardt_InvalidInitialGuessOutOfBounds)
{
  Unfit::LevenbergMarquardt test;
  test.bounds.SetBounds({0.0, 0.0, 0.0}, {0.9, 1.9, 2.9});
  std::vector<double> x {-5.0, -4.0, -3.0, -2.0, -1.0, 0.0,
                         1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> y {18.0, 11.0, 6.0, 3.0, 2.0, 3.0,
                         6.0, 11.0, 18.0, 27.0, 38.0};
  Unfit::Examples::Parabolic para_eq(x, y);
  std::vector<double> init_guess {1.0, 1.0, 1.0};
  CHECK_EQUAL(2, test.FindMin(para_eq, init_guess));
}

TEST(LevenbergMarquardt_FindMinParabolicHitMaxIterations)
{
  Unfit::LevenbergMarquardt test;
  test.options.SetMaxIterations(1);
  std::vector<double> x {-5.0, -4.0, -3.0, -2.0, -1.0, 0.0,
                         1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> y {18.0, 11.0, 6.0, 3.0, 2.0, 3.0,
                         6.0, 11.0, 18.0, 27.0, 38.0};
  Unfit::Examples::Parabolic para_eq(x, y);
  std::vector<double> init_guess {0, 0, 0};
  int rc = test.FindMin(para_eq, init_guess);
  CHECK_EQUAL(5, rc);
}

TEST(LevenbergMarquardt_SampleEquation1)
{
  Unfit::LevenbergMarquardt object;
  Unfit::Examples::Equation1 cost_func;
  // Initial guess
  std::vector<double> min_point = {10.0};
  // Check the function calculates the correct cost at the initial guess
  std::vector<double> residual = cost_func(min_point);
  double sum {0.0};
  for (unsigned i = 0; i < residual.size(); ++i) sum += residual[i]*residual[i];
  CHECK_CLOSE(14641, sum, 1e-4);
  // Minimise
  auto rc = object.FindMin(cost_func, min_point);
  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(-1, min_point[0], 1e-3);
}

TEST(LevenbergMarquardt_SampleEquation2)
{
  Unfit::LevenbergMarquardt object;
  Unfit::Examples::Equation2 cost_func;
  object.options.SetCostTolerance(1e-21);
  object.options.SetEpsilon(1e-21);
  // Initial guess
  std::vector<double> min_point = {-2.0, 2.0};
  // Check the function calculates the correct cost at the initial guess
  std::vector<double> residual = cost_func(min_point);
  double sum {0.0};
  for (unsigned i = 0; i < residual.size(); ++i) sum += residual[i]*residual[i];
  CHECK_CLOSE(1536, sum, 1e-4);
  // Minimise
  auto rc = object.FindMin(cost_func, min_point);
  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0, min_point[0], 1e-3);
  CHECK_CLOSE(0, min_point[1], 1e-3);
}

TEST(LevenbergMarquardt_SampleEquation2NoBroydenUpdates)
{
  Unfit::LevenbergMarquardt object;
  Unfit::Examples::Equation2 cost_func;
  object.options.SetCostTolerance(1e-21);
  object.options.SetEpsilon(1e-21);
  object.options.SetUseBroydenUpdates(false);
  // Initial guess
  std::vector<double> min_point = {-2.0, 2.0};
  // Check the function calculates the correct cost at the initial guess
  std::vector<double> residual = cost_func(min_point);
  double sum {0.0};
  for (unsigned i = 0; i < residual.size(); ++i) sum += residual[i]*residual[i];
  CHECK_CLOSE(1536, sum, 1e-4);
  // Minimise
  auto rc = object.FindMin(cost_func, min_point);
  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0, min_point[0], 1e-3);
  CHECK_CLOSE(0, min_point[1], 1e-3);
}

TEST(LevenbergMarquardt_SampleEquation3)
{
  Unfit::LevenbergMarquardt object;
  Unfit::Examples::Equation3 cost_func;
  // Initial guess
  std::vector<double> min_point = {-2.0, 2.0};
  // Check the function calculates the correct cost at the initial guess
  std::vector<double> residual = cost_func(min_point);
  double sum {0.0};
  for (unsigned i = 0; i < residual.size(); ++i) sum += residual[i]*residual[i];
  CHECK_CLOSE(82, sum, 1e-4);
  // Minimise
  auto rc = object.FindMin(cost_func, min_point);
  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1.0, min_point[0], 1e-3);
  CHECK_CLOSE(1.0, min_point[1], 1e-3);
}

TEST(LevenbergMarquardt_SampleEquation4)
{
  Unfit::LevenbergMarquardt object;
  Unfit::Examples::Equation4 cost_func;
  object.options.SetCostTolerance(1e-21);
  object.options.SetEpsilon(1e-21);
  // Initial guess
  std::vector<double> min_point = {-2.0};
  // Check the function calculates the correct cost at the initial guess
  std::vector<double> residual = cost_func(min_point);
  double sum {0.0};
  for (unsigned i = 0; i < residual.size(); ++i) sum += residual[i]*residual[i];
  CHECK_CLOSE(4096, sum, 1e-4);
  // Minimise
  auto rc = object.FindMin(cost_func, min_point);
  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 2e-3);
}

TEST(LevenbergMarquardt_GetCost)
{
  std::vector<double> min_point {5.0, 1.0, -1.0};
  std::vector<double> x {-1.0, -0.5, 0.0, 0.5};
  std::vector<double> y {-0.259240, -0.308818, -0.367879, -0.438234};
  Unfit::Examples::Exponential cost_func(x, y);
  Unfit::LevenbergMarquardt lm;
  lm.options.SetMaxIterations(1);
  auto rc = lm.FindMin(cost_func, min_point);
  CHECK_EQUAL(5, rc);
  CHECK_CLOSE(0.00369037, lm.GetCost(), 1.0e-6);
}

TEST(LevenbergMarquardt_Reset)
{
  std::vector<double> min_point {5.0, 1.0, -1.0};
  std::vector<double> x {-1.0, -0.5, 0.0, 0.5};
  std::vector<double> y {-0.259240, -0.308818, -0.367879, -0.438234};
  Unfit::Examples::Exponential cost_func(x, y);
  Unfit::LevenbergMarquardt lm;
  CHECK_CLOSE(1e-10, lm.options.GetDelta(), 1.0e-12);
  CHECK_CLOSE(0.001, lm.options.GetTau(), 1.0e-5);
  CHECK_CLOSE(1e-12, lm.options.GetCostTolerance(), 1.0e-14);
  lm.options.SetDelta(1.0e-6);
  lm.options.SetTau(0.1);
  lm.options.SetCostTolerance(1e-15);
  CHECK_CLOSE(1e-6, lm.options.GetDelta(), 1.0e-8);
  CHECK_CLOSE(0.1, lm.options.GetTau(), 1.0e-3);
  CHECK_CLOSE(1e-15, lm.options.GetCostTolerance(), 1.0e-17);
  auto rc = lm.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(2.81916e-13, lm.GetCost(), 1.0e-16);

  lm.Reset();
  CHECK_CLOSE(1e-10, lm.options.GetDelta(), 1.0e-12);
  CHECK_CLOSE(0.001, lm.options.GetTau(), 1.0e-5);
  CHECK_CLOSE(1e-12, lm.options.GetCostTolerance(), 1.0e-14);
  rc = lm.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(2.81916e-13, lm.GetCost(), 1.0e-16);
}

TEST_FIXTURE(TestLevenbergMarquardt, BroydenUpdate)
{
  std::vector<double> jacobianT_data {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  Matrix jacobianT(2, 3, jacobianT_data);
  std::vector<double> hlm {-0.2, -0.1};
  std::vector<double> residuals_new {0.3, 0.4, 0.5};
  std::vector<double> residuals {0.1, 0.1, 0.2};

  // JT = [1, 2, 3]
  //      [4, 5, 6]
  // J = [1, 4]
  //     [2, 5]
  //     [3, 6]
  // J*h = [-0.6, -0.9, -1.2]
  // (rnew - r - Jh) = [0.8, 1.2, 1.5]
  // hTh = 0.05
  // u = (1/hth)*(rnew - r - Jh) = [16, 24, 30]
  // u*hT = [-3.2, -1.6]
  //        [-4.8, -2.4]
  //        [-6.0, -3.0]
  // (u*hT)T  = [-3.2, -4.8, -6.0]
  //            [-1.6, -2.4, -3.0]
  // JT(new) = JT + (u*hT)T
  //         = [-2.2, -2.8, -3.0]
  //           [ 2.4,  2.6,  3.0]

  AccessBroydenUpdate(residuals_new, residuals, hlm, jacobianT);
  CHECK_CLOSE(-2.2, jacobianT.values_[0], 1.0e-8);
  CHECK_CLOSE(-2.8, jacobianT.values_[1], 1.0e-8);
  CHECK_CLOSE(-3.0, jacobianT.values_[2], 1.0e-8);
  CHECK_CLOSE(2.4, jacobianT.values_[3], 1.0e-8);
  CHECK_CLOSE(2.6, jacobianT.values_[4], 1.0e-8);
  CHECK_CLOSE(3.0, jacobianT.values_[5], 1.0e-8);
}

TEST(LevenbergMarquardt_GetInvalidCostDuringIterations)
{
  Unfit::Examples::ParabolaWithHole cost_func;
  Unfit::LevenbergMarquardt lm;
  std::vector<double> min_point {1.0};
  auto rc = lm.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.1, min_point[0], 1.0e-6);
}
}  // suite UnitTestLevenbergMarquardt
}  // namespace UnitTests
}  // namespace Unfit
