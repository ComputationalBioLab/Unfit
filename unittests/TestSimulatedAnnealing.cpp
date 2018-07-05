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
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include "GenericCostFunction.hpp"
#include "SimulatedAnnealing.hpp"
#include "TestFunctions.hpp"
#include "UnitTest++.h"

namespace Unfit
{
/**
 * This class is designed solely to provide the functionality to access the
 * private methods and private members of the SimulatedAnnealing class.
 */
class TestSimulatedAnnealing : public SimulatedAnnealing
{
 public:
  /**
   * This method provides access to the private method InitialiseParameters
   * from the SimulatedAnnealing class for testing purposes.
   */
  void AccessInitialiseParameters()
  {
    InitialiseParameters();
  }

  /**
   * This method provides access to the private method GenerateTrialPoint
   * from the SimulatedAnnealing class for testing purposes.
   *
   * \param trial_point The coordinate vector to be perturbed
   * \param i The index (coordinate) to perturb
   */
  void AccessGenerateTrialPoint(std::vector<double> &trial_point, int i)
  {
    GenerateTrialPoint(trial_point, i);
  }

  /**
   * This method provides access to the private method UpdateStepSizes
   * from the SimulatedAnnealing class for testing purposes.
   */
  void AccessUpdateStepSizes()
  {
    UpdateStepSizes();
  }

  /**
   * This method provides access to the private method ResetStepSizes
   * from the SimulatedAnnealing class for testing purposes.
   *
   * \param step_size A value used for scaling the step size
   */
  void AccessResetStepSizes(double step_size)
  {
    ResetStepSizes(step_size);
  }

  /**
   * A function to access the private random number generator in the
   * SimulatedAnnealing class and change its seed.
   *
   * \param seed The new seed
   */
  void AccessSetGeneratorSeed(unsigned seed)
  {
    generator_.seed(seed);
  }

  /**
   * This method allows the number of dimensions of the problem to be set
   * directly for testing purposes. dimensions_ is a private
   * member of the SimulatedAnnealing Class.
   *
   * \param dimensions The number of dimensions required
   */
  void SetDimensions(std::size_t dimensions)
  {
    dimensions_ = dimensions;
  }

  /**
   * Get the private member dimensions_ of the SimulatedAnnealing class
   *
   * \return The number of dimensions
   */
  std::size_t GetDimensions()
  {
    return dimensions_;
  }

  /**
   * Get the private member cost_ of the SimulatedAnnealing class
   *
   * \return The cost index
   */
  std::size_t GetCostIndex()
  {
    return cost_;
  }

  /**
   * Get the private member previous_best_cost_ of the SimulatedAnnealing class
   *
   * \return The cost index
   */
  double GetPreviousBestCost()
  {
    return previous_best_cost_;
  }

  /**
   * Set the private member step_sizes_ of the SimulatedAnnealing class
   *
   * \param step_sizes The vector of step sizes
   */
  void SetStepSizes(std::vector<double> &step_sizes)
  {
    step_sizes_ = step_sizes;
  }

  /**
   * Get the private member step_sizes_ of the SimulatedAnnealing class
   *
   * \return The vector of step sizes
   */
  std::vector<double> GetStepSizes()
  {
    return step_sizes_;
  }

  /**
   * Set the private member acceptance_ratios_ of the SimulatedAnnealing class
   *
   * \param ratios The vector of acceptance ratios
   */
  void SetAcceptanceRatios(std::vector<double> &ratios)
  {
    acceptance_ratios_ = ratios;
  }

  /**
   * Get the private member acceptance_ratios_ of the SimulatedAnnealing class
   *
   * \return The vector of acceptance ratios
   */
  std::vector<double> GetAcceptanceRatios()
  {
    return acceptance_ratios_;
  }

  /**
   * A function to access the private random number generator in the
   * SimulatedAnnealing class, and the uniform distribution.
   *
   * \return A uniform random number on the interval [0, 1]
   */
  double GetRandomNumber()
  {
    return uniform_dist_(generator_);
  }
};

namespace UnitTests
{
SUITE(UnitTestSimulatedAnnealing)
{
TEST_FIXTURE(TestSimulatedAnnealing, TestConstructor)
{
  // Cost index
  CHECK_EQUAL(0u, GetCostIndex());
  // Number of dimensions
  CHECK_EQUAL(0u, GetDimensions());
  // Previous best cost
  CHECK_CLOSE(0.0, GetPreviousBestCost(), 1e-12);
  // Step size
  auto step_sizes = GetStepSizes();
  CHECK(step_sizes.empty());
  // Acceptance ratios
  auto acceptance_ratios = GetAcceptanceRatios();
  CHECK(acceptance_ratios.empty());
  // Random number generator
  for (auto i = 0; i < 10; ++i) {
    CHECK(GetRandomNumber() >= 0.0);
    CHECK(GetRandomNumber() <= 1.0);
  }
}

TEST_FIXTURE(TestSimulatedAnnealing, TestReset)
{
  // A 2D problem with default parameters
  SimpleCostFunction cost_func;
  bounds.SetBounds(0, -10.0, 10.0);
  bounds.SetBounds(1, -10.0, 10.0);
  std::vector<double> min_point {1.0, 1.0};
  auto rc = FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 0.01);
  CHECK_CLOSE(0.0, min_point[1], 0.01);

  // Check the state variables have been set correctly
  // Cost index
  CHECK_EQUAL(2u, GetCostIndex());
  // Number of dimensions
  CHECK_EQUAL(2u, GetDimensions());
  // Previous best cost
  CHECK_CLOSE(0.0, GetPreviousBestCost(), 1e-6);
  // Step size
  auto step_sizes = GetStepSizes();
  CHECK_EQUAL(2u, step_sizes.size());
  // Acceptance ratios
  auto acceptance_ratios = GetAcceptanceRatios();
  CHECK_EQUAL(2u, acceptance_ratios.size());
  // Random number generator
  for (auto i = 0; i < 10; ++i) {
    CHECK(GetRandomNumber() >= 0.0);
    CHECK(GetRandomNumber() <= 1.0);
  }

  // Reset, check the state variables have been reset
  Reset();
  // Cost index
  CHECK_EQUAL(0u, GetCostIndex());
  // Number of dimensions
  CHECK_EQUAL(0u, GetDimensions());
  // Previous best cost
  CHECK_CLOSE(0.0, GetPreviousBestCost(), 1e-12);
  // Step size
  step_sizes = GetStepSizes();
  CHECK(step_sizes.empty());
  // Acceptance ratios
  acceptance_ratios = GetAcceptanceRatios();
  CHECK(acceptance_ratios.empty());
  // Random number generator
  for (auto i = 0; i < 10; ++i) {
    CHECK(GetRandomNumber() >= 0.0);
    CHECK(GetRandomNumber() <= 1.0);
  }
}

TEST(SimulatedAnnealing_ResetAndSolve)
{
  // A 2D problem with default parameters
  SimpleCostFunction cost_func;
  SimulatedAnnealing sa;
  sa.bounds.SetBounds(0, -10.0, 10.0);
  sa.bounds.SetBounds(1, -10.0, 10.0);
  std::vector<double> min_point {1.0, 1.0};
  auto rc = sa.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 0.01);
  CHECK_CLOSE(0.0, min_point[1], 0.01);

  // Reset, then solve a 3D problem with non-default parameter values
  sa.Reset();
  sa.options.SetRandomSeed(5u);
  sa.bounds.SetBounds(0, -5.0, 5.0);
  sa.bounds.SetBounds(1, -5.0, 5.0);
  sa.bounds.SetBounds(2, -5.0, 5.0);
  sa.options.SetNumberOfTemperatureLoops(10);
  sa.options.SetNumberOfCycles(10);
  sa.options.SetStepReductionFactor(0.8);
  sa.options.SetTemperatureReductionFactor(0.8);
  sa.options.SetTemperature(373.0);
  min_point = {1.0, 1.0, 1.0};
  rc = sa.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 0.01);
  CHECK_CLOSE(0.0, min_point[1], 0.01);
  CHECK_CLOSE(0.0, min_point[2], 0.01);

  // Reset and re-solve original problem
  sa.Reset();
  sa.bounds.SetBounds(0, -10.0, 10.0);
  sa.bounds.SetBounds(1, -10.0, 10.0);
  min_point = {1.0, 1.0};
  rc = sa.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 0.01);
  CHECK_CLOSE(0.0, min_point[1], 0.01);
}

TEST_FIXTURE(TestSimulatedAnnealing, TestInitialiseParameters)
{
  SetDimensions(4);
  AccessInitialiseParameters();
  // Step size
  auto step_sizes = GetStepSizes();
  CHECK_EQUAL(4u, step_sizes.size());
  for (auto step : step_sizes) {
    CHECK_CLOSE(1.0, step, 1e-12);
  }
  // Acceptance ratios
  auto acceptance_ratios = GetAcceptanceRatios();
  CHECK_EQUAL(4u, acceptance_ratios.size());
  for (auto accept : acceptance_ratios) {
    CHECK_CLOSE(1.0, accept, 1e-12);
  }
  // Bounds
  CHECK_EQUAL(4u, bounds.GetNumberOfBounds());
  // Random number generator
  for (auto i = 0; i < 10; ++i) {
    CHECK(GetRandomNumber() >= 0.0);
    CHECK(GetRandomNumber() <= 1.0);
  }
}

TEST_FIXTURE(TestSimulatedAnnealing, GenerateTrialPointsTwoDimensions)
{
  SetDimensions(2);
  AccessInitialiseParameters();
  AccessSetGeneratorSeed(5u);
  std::vector<double> trial_point {1.0, 1.0};
  std::vector<double> step_sizes {0.5, 0.5};
  auto reference_point = trial_point;
  AccessGenerateTrialPoint(trial_point, 0);
  CHECK_EQUAL(2u, trial_point.size());
  // Only the 0th index has changed, check the first is still the same
  CHECK_CLOSE(reference_point[1], trial_point[1], 1e-12);
  // Check the 0th point has changed
  CHECK(fabs(reference_point[0] - trial_point[0]) > 1e-12);
  // Check the step is within the allowable step size
  CHECK(fabs(reference_point[0] - trial_point[0]) <
      (reference_point[0] + step_sizes[0]));
  reference_point = trial_point;
  AccessGenerateTrialPoint(trial_point, 1);
  // Now only the 1st index has changed, check the 0th is still the same
  CHECK_CLOSE(reference_point[0], trial_point[0], 1e-12);
  // Check the first point has changed
  CHECK(fabs(reference_point[1] - trial_point[1]) > 1e-12);
  // Check the step is within the allowable step size
  CHECK(fabs(reference_point[1] - trial_point[1]) <
      (reference_point[1] + step_sizes[1]));
  // Check we are still within the bounds
  CHECK(bounds.IsWithinBounds(trial_point));
}

TEST_FIXTURE(TestSimulatedAnnealing, GenerateTrialPointsFourDimensions)
{
  SetDimensions(4);
  AccessInitialiseParameters();
  bounds.SetBounds(0, -5.0, 5.0);
  bounds.SetBounds(1, -5.0, 5.0);
  bounds.SetBounds(2, -5.0, 5.0);
  bounds.SetBounds(3, -5.0, 5.0);
  std::vector<double> trial_point {1.0, 2.0, 3.0, 4.0};
  std::vector<double> step_sizes {0.1, 0.2, 0.3, 0.4};
  auto reference_point = trial_point;
  AccessGenerateTrialPoint(trial_point, 0);
  CHECK_EQUAL(4u, trial_point.size());
  // Only the 0th index has changed, check the others are still the same
  CHECK_CLOSE(reference_point[1], trial_point[1], 1e-12);
  CHECK_CLOSE(reference_point[2], trial_point[2], 1e-12);
  CHECK_CLOSE(reference_point[3], trial_point[3], 1e-12);
  // Check the 0th point has changed
  CHECK(fabs(reference_point[0] - trial_point[0]) > 1e-12);
  // Check the step is within the allowable step size
  CHECK(fabs(reference_point[0] - trial_point[0]) <
      (reference_point[0] + step_sizes[0]));
  // Check we are still within the bounds
  AccessGenerateTrialPoint(trial_point, 3);
  // The 0th & 3rd index has changed, check the others are still the same
  CHECK_CLOSE(reference_point[1], trial_point[1], 1e-12);
  CHECK_CLOSE(reference_point[2], trial_point[2], 1e-12);
  // Check the 3rd point has changed
  CHECK(fabs(reference_point[3] - trial_point[3]) > 1e-12);
  // Check the step is within the allowable step size
  CHECK(fabs(reference_point[3] - trial_point[3]) <
      (reference_point[3] + step_sizes[3]));
  CHECK(bounds.IsWithinBounds(trial_point));
}

TEST_FIXTURE(TestSimulatedAnnealing, TestUpdateStepSizes)
{
  SetDimensions(4);
  AccessInitialiseParameters();
  bounds.SetBounds(0, -5.0, 5.0);
  bounds.SetBounds(1, -5.0, 5.0);
  bounds.SetBounds(2, -5.0, 5.0);
  bounds.SetBounds(3, -5.0, 5.0);
  std::vector<double> acceptance_ratios {9.0, 0.95, 0.2, 0.5};
  SetAcceptanceRatios(acceptance_ratios);
  AccessUpdateStepSizes();
  std::vector<double> step_sizes = GetStepSizes();
  CHECK_CLOSE(10.0, step_sizes[0], 1e-12);
  CHECK_CLOSE(2.75, step_sizes[1], 1e-12);
  CHECK_CLOSE(0.5, step_sizes[2], 1e-12);
  CHECK_CLOSE(1.0, step_sizes[3], 1e-12);
  // Updating the step sizes should reset the acceptance ratios to 1
  acceptance_ratios = GetAcceptanceRatios();
  CHECK_CLOSE(1.0, acceptance_ratios[0], 1e-12);
  CHECK_CLOSE(1.0, acceptance_ratios[1], 1e-12);
  CHECK_CLOSE(1.0, acceptance_ratios[2], 1e-12);
  CHECK_CLOSE(1.0, acceptance_ratios[3], 1e-12);
}

TEST_FIXTURE(TestSimulatedAnnealing, TestResetStepSizes)
{
  SetDimensions(2);
  AccessInitialiseParameters();
  bounds.SetBounds(0, -10.0, 10.0);
  bounds.SetBounds(1, -5.0, 20.0);
  double step_factor = 1.0;
  AccessResetStepSizes(step_factor);
  auto step_sizes = GetStepSizes();
  CHECK_CLOSE(20.0, step_sizes[0], 1e-12);
  CHECK_CLOSE(25.0, step_sizes[1], 1e-12);
}

TEST(SimulatedAnnealing_GetSetParameters)
{
  SimulatedAnnealing sa;
  // Defaults
  CHECK_EQUAL(5, sa.options.GetNumberOfTemperatureLoops());
  CHECK_EQUAL(20, sa.options.GetNumberOfCycles());
  CHECK_CLOSE(0.9, sa.options.GetStepReductionFactor(), 1e-5);
  CHECK_CLOSE(0.5, sa.options.GetTemperatureReductionFactor(), 1e-5);
  CHECK_CLOSE(1000.0, sa.options.GetTemperature(), 1e-5);
  // Set non-default values to parameters
  sa.options.SetNumberOfTemperatureLoops(10);
  CHECK_EQUAL(10, sa.options.GetNumberOfTemperatureLoops());
  sa.options.SetNumberOfCycles(10);
  CHECK_EQUAL(10, sa.options.GetNumberOfCycles());
  sa.options.SetStepReductionFactor(0.8);
  CHECK_CLOSE(0.8, sa.options.GetStepReductionFactor(), 1e-5);
  sa.options.SetTemperatureReductionFactor(0.8);
  CHECK_CLOSE(0.8, sa.options.GetTemperatureReductionFactor(), 1e-5);
  sa.options.SetTemperature(373.0);
  CHECK_CLOSE(373.0, sa.options.GetTemperature(), 1e-5);
  // Reset and check for parameters' default values
  sa.Reset();
  CHECK_EQUAL(5, sa.options.GetNumberOfTemperatureLoops());
  CHECK_EQUAL(20, sa.options.GetNumberOfCycles());
  CHECK_CLOSE(0.9, sa.options.GetStepReductionFactor(), 1e-5);
  CHECK_CLOSE(0.5, sa.options.GetTemperatureReductionFactor(), 1e-5);
  CHECK_CLOSE(1000.0, sa.options.GetTemperature(), 1e-5);
}

TEST(SimulatedAnnealing_FindMinDefaultParameters)
{
  SimulatedAnnealing sa;
  sa.bounds.SetBounds(0, -10.0, 10.0);
  sa.bounds.SetBounds(1, -10.0, 10.0);
  sa.options.SetRandomSeed(5u);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = sa.FindMin(cost_func, min_point);
  // Note: using default bounds would lead to infinite cost
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 0.01);
  CHECK_CLOSE(0.0, min_point[1], 0.01);
}

TEST(SimulatedAnnealing_FindMinUserSetInitialGuess)
{
  SimulatedAnnealing sa;
  sa.bounds.SetBounds(0, -10.0, 10.0);
  sa.bounds.SetBounds(1, -10.0, 10.0);
  sa.options.SetRandomSeed(5u);
  sa.options.SetAddInitialToPopulation(true);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = sa.FindMin(cost_func, min_point);
  // Note: using default bounds would lead to infinite cost
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 0.01);
  CHECK_CLOSE(0.0, min_point[1], 0.01);
}

TEST(SimulatedAnnealing_FindMinUserSetParameters)
{
  SimulatedAnnealing sa;
  sa.bounds.SetBounds(0, -10.0, 10.0);
  sa.bounds.SetBounds(1, -10.0, 10.0);
  sa.options.SetRandomSeed(5u);
  sa.options.SetNumberOfCycles(10);
  sa.options.SetNumberOfTemperatureLoops(6);
  sa.options.SetStepReductionFactor(0.8);
  sa.options.SetTemperatureReductionFactor(0.8);
  sa.options.SetTemperature(343.0);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = sa.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 0.01);
  CHECK_CLOSE(0.0, min_point[1], 0.01);
}

TEST(SimulatedAnnealing_FindMinGenerateInitialGuessWithTrickyCostFunction)
{
  SimulatedAnnealing sa;
  sa.bounds.SetBounds(0, -10.0, 10.0);
  sa.bounds.SetBounds(1, -10.0, 10.0);
  sa.options.SetRandomSeed(5u);
  sa.options.SetAddInitialToPopulation(true);
  TableTopCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = sa.FindMin(cost_func, min_point);
  // Note: using default bounds would lead to infinite cost
  CHECK_EQUAL(1, rc);
  CHECK_CLOSE(1.0, min_point[0], 0.01);
  CHECK_CLOSE(1.0, min_point[1], 0.01);
}

TEST(SimulatedAnnealing_FindMinEmptyCoordinates)
{
  SimulatedAnnealing sa;
  sa.options.SetRandomSeed(5u);
  sa.bounds.SetBounds(0, -10.0, 10.0);
  sa.bounds.SetBounds(1, -10.0, 10.0);
  SimpleCostFunction cost_func;
  std::vector<double> min_point;
  auto rc = sa.FindMin(cost_func, min_point);
  CHECK_EQUAL(-1, rc);
}

TEST(SimulatedAnnealing_FindMinHitMaxFunctionEvaluations)
{
  SimulatedAnnealing sa;
  sa.options.SetRandomSeed(5u);
  sa.bounds.SetBounds(0, -10.0, 10.0);
  sa.bounds.SetBounds(1, -10.0, 10.0);
  sa.options.SetMaxFunctionEvaluations(10u);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = sa.FindMin(cost_func, min_point);
  CHECK_EQUAL(1, rc);
}

TEST(SimulatedAnnealing_FindMinHitMaxIterations)
{
  SimulatedAnnealing sa;
  sa.options.SetRandomSeed(5u);
  sa.bounds.SetBounds(0, -10.0, 10.0);
  sa.bounds.SetBounds(1, -10.0, 10.0);
  sa.options.SetMaxIterations(10u);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = sa.FindMin(cost_func, min_point);
  CHECK_EQUAL(2, rc);
}

TEST(SimulatedAnnealing_FindMinInitialGuessOutOfBounds)
{
  SimulatedAnnealing sa;
  sa.options.SetRandomSeed(5u);
  sa.bounds.SetBounds(0, -10.0, 10.0);
  sa.bounds.SetBounds(1, -10.0, 10.0);
  sa.options.SetAddInitialToPopulation(true);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {-25.0, 60.0};
  auto rc = sa.FindMin(cost_func, min_point);
  CHECK_EQUAL(-2, rc);
}

TEST(SimulatedAnnealing_FindMinInitialGuessInvalidCost)
{
  SimulatedAnnealing sa;
  sa.options.SetRandomSeed(5u);
  sa.bounds.SetBounds(0, -10.0, 10.0);
  sa.bounds.SetBounds(1, -10.0, 10.0);
  sa.options.SetAddInitialToPopulation(true);
  FirstNanCostFunction cost_func;
  std::vector<double> min_point {-1.0, 1.0};
  auto rc = sa.FindMin(cost_func, min_point);
  CHECK_EQUAL(-3, rc);
}


TEST(SimulatedAnnealing_OutputLevelOne)
{
  SimpleCostFunction cost_func;
  SimulatedAnnealing sa;
  sa.options.SetRandomSeed(0u);
  sa.bounds.SetBounds(0, -10.0, 10.0);
  sa.bounds.SetBounds(1, -10.0, 10.0);
  sa.options.SetMaxIterations(5u);
  sa.options.SetOutputLevel(1);
  std::vector<double> min_point {1.0, 1.0};
  auto rc = sa.FindMin(cost_func, min_point);
  CHECK_EQUAL(2, rc);
}

TEST(SimulatedAnnealing_OutputLevelTwo)
{
  SimpleCostFunction cost_func;
  SimulatedAnnealing sa;
  sa.options.SetRandomSeed(0u);
  sa.bounds.SetBounds(0, -10.0, 10.0);
  sa.bounds.SetBounds(1, -10.0, 10.0);
  sa.options.SetMaxIterations(5u);
  sa.options.SetOutputLevel(2);
  std::vector<double> min_point {1.0, 1.0};
  auto rc = sa.FindMin(cost_func, min_point);
  CHECK_EQUAL(2, rc);
}

TEST(SimulatedAnnealing_OutputLevelThree)
{
  SimpleCostFunction cost_func;
  SimulatedAnnealing sa;
  sa.options.SetRandomSeed(0u);
  sa.bounds.SetBounds(0, -10.0, 10.0);
  sa.bounds.SetBounds(1, -10.0, 10.0);
  sa.options.SetMaxIterations(5u);
  sa.options.SetOutputLevel(3);
  std::vector<double> min_point {1.0, 1.0};
  auto rc = sa.FindMin(cost_func, min_point);
  CHECK_EQUAL(2, rc);
}
}  // suite UnitTestSimulatedAnnealing
}  // namespace UnitTests
}  // namspace Unfit

