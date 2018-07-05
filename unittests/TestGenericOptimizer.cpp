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
#include <vector>
#include "GenericOptimizer.hpp"
#include "TestFunctions.hpp"
#include "UnitTest++.h"

namespace Unfit
{
/**
 * A dummy optimizer class that is used to test GenericOptimizer. This is needed
 * because GenericOptimizer has pure virtual methods and therefore cannot be
 * instantiated directly.
 */
class TestGenericOptimizer : public GenericOptimizer
{
 public:
  /**
   * A dummy implementation of FindMin, which is a pure virtual method in the
   * base GenericOptimizer class.
   *
   * \param cost_function The function used to calculate the residuals
   * \param coordinates The initial guess on entry, the result on exit
   * \return an integer; zero for success, one for failure.
   */
  int FindMin(GenericCostFunction &cost_function,
      std::vector<double> &coordinates)
  {
    ++iterations_;
    if (CalculateCost(cost_function, coordinates)) return 0;
    return 1;
  }

  /**
   * A dummy implementation of Reset, which is a pure virtual method in the
   * base GenericOptimizer class.
   */
  void Reset()
  {
    ResetGenericOptimizer();
  }
};

namespace UnitTests
{
SUITE(UnitTestGenericOptimizer)
{
TEST_FIXTURE(TestGenericOptimizer, ConstructionAndReset)
{
  // Check the construction
  CHECK_EQUAL(0u, bounds.GetNumberOfBounds());
  CHECK_EQUAL(10000u, options.GetMaxIterations());
  CHECK(population_.empty());
  CHECK(random_engines_.empty());
  CHECK_EQUAL(0u, function_evaluations_);
  CHECK_EQUAL(0u, iterations_);
  // Change the defaults
  bounds.SetBounds(0, -1.0, 1.0);
  options.SetMaxIterations(100);
  population_ = {{1.0, 1.0}, {2.0, 2.0}};
  GenerateRandomEngines();
  ++function_evaluations_;
  ++iterations_;
  // Check the defaults have changed
  CHECK_EQUAL(1u, bounds.GetNumberOfBounds());
  CHECK_EQUAL(100u, options.GetMaxIterations());
  CHECK_EQUAL(2u, population_.size());
  CHECK_EQUAL(20u, random_engines_.size());
  CHECK_EQUAL(1u, function_evaluations_);
  CHECK_EQUAL(1u, iterations_);
  // Reset and check we are back at the defaults
  ResetGenericOptimizer();
  CHECK_EQUAL(1u, bounds.GetNumberOfBounds());  // Resets the values only
  CHECK_EQUAL(10000u, options.GetMaxIterations());
  CHECK(population_.empty());
  CHECK(random_engines_.empty());
  CHECK_EQUAL(0u, function_evaluations_);
  CHECK_EQUAL(0u, iterations_);
  // Change the defaults again
  bounds.SetBounds(0, -1.0, 1.0);
  options.SetMaxIterations(100);
  population_ = {{1.0, 1.0}, {2.0, 2.0}};
  GenerateRandomEngines();
  ++function_evaluations_;
  ++iterations_;
  // Check the defaults have changed
  CHECK_EQUAL(1u, bounds.GetNumberOfBounds());
  CHECK_EQUAL(100u, options.GetMaxIterations());
  CHECK_EQUAL(2u, population_.size());
  CHECK_EQUAL(20u, random_engines_.size());
  CHECK_EQUAL(1u, function_evaluations_);
  CHECK_EQUAL(1u, iterations_);
  // Reset and check we are back at the defaults
  Reset();
  CHECK_EQUAL(1u, bounds.GetNumberOfBounds());  // Resets the values only
  CHECK_EQUAL(10000u, options.GetMaxIterations());
  CHECK(population_.empty());
  CHECK(random_engines_.empty());
  CHECK_EQUAL(0u, function_evaluations_);
  CHECK_EQUAL(0u, iterations_);
}

TEST(GenericOptimizer_FindMin)
{
  TestGenericOptimizer opt;
  SimpleCostFunction cost_func;
  std::vector<double> coordinates {1.0, 2.0};
  CHECK_EQUAL(0, opt.FindMin(cost_func, coordinates));

  FirstInfCostFunction cost_func2;
  std::vector<double> coordinates2 {3.0, 4.0};
  CHECK_EQUAL(1, opt.FindMin(cost_func2, coordinates2));
}

TEST(GenericOptimizer_GetCost)
{
  TestGenericOptimizer opt;
  std::vector<std::vector<double>> population = {{1.0, 1.0, 1.0},
      {2.0, 2.0, 4.0}, {3.0, 3.0, 9.0}};
  opt.SetPopulation(population);
  CHECK_CLOSE(1.0, opt.GetCost(0), 1e-12);
  CHECK_CLOSE(4.0, opt.GetCost(1), 1e-12);
  CHECK_CLOSE(9.0, opt.GetCost(2), 1e-12);
  CHECK(!std::isfinite(opt.GetCost(3)));
}

TEST(GenericOptimizer_GetPopulationSetPopulation)
{
  TestGenericOptimizer opt;
  std::vector<std::vector<double>> population = {{1.0, 1.0, 1.0},
      {2.0, 2.0, 4.0}, {3.0, 3.0, 9.0}};
  opt.SetPopulation(population);
  CHECK_EQUAL(3u, opt.options.GetPopulationSize());
  auto new_population = opt.GetPopulation();
  CHECK(population.size() == new_population.size());
  CHECK(population[0].size() == new_population[0].size());
  for (auto i = 0u; i < population.size(); ++i) {
    for (auto j = 0u; j < population[i].size(); ++j) {
      CHECK_CLOSE(population[i][j], new_population[i][j], 1e-12);
    }
  }
}

TEST(GenericOptimizer_GetSolution)
{
  TestGenericOptimizer opt;
  std::vector<std::vector<double>> population = {{1.0, 1.0, 1.0},
      {2.0, 2.0, 4.0}, {3.0, 3.0, 9.0}};
  opt.SetPopulation(population);
  for (auto i = 0u; i < population.size(); ++i) {
    auto member = opt.GetSolution(i);
    CHECK_EQUAL(2u, member.size());
    for (auto j = 0u; j < member.size(); ++j) {
      CHECK_CLOSE(population[i][j], member[j], 1e-12);
    }
  }
  CHECK(opt.GetSolution(3).empty());
}

TEST_FIXTURE(TestGenericOptimizer, CalculateCost)
{
  SimpleCostFunction cost_func;
  std::vector<double> coordinates {1.0, 2.0, 3.0, 0.0};
  CHECK(CalculateCost(cost_func, coordinates));
  CHECK_CLOSE(14.0, coordinates.back(), 1e-12);
  coordinates[1] = std::numeric_limits<double>::infinity();
  CHECK(!CalculateCost(cost_func, coordinates));
  FirstNanCostFunction cost_func2;
  CHECK(!CalculateCost(cost_func2, coordinates));
  LastNanCostFunction cost_func3;
  CHECK(!CalculateCost(cost_func3, coordinates));
  FirstInfCostFunction cost_func4;
  CHECK(!CalculateCost(cost_func4, coordinates));
  LastInfCostFunction cost_func5;
  CHECK(!CalculateCost(cost_func5, coordinates));
}

TEST_FIXTURE(TestGenericOptimizer, GeneratePopulationSimpleCostFunction)
{
  unsigned dimensions = 3u;
  options.SetPopulationSize(10);
  bounds.SetBounds(0, 0.0, 10.0);
  bounds.SetBounds(1, 5.0, 10.0);
  bounds.SetBounds(2, -5.0, 0.0);
  SimpleCostFunction cost_func;

  GeneratePopulation(cost_func, dimensions);

  CHECK_EQUAL(10u, options.GetPopulationSize());
  auto population = GetPopulation();
  // Population is the correct size
  CHECK_EQUAL(dimensions + 1, population[0].size());
  CHECK_EQUAL(10u, population.size());
  // Members are within the desired boundaries
  for (auto member : population) {
    CHECK(bounds.IsWithinBounds(member));
  }
  // Members have valid costs
  for (auto member : population) {
    CHECK(std::isfinite(member.back()));
  }
  // Members are not identical
  for (auto i = 0u; i < dimensions; ++i) {
    CHECK(fabs(population[0][i] - population[1][i]) > 1.0e-12);
  }
}

TEST_FIXTURE(TestGenericOptimizer, GeneratePopulationMultiThreadedExecution)
{
  unsigned dimensions = 3u;
  options.SetUseMultiThreaded(true);
  options.SetPopulationSize(10);
  bounds.SetBounds(0, 0.0, 10.0);
  bounds.SetBounds(1, 5.0, 10.0);
  bounds.SetBounds(2, -5.0, 0.0);
  SimpleCostFunction cost_func;

  GeneratePopulation(cost_func, dimensions);

  CHECK_EQUAL(10u, options.GetPopulationSize());
  auto population = GetPopulation();
  // Population is the correct size
  CHECK_EQUAL(dimensions + 1, population[0].size());
  CHECK_EQUAL(10u, population.size());
  // Members are within the desired boundaries
  for (auto member : population) {
    CHECK(bounds.IsWithinBounds(member));
  }
  // Members have valid costs
  for (auto member : population) {
    CHECK(std::isfinite(member.back()));
  }
  // Members are not identical
  for (auto i = 0u; i < dimensions; ++i) {
    CHECK(fabs(population[0][i] - population[1][i]) > 1.0e-12);
  }
}

TEST_FIXTURE(TestGenericOptimizer, GeneratePopulationZeroDimensions)
{
  unsigned dimensions = 0u;
  options.SetPopulationSize(10);
  SimpleCostFunction cost_func;

  GeneratePopulation(cost_func, dimensions);

  CHECK_EQUAL(0u, options.GetPopulationSize());
  auto population = GetPopulation();
  CHECK(population.empty());
}

TEST_FIXTURE(TestGenericOptimizer, GeneratePopulationZeroPopulationSize)
{
  unsigned dimensions = 3u;
  options.SetPopulationSize(0);
  SimpleCostFunction cost_func;

  GeneratePopulation(cost_func, dimensions);

  CHECK_EQUAL(0u, options.GetPopulationSize());
  auto population = GetPopulation();
  CHECK(population.empty());
}

TEST_FIXTURE(TestGenericOptimizer, GeneratePopulationZeroSizeZeroDimensions)
{
  unsigned dimensions = 0u;
  options.SetPopulationSize(0);
  SimpleCostFunction cost_func;

  GeneratePopulation(cost_func, dimensions);

  CHECK_EQUAL(0u, options.GetPopulationSize());
  auto population = GetPopulation();
  CHECK(population.empty());
}

TEST_FIXTURE(TestGenericOptimizer, IsConvergedBestCostConverged)
{
  // Best cost is below cost_tolerance
  std::vector<std::vector<double>> population = {{1.0, 1.0, 0.0},
      {2.0, 2.0, 4.0}, {3.0, 3.0, 9.0}};
  SetPopulation(population);
  CHECK(IsConverged(population[0]));

  // Best cost is above cost_tolerance
  population[0][2] = 1.0;
  SetPopulation(population);
  CHECK(!IsConverged(population[0]));

  // Best cost is negative and above cost_tolerance
  population[0][2] = -1.0;
  SetPopulation(population);
  CHECK(!IsConverged(population[0]));
}

TEST_FIXTURE(TestGenericOptimizer, IsConvergedCostAndGeometryConverged)
{
  // Costs converged but geometry is not
  std::vector<std::vector<double>> population = {{1.0, 1.0, 1.0},
      {2.0, 2.0, 1.0}, {3.0, 3.0, 1.0}};
  SetPopulation(population);
  CHECK(!IsConverged(population[0]));

  // Costs not converged but geometry is
  population = {{1.0, 1.0, 1.0}, {1.0, 1.0, 2.0}, {1.0, 1.0, 3.0}};
  SetPopulation(population);
  CHECK(!IsConverged(population[0]));

  // Costs converged and geometry converged
  population = {{1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}};
  SetPopulation(population);
  CHECK(IsConverged(population[0]));
}

TEST_FIXTURE(TestGenericOptimizer, IsConvergedTruncatedPopulation)
{
  // Costs converged but geometry is not
  std::vector<std::vector<double>> population = {{1.0, 1.0, 1.0},
      {2.0, 2.0, 1.0}, {3.0, 3.0, 5.0}};
  SetPopulation(population);
  CHECK(!IsConverged(population[0], 1));

  // Costs not converged but geometry is
  population = {{1.0, 1.0, 1.0}, {1.0, 1.0, 2.0}, {3.0, 3.0, 5.0}};
  SetPopulation(population);
  CHECK(!IsConverged(population[0], 1));

  // Costs converged and geometry converged
  population = {{1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {3.0, 3.0, 5.0}};
  SetPopulation(population);
  CHECK(IsConverged(population[0], 1));
}

TEST_FIXTURE(TestGenericOptimizer, PrintInitialOutput)
{
  options.SetOutputLevel(0);
  PrintInitialOutput(35.0);
  options.SetOutputLevel(1);
  PrintInitialOutput(36.0);
}

TEST_FIXTURE(TestGenericOptimizer, PrintIterationOutput)
{
  function_evaluations_ = 22;
  iterations_ = 17;
  options.SetOutputLevel(0);
  PrintInitialOutput(35.0);
  options.SetOutputLevel(1);
  PrintInitialOutput(36.0);
  options.SetOutputLevel(2);
  PrintInitialOutput(37.0);
}

TEST_FIXTURE(TestGenericOptimizer, PrintFinalOutput)
{
  std::vector<std::vector<double>> population = {{1.0, 1.0, 1.0},
      {2.0, 2.0, 4.0}, {3.0, 3.0, 9.0}};
  SetPopulation(population);
  options.SetOutputLevel(0);
  PrintFinalOutput();
  options.SetOutputLevel(1);
  PrintFinalOutput();
  options.SetOutputLevel(2);
  PrintFinalOutput();
  options.SetOutputLevel(3);
  PrintFinalOutput();
}

TEST_FIXTURE(TestGenericOptimizer, SortPopulation)
{
  // Already sorted population is unchanged
  std::vector<std::vector<double>> population = {{1.0, 1.0, 1.0},
      {2.0, 2.0, 4.0}, {3.0, 3.0, 9.0}};
  SetPopulation(population);
  SortPopulation();
  auto new_population = GetPopulation();
  for (auto i = 0u; i < population.size(); ++i) {
    for (auto j = 0u; j < population[i].size(); ++j) {
      CHECK_CLOSE(population[i][j], new_population[i][j], 1e-12);
    }
  }

  // Unsorted population with positive costs
  population = {{1.0, 1.0, 9.0}, {2.0, 2.0, 1.0}, {3.0, 3.0, 4.0}};
  SetPopulation(population);
  SortPopulation();
  new_population = GetPopulation();
  CHECK_CLOSE(2.0, new_population[0][0], 1e-12);
  CHECK_CLOSE(2.0, new_population[0][1], 1e-12);
  CHECK_CLOSE(1.0, new_population[0][2], 1e-12);
  CHECK_CLOSE(3.0, new_population[1][0], 1e-12);
  CHECK_CLOSE(3.0, new_population[1][1], 1e-12);
  CHECK_CLOSE(4.0, new_population[1][2], 1e-12);
  CHECK_CLOSE(1.0, new_population[2][0], 1e-12);
  CHECK_CLOSE(1.0, new_population[2][1], 1e-12);
  CHECK_CLOSE(9.0, new_population[2][2], 1e-12);

  // Unsorted population with negative costs
  population = {{1.0, 1.0, -1.0}, {2.0, 2.0, -9.0}, {3.0, 3.0, -4.0}};
  SetPopulation(population);
  SortPopulation();
  new_population = GetPopulation();
  CHECK_CLOSE(2.0, new_population[0][0], 1e-12);
  CHECK_CLOSE(2.0, new_population[0][1], 1e-12);
  CHECK_CLOSE(-9.0, new_population[0][2], 1e-12);
  CHECK_CLOSE(3.0, new_population[1][0], 1e-12);
  CHECK_CLOSE(3.0, new_population[1][1], 1e-12);
  CHECK_CLOSE(-4.0, new_population[1][2], 1e-12);
  CHECK_CLOSE(1.0, new_population[2][0], 1e-12);
  CHECK_CLOSE(1.0, new_population[2][1], 1e-12);
  CHECK_CLOSE(-1.0, new_population[2][2], 1e-12);
}

TEST_FIXTURE(TestGenericOptimizer, IsPopulationBased)
{
  CHECK(GetIsPopulationBased());
  is_population_based_ = false;
  CHECK(!GetIsPopulationBased());
  ResetGenericOptimizer();
  CHECK(!GetIsPopulationBased());
  is_population_based_ = true;
  CHECK(GetIsPopulationBased());
  ResetGenericOptimizer();
  CHECK(GetIsPopulationBased());
}

TEST(GenericOptimizer_GetIterationsFunctionEvaluations)
{
  TestGenericOptimizer opt;
  SimpleCostFunction cost_func;
  std::vector<double> coordinates {1.0, 2.0};

  CHECK_EQUAL(0u, opt.GetNumberOfIterations());
  CHECK_EQUAL(0u, opt.GetNumberOfFunctionEvaluations());
  CHECK_EQUAL(0, opt.FindMin(cost_func, coordinates));
  CHECK_EQUAL(1u, opt.GetNumberOfIterations());
  CHECK_EQUAL(1u, opt.GetNumberOfFunctionEvaluations());
  opt.ResetGenericOptimizer();
  CHECK_EQUAL(0u, opt.GetNumberOfIterations());
  CHECK_EQUAL(0u, opt.GetNumberOfFunctionEvaluations());
}
}  // suite UnitTestGenericOptimizer
}  // namespace UnitTests
}  // namespace Unfit

