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
#include <limits>
#include <numeric>
#include <vector>
#include "Bounds.hpp"
#include "GenericCostFunction.hpp"
#include "DifferentialEvolution.hpp"
#include "TestFunctions.hpp"
#include "UnitTest++.h"

namespace Unfit
{
namespace UnitTests
{
SUITE(UnitTestDifferentialEvolution)
{
TEST(DifferentialEvolution_Constructor)
{
  DifferentialEvolution de;
  CHECK_EQUAL(0u, de.bounds.GetNumberOfBounds());
  CHECK_EQUAL(20u, de.options.GetPopulationSize());
  CHECK_EQUAL(0u, de.options.GetRandomSeed());
  CHECK_EQUAL(1u, de.options.GetStrategy());
  CHECK_CLOSE(0.8, de.options.GetWeightingFactor(), 1e-8);
  CHECK_CLOSE(0.9, de.options.GetCrossOver(), 1e-8);
}

TEST(DifferentialEvolution_FindMinDefaultParams)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(DifferentialEvolution_FindMinUserSetPopulationSize)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  de.options.SetPopulationSize(15u);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
  CHECK_EQUAL(15u, de.options.GetPopulationSize());
}

TEST(DifferentialEvolution_FindMinUserSetInvalidPopulationSize)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  de.options.SetPopulationSize(1u);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_EQUAL(6u, de.options.GetPopulationSize());
}

TEST(DifferentialEvolution_FindMinEmptyCoordinates)
{
  DifferentialEvolution de;
  SimpleCostFunction cost_func;
  std::vector<double> min_point;
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(-1, rc);
}

TEST(DifferentialEvolution_FindMinHitMaxIterations)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  de.options.SetMaxIterations(1);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(2, rc);
}

TEST(DifferentialEvolution_FindMinHitMaxFunctionEvaluations)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  de.options.SetMaxFunctionEvaluations(10);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(1, rc);
}

TEST(DifferentialEvolution_FindMinStrategy1)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  de.options.SetStrategy(1);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(DifferentialEvolution_FindMinStrategy2)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  de.options.SetStrategy(2);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(DifferentialEvolution_FindMinStrategy3)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  de.options.SetStrategy(3);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(DifferentialEvolution_FindMinStrategy4)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  de.options.SetStrategy(4);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(DifferentialEvolution_FindMinStrategy5)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  de.options.SetStrategy(5);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(DifferentialEvolution_FindMinStrategy6)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  de.options.SetStrategy(6);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(DifferentialEvolution_FindMinStrategy7)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  de.options.SetStrategy(7);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(DifferentialEvolution_FindMinStrategy8)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  de.options.SetStrategy(8);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(DifferentialEvolution_FindMinStrategy9)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  de.options.SetStrategy(9);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(DifferentialEvolution_FindMinStrategy10)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  de.options.SetStrategy(10);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(DifferentialEvolution_FindMinConvergedOnCost)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  de.options.SetGeometricTolerance(1e-16);
  de.options.SetCostTolerance(1e-6);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {0.0, 0.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(DifferentialEvolution_OutputLevel1)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  de.options.SetMaxIterations(5);
  de.options.SetPopulationSize(10);
  de.options.SetOutputLevel(1);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(2, rc);
}

TEST(DifferentialEvolution_OutputLevel2)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  de.options.SetMaxIterations(5);
  de.options.SetPopulationSize(10);
  de.options.SetOutputLevel(2);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(2, rc);
}

TEST(DifferentialEvolution_OutputLevel3)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  de.options.SetMaxIterations(5);
  de.options.SetPopulationSize(10);
  de.options.SetOutputLevel(3);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(2, rc);
}

TEST(DifferentialEvolution_GetCostGetSolution)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);

  // Best cost
  CHECK_CLOSE(0.0, de.GetCost(), 1e-4);
  CHECK_CLOSE(0.0, de.GetCost(0), 1e-4);
  // Last cost
  CHECK_CLOSE(0.0, de.GetCost(19), 1e-4);
  // Invalid index
  CHECK(!std::isfinite(de.GetCost(20)));

  // First solution
  auto first = de.GetSolution(0);
  CHECK_CLOSE(0.0, first[0], 1e-2);
  CHECK_CLOSE(0.0, first[1], 1e-2);
  // Last solution
  auto last = de.GetSolution(19);
  CHECK_CLOSE(0.0, last[0], 1e-2);
  CHECK_CLOSE(0.0, last[1], 1e-2);
  // Invalid index
  auto invalid = de.GetSolution(20);
  CHECK(invalid.empty());
}

TEST(DifferentialEvolution_Reset)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  CHECK_EQUAL(2u, de.bounds.GetNumberOfBounds());
  de.options.SetPopulationSize(10u);
  CHECK_EQUAL(10u, de.options.GetPopulationSize());
  de.options.SetRandomSeed(42u);
  CHECK_EQUAL(42u, de.options.GetRandomSeed());
  de.options.SetStrategy(6u);
  CHECK_EQUAL(6u, de.options.GetStrategy());
  de.options.SetWeightingFactor(0.75);
  CHECK_CLOSE(0.75, de.options.GetWeightingFactor(), 1e-8);
  de.options.SetCrossOver(0.85);
  CHECK_CLOSE(0.85, de.options.GetCrossOver(), 1e-8);
  de.options.SetAddInitialToPopulation(true);
  CHECK(de.options.GetAddInitialToPopulation());
  de.options.SetUseHardBounds(true);
  CHECK(de.options.GetUseHardBounds());
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);

  de.Reset();
  CHECK_EQUAL(2u, de.bounds.GetNumberOfBounds());
  CHECK_EQUAL(20u, de.options.GetPopulationSize());
  CHECK_EQUAL(0u, de.options.GetRandomSeed());
  CHECK_EQUAL(1u, de.options.GetStrategy());
  CHECK_CLOSE(0.8, de.options.GetWeightingFactor(), 1e-8);
  CHECK_CLOSE(0.9, de.options.GetCrossOver(), 1e-8);
  CHECK(!de.options.GetAddInitialToPopulation());
  CHECK(!de.options.GetUseHardBounds());
}

TEST(DifferentialEvolution_AddInitialGuess)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  de.options.SetAddInitialToPopulation(true);
  SimpleCostFunction cost_func;
  std::vector<double> minpt {1.0, 1.0};
  auto rc = de.FindMin(cost_func, minpt);
  CHECK_EQUAL(0, rc);
}

TEST(DifferentialEvolution_AddInvalidInitialGuess)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  de.options.SetAddInitialToPopulation(true);
  SimpleCostFunction cost_func;
  std::vector<double> minpt {1.0, std::numeric_limits<double>::signaling_NaN()};
  auto rc = de.FindMin(cost_func, minpt);
  CHECK_EQUAL(0, rc);
}

TEST(DifferentialEvolution_FindMinHardBounds)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, 5.0, 10.0);
  de.bounds.SetBounds(1, 5.0, 10.0);
  de.options.SetUseHardBounds(true);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(5.0, min_point[0], 1e-2);
  CHECK_CLOSE(5.0, min_point[1], 1e-2);
}

TEST(DifferentialEvolution_FindMinUserSetPopulation)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> population = {{1.0, 1.0}, {2.0, 2.0},
      {3.0, 3.0}, {4.0, 4.0}, {-1.0, -1.0}, {-2.0, -2.0}, {-3.0, -3.0},
      {-4.0, -4.0}};
  for (auto &member : population) {
    auto residuals = cost_func(member);
    member.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  DifferentialEvolution de;
  de.SetPopulation(population);
  CHECK_EQUAL(8u, de.options.GetPopulationSize());
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 0.1);
  CHECK_CLOSE(0.0, min_point[1], 0.1);
  CHECK_EQUAL(8u, de.options.GetPopulationSize());
}

TEST(DifferentialEvolution_FindMinUserSetPopulationInvalidInitialGuess)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> population = {{1.0, 1.0}, {2.0, 2.0},
      {3.0, 3.0}, {4.0, 4.0}, {-1.0, -1.0}, {-2.0, -2.0}, {-3.0, -3.0},
      {-4.0, -4.0}};
  for (auto &member : population) {
    auto residuals = cost_func(member);
    member.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  DifferentialEvolution de;
  de.SetPopulation(population);
  std::vector<double> min_point {1.0, 1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(-2, rc);
}

TEST(DifferentialEvolution_FindMinUserSetPopulationInvalidPopulation)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> population = {{1.0, 1.0}, {2.0, 2.0},
      {3.0, 3.0}, {4.0, 4.0}, {-1.0, -1.0}, {-2.0, -2.0}, {-3.0, -3.0},
      {-4.0}};
  for (auto &member : population) {
    auto residuals = cost_func(member);
    member.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  DifferentialEvolution de;
  de.SetPopulation(population);
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(-2, rc);
}

TEST(DifferentialEvolution_FindMinUserSetPopulationTooSmall)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> population = {{1.0, 1.0}, {2.0, 2.0}};
  for (auto &member : population) {
    auto residuals = cost_func(member);
    member.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  DifferentialEvolution de;
  de.SetPopulation(population);
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(-2, rc);
}

TEST(DifferentialEvolution_FindMinUserSetPopulationEmptyPopulation)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> population = {};
  DifferentialEvolution de;
  de.SetPopulation(population);
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(-2, rc);
}

TEST(DifferentialEvolution_FindMinMultiThreadedExecution)
{
  DifferentialEvolution de;
  de.bounds.SetBounds(0, -10.0, 10.0);
  de.bounds.SetBounds(1, -10.0, 10.0);
  de.options.SetUseMultiThreaded(true);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = de.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-4);
  CHECK_CLOSE(0.0, min_point[1], 1e-4);
}
}  // suite UnitTestDifferentialEvolution
}  // namespace UnitTests
}  // namespace Unfit

