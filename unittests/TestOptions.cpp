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
#include "UnitTest++.h"
#include "Options.hpp"

static const double options_test_tol {1e-16};

namespace Unfit
{
namespace UnitTests
{
SUITE(UnitTestOptions)
{
TEST(Options_GetSetMaximumFunctionEvaluations)
{
  Unfit::Options options;
  CHECK_EQUAL(100000u, options.GetMaxFunctionEvaluations());
  options.SetMaxFunctionEvaluations(10);
  CHECK_EQUAL(10u, options.GetMaxFunctionEvaluations());
  options.SetMaxFunctionEvaluations(5000000);
  CHECK_EQUAL(5000000u, options.GetMaxFunctionEvaluations());
  options.SetMaxFunctionEvaluations(0);
  CHECK_EQUAL(0u, options.GetMaxFunctionEvaluations());
  options.ResetOptions();
  CHECK_EQUAL(100000u, options.GetMaxFunctionEvaluations());
}

TEST(Options_GetSetMaximumIterations)
{
  Unfit::Options options;
  CHECK_EQUAL(10000u, options.GetMaxIterations());
  options.SetMaxIterations(10);
  CHECK_EQUAL(10u, options.GetMaxIterations());
  options.SetMaxIterations(5000000);
  CHECK_EQUAL(5000000u, options.GetMaxIterations());
  options.SetMaxIterations(0);
  CHECK_EQUAL(0u, options.GetMaxIterations());
  options.ResetOptions();
  CHECK_EQUAL(10000u, options.GetMaxIterations());
}

TEST(Options_GetSetOutputLevel)
{
  Unfit::Options options;
  CHECK_EQUAL(0u, options.GetOutputLevel());
  options.SetOutputLevel(1);
  CHECK_EQUAL(1u, options.GetOutputLevel());
  options.SetOutputLevel(5);
  CHECK_EQUAL(5u, options.GetOutputLevel());
  options.ResetOptions();
  CHECK_EQUAL(0u, options.GetOutputLevel());
}

TEST(Options_GetSetTolerances)
{
  Unfit::Options options;
  CHECK_CLOSE(1e-12, options.GetCostTolerance(), options_test_tol);
  CHECK_CLOSE(1e-8, options.GetDegenerateTolerance(), options_test_tol);
  CHECK_CLOSE(1e-4, options.GetGeometricTolerance(), options_test_tol);

  // Check allowed changes work
  options.SetCostTolerance(1e-16);
  options.SetDegenerateTolerance(1e-12);
  options.SetGeometricTolerance(1e-8);
  CHECK_CLOSE(1e-16, options.GetCostTolerance(), options_test_tol);
  CHECK_CLOSE(1e-12, options.GetDegenerateTolerance(), options_test_tol);
  CHECK_CLOSE(1e-8, options.GetGeometricTolerance(), options_test_tol);

  // Check invalid changes do nothing
  options.SetCostTolerance(-1e-16);
  options.SetDegenerateTolerance(-1e-12);
  options.SetGeometricTolerance(-1e-8);
  CHECK_CLOSE(1e-16, options.GetCostTolerance(), options_test_tol);
  CHECK_CLOSE(1e-12, options.GetDegenerateTolerance(), options_test_tol);
  CHECK_CLOSE(1e-8, options.GetGeometricTolerance(), options_test_tol);
  options.SetCostTolerance(0.0);
  options.SetDegenerateTolerance(0.0);
  options.SetGeometricTolerance(0.0);
  CHECK_CLOSE(1e-16, options.GetCostTolerance(), options_test_tol);
  CHECK_CLOSE(1e-12, options.GetDegenerateTolerance(), options_test_tol);
  CHECK_CLOSE(1e-8, options.GetGeometricTolerance(), options_test_tol);

  // Check reset works
  options.ResetOptions();
  CHECK_CLOSE(1e-12, options.GetCostTolerance(), options_test_tol);
  CHECK_CLOSE(1e-8, options.GetDegenerateTolerance(), options_test_tol);
  CHECK_CLOSE(1e-4, options.GetGeometricTolerance(), options_test_tol);
}

TEST(Options_AlphaBetaDeltaGamma)
{
  Unfit::Options options;
  // Check the defaults
  CHECK_CLOSE(1.0, options.GetAlpha(), options_test_tol);
  CHECK_CLOSE(2.0, options.GetBeta(), options_test_tol);
  CHECK_CLOSE(0.5, options.GetDelta(), options_test_tol);
  CHECK_CLOSE(0.5, options.GetGamma(), options_test_tol);

  // Set some different values and check they are set
  options.SetAlpha(1.1);
  options.SetBeta(2.1);
  options.SetDelta(0.6);
  options.SetGamma(0.4);
  CHECK_CLOSE(1.1, options.GetAlpha(), options_test_tol);
  CHECK_CLOSE(2.1, options.GetBeta(), options_test_tol);
  CHECK_CLOSE(0.6, options.GetDelta(), options_test_tol);
  CHECK_CLOSE(0.4, options.GetGamma(), options_test_tol);
}

TEST(Options_NelderMeadAlphaBetaDeltaGamma)
{
  Unfit::Options options;
  double alpha, beta, delta, gamma;
  // Check the defaults
  options.GetNelderMeadStepSizes(alpha, beta, delta, gamma);
  CHECK_CLOSE(1.0, alpha, options_test_tol);
  CHECK_CLOSE(2.0, beta, options_test_tol);
  CHECK_CLOSE(0.5, delta, options_test_tol);
  CHECK_CLOSE(0.5, gamma, options_test_tol);

  // Set some different values and check they are set
  options.SetNelderMeadStepSizes(1.1, 2.1, 0.6, 0.4);
  options.GetNelderMeadStepSizes(alpha, beta, delta, gamma);
  CHECK_CLOSE(1.1, alpha, options_test_tol);
  CHECK_CLOSE(2.1, beta, options_test_tol);
  CHECK_CLOSE(0.6, delta, options_test_tol);
  CHECK_CLOSE(0.4, gamma, options_test_tol);

  // Reset and check they are reset
  options.ResetOptions();
  options.GetNelderMeadStepSizes(alpha, beta, delta, gamma);
  CHECK_CLOSE(1.0, alpha, options_test_tol);
  CHECK_CLOSE(2.0, beta, options_test_tol);
  CHECK_CLOSE(0.5, delta, options_test_tol);
  CHECK_CLOSE(0.5, gamma, options_test_tol);
}

TEST(Options_NelderMeadAlphaBetaDeltaGammaConstraints)
{
  Unfit::Options options;
  double alpha, beta, delta, gamma;
  // Check the defaults
  options.GetNelderMeadStepSizes(alpha, beta, delta, gamma);
  CHECK_CLOSE(1.0, alpha, options_test_tol);
  CHECK_CLOSE(2.0, beta, options_test_tol);
  CHECK_CLOSE(0.5, delta, options_test_tol);
  CHECK_CLOSE(0.5, gamma, options_test_tol);

  // Gamma to high (>=1)or too low (<=0) should result in no changes
  options.SetNelderMeadStepSizes(1.0, 2.0, 0.5, 0.0);
  options.GetNelderMeadStepSizes(alpha, beta, delta, gamma);
  CHECK_CLOSE(1.0, alpha, options_test_tol);
  CHECK_CLOSE(2.0, beta, options_test_tol);
  CHECK_CLOSE(0.5, delta, options_test_tol);
  CHECK_CLOSE(0.5, gamma, options_test_tol);
  options.SetNelderMeadStepSizes(1.0, 2.0, 0.5, -1.0);
  options.GetNelderMeadStepSizes(alpha, beta, delta, gamma);
  CHECK_CLOSE(1.0, alpha, options_test_tol);
  CHECK_CLOSE(2.0, beta, options_test_tol);
  CHECK_CLOSE(0.5, delta, options_test_tol);
  CHECK_CLOSE(0.5, gamma, options_test_tol);
  options.SetNelderMeadStepSizes(1.0, 2.0, 0.5, 1.0);
  options.GetNelderMeadStepSizes(alpha, beta, delta, gamma);
  CHECK_CLOSE(1.0, alpha, options_test_tol);
  CHECK_CLOSE(2.0, beta, options_test_tol);
  CHECK_CLOSE(0.5, delta, options_test_tol);
  CHECK_CLOSE(0.5, gamma, options_test_tol);
  options.SetNelderMeadStepSizes(1.0, 2.0, 0.5, 2.0);
  options.GetNelderMeadStepSizes(alpha, beta, delta, gamma);
  CHECK_CLOSE(1.0, alpha, options_test_tol);
  CHECK_CLOSE(2.0, beta, options_test_tol);
  CHECK_CLOSE(0.5, delta, options_test_tol);
  CHECK_CLOSE(0.5, gamma, options_test_tol);

  // Delta to high (>=1)or too low (<=0) should result in no changes
  options.SetNelderMeadStepSizes(1.0, 2.0, 0.0, 0.5);
  options.GetNelderMeadStepSizes(alpha, beta, delta, gamma);
  CHECK_CLOSE(1.0, alpha, options_test_tol);
  CHECK_CLOSE(2.0, beta, options_test_tol);
  CHECK_CLOSE(0.5, delta, options_test_tol);
  CHECK_CLOSE(0.5, gamma, options_test_tol);
  options.SetNelderMeadStepSizes(1.0, 2.0, -1.0, 0.5);
  options.GetNelderMeadStepSizes(alpha, beta, delta, gamma);
  CHECK_CLOSE(1.0, alpha, options_test_tol);
  CHECK_CLOSE(2.0, beta, options_test_tol);
  CHECK_CLOSE(0.5, delta, options_test_tol);
  CHECK_CLOSE(0.5, gamma, options_test_tol);
  options.SetNelderMeadStepSizes(1.0, 2.0, 1.0, 0.5);
  options.GetNelderMeadStepSizes(alpha, beta, delta, gamma);
  CHECK_CLOSE(1.0, alpha, options_test_tol);
  CHECK_CLOSE(2.0, beta, options_test_tol);
  CHECK_CLOSE(0.5, delta, options_test_tol);
  CHECK_CLOSE(0.5, gamma, options_test_tol);
  options.SetNelderMeadStepSizes(1.0, 2.0, 2.0, 0.5);
  options.GetNelderMeadStepSizes(alpha, beta, delta, gamma);
  CHECK_CLOSE(1.0, alpha, options_test_tol);
  CHECK_CLOSE(2.0, beta, options_test_tol);
  CHECK_CLOSE(0.5, delta, options_test_tol);
  CHECK_CLOSE(0.5, gamma, options_test_tol);

  // Gamma >= alpha should result in no changes
  options.SetNelderMeadStepSizes(0.5, 2.0, 0.5, 0.5);
  options.GetNelderMeadStepSizes(alpha, beta, delta, gamma);
  CHECK_CLOSE(1.0, alpha, options_test_tol);
  CHECK_CLOSE(2.0, beta, options_test_tol);
  CHECK_CLOSE(0.5, delta, options_test_tol);
  CHECK_CLOSE(0.5, gamma, options_test_tol);
  options.SetNelderMeadStepSizes(0.2, 2.0, 0.5, 0.5);
  options.GetNelderMeadStepSizes(alpha, beta, delta, gamma);
  CHECK_CLOSE(1.0, alpha, options_test_tol);
  CHECK_CLOSE(2.0, beta, options_test_tol);
  CHECK_CLOSE(0.5, delta, options_test_tol);
  CHECK_CLOSE(0.5, gamma, options_test_tol);

  // alpha >= beta should result in no changes
  options.SetNelderMeadStepSizes(2.0, 2.0, 0.5, 0.5);
  options.GetNelderMeadStepSizes(alpha, beta, delta, gamma);
  CHECK_CLOSE(1.0, alpha, options_test_tol);
  CHECK_CLOSE(2.0, beta, options_test_tol);
  CHECK_CLOSE(0.5, delta, options_test_tol);
  CHECK_CLOSE(0.5, gamma, options_test_tol);
  options.SetNelderMeadStepSizes(5.0, 2.0, 0.5, 0.5);
  options.GetNelderMeadStepSizes(alpha, beta, delta, gamma);
  CHECK_CLOSE(1.0, alpha, options_test_tol);
  CHECK_CLOSE(2.0, beta, options_test_tol);
  CHECK_CLOSE(0.5, delta, options_test_tol);
  CHECK_CLOSE(0.5, gamma, options_test_tol);
}

TEST(Options_GetSetTau)
{
  Unfit::Options options;
  // Check the default
  CHECK_CLOSE(1e-3, options.GetTau(), options_test_tol);
  // Try to set an invalid value, should be ignored
  options.SetTau(-1.0);
  CHECK_CLOSE(1e-3, options.GetTau(), options_test_tol);
  // Try to set a valid value
  options.SetTau(1.0);
  CHECK_CLOSE(1.0, options.GetTau(), options_test_tol);
  // Try a reset, should go back to default
  options.ResetOptions();
  CHECK_CLOSE(1e-3, options.GetTau(), options_test_tol);
}

TEST(Options_GetSetUseAdaptive)
{
  Unfit::Options options;
  CHECK(!options.GetUseAdaptiveParameters());
  options.SetUseAdaptiveParameters(true);
  CHECK(options.GetUseAdaptiveParameters());
  options.SetUseAdaptiveParameters(false);
  CHECK(!options.GetUseAdaptiveParameters());
  options.SetUseAdaptiveParameters(true);
  CHECK(options.GetUseAdaptiveParameters());
  options.ResetOptions();
  CHECK(!options.GetUseAdaptiveParameters());
}

TEST(Options_GetSetPopulationSize)
{
  Unfit::Options options;
  CHECK_EQUAL(20u, options.GetPopulationSize());
  options.SetPopulationSize(10u);
  CHECK_EQUAL(10u, options.GetPopulationSize());
  options.ResetOptions();
  CHECK_EQUAL(20u, options.GetPopulationSize());
}

TEST(Options_GetSetRandomSeed)
{
  Unfit::Options options;
  CHECK_EQUAL(0u, options.GetRandomSeed());
  options.SetRandomSeed(10u);
  CHECK_EQUAL(10u, options.GetRandomSeed());
  options.ResetOptions();
  CHECK_EQUAL(0u, options.GetRandomSeed());
}

TEST(Options_GetSetStrategy)
{
  Unfit::Options options;
  CHECK_EQUAL(1u, options.GetStrategy());
  options.SetStrategy(5u);
  CHECK_EQUAL(5u, options.GetStrategy());
  options.ResetOptions();
  CHECK_EQUAL(1u, options.GetStrategy());
}

TEST(Options_GetSetElitism)
{
  Unfit::Options options;
  CHECK_EQUAL(1u, options.GetElitism());
  // Valid (default population size is 20)
  options.SetElitism(5u);
  CHECK_EQUAL(5u, options.GetElitism());
  // Invalid (above 20 = ignored)
  options.SetElitism(25u);
  CHECK_EQUAL(5u, options.GetElitism());
  options.ResetOptions();
  CHECK_EQUAL(1u, options.GetElitism());
}

TEST(Options_GetSetWeightingFactor)
{
  Unfit::Options options;
  CHECK_CLOSE(0.8, options.GetWeightingFactor(), options_test_tol);
  options.SetWeightingFactor(0.5);
  CHECK_CLOSE(0.5, options.GetWeightingFactor(), options_test_tol);
  options.ResetOptions();
  CHECK_CLOSE(0.8, options.GetWeightingFactor(), options_test_tol);
}

TEST(Options_GetSetCrossOver)
{
  Unfit::Options options;
  CHECK_CLOSE(0.9, options.GetCrossOver(), options_test_tol);
  // Within valid range
  options.SetCrossOver(0.5);
  CHECK_CLOSE(0.5, options.GetCrossOver(), options_test_tol);
  // Below valid range
  options.SetCrossOver(-1.0);
  CHECK_CLOSE(0.0, options.GetCrossOver(), options_test_tol);
  // Above valid range
  options.SetCrossOver(1.5);
  CHECK_CLOSE(1.0, options.GetCrossOver(), options_test_tol);
  options.ResetOptions();
  CHECK_CLOSE(0.9, options.GetCrossOver(), options_test_tol);
}

TEST(Options_GetSetSurvivalRate)
{
  Unfit::Options options;
  CHECK_CLOSE(0.5, options.GetSurvivalRate(), options_test_tol);
  // Within valid range
  options.SetSurvivalRate(0.9);
  CHECK_CLOSE(0.9, options.GetSurvivalRate(), options_test_tol);
  // Below valid range
  options.SetSurvivalRate(-1.0);
  CHECK_CLOSE(0.0, options.GetSurvivalRate(), options_test_tol);
  // Above valid range
  options.SetSurvivalRate(1.5);
  CHECK_CLOSE(1.0, options.GetSurvivalRate(), options_test_tol);
  options.ResetOptions();
  CHECK_CLOSE(0.5, options.GetSurvivalRate(), options_test_tol);
}

TEST(Options_GetSetAddInitialToPopulation)
{
  Unfit::Options options;
  CHECK(!options.GetAddInitialToPopulation());
  options.SetAddInitialToPopulation(true);
  CHECK(options.GetAddInitialToPopulation());
  options.ResetOptions();
  CHECK(!options.GetAddInitialToPopulation());
}

TEST(Options_GetSetUseHardBounds)
{
  Unfit::Options options;
  CHECK(!options.GetUseHardBounds());
  options.SetUseHardBounds(true);
  CHECK(options.GetUseHardBounds());
  options.ResetOptions();
  CHECK(!options.GetUseHardBounds());
}

TEST(Options_GetSetUseMultiThreaded)
{
  Unfit::Options options;
  CHECK(!options.GetUseMultiThreaded());
  options.SetUseMultiThreaded(true);
  CHECK(options.GetUseMultiThreaded());
  options.ResetOptions();
  CHECK(!options.GetUseMultiThreaded());
}

TEST(Options_SimulatedAnnealingParameters)
{
  Unfit::Options options;
  // Check the defaults
  CHECK_CLOSE(1000.0, options.GetTemperature(), options_test_tol);
  CHECK_CLOSE(0.9, options.GetStepReductionFactor(), options_test_tol);
  CHECK_CLOSE(0.5, options.GetTemperatureReductionFactor(), options_test_tol);
  CHECK_EQUAL(20, options.GetNumberOfCycles());
  CHECK_EQUAL(5, options.GetNumberOfTemperatureLoops());

  // Set non-default values to parameters
  options.SetTemperature(373.0);
  CHECK_CLOSE(373.0, options.GetTemperature(), options_test_tol);
  options.SetStepReductionFactor(0.8);
  CHECK_CLOSE(0.8, options.GetStepReductionFactor(), options_test_tol);
  options.SetTemperatureReductionFactor(0.8);
  CHECK_CLOSE(0.8, options.GetTemperatureReductionFactor(), options_test_tol);
  options.SetNumberOfCycles(10);
  CHECK_EQUAL(10, options.GetNumberOfCycles());
  options.SetNumberOfTemperatureLoops(10);
  CHECK_EQUAL(10, options.GetNumberOfTemperatureLoops());

  // Reset and check for parameters' default values
  options.ResetOptions();
  CHECK_CLOSE(1000.0, options.GetTemperature(), options_test_tol);
  CHECK_CLOSE(0.9, options.GetStepReductionFactor(), options_test_tol);
  CHECK_CLOSE(0.5, options.GetTemperatureReductionFactor(), options_test_tol);
  CHECK_EQUAL(20, options.GetNumberOfCycles());
  CHECK_EQUAL(5, options.GetNumberOfTemperatureLoops());
}
}  // suite UnitTestOptions
}  // namespace UnitTests
}  // namespace Unfit
