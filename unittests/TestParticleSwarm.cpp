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
#include "Bounds.hpp"
#include "GenericCostFunction.hpp"
#include "ParticleSwarm.hpp"
#include "TestFunctions.hpp"
#include "UnitTest++.h"

namespace Unfit
{
/**
 * This class is designed solely to provide the functionality to access the
 * private methods and member variables of the ParticleSwarm class.
 */
class TestParticleSwarm : public ParticleSwarm
{
 public:
  /**
   * This method provides access to the private method ChaosEnhancement from
   * the ParticleSwarm class for testing purposes.
   *
   * \param enhancement_strategy An index (1-12) to select the chaotic map
   */
  void AccessChaosEnhancement(unsigned enhancement_strategy)
  {
    ChaosEnhancement(enhancement_strategy);
  }

  /**
   * This method provides access to the private method GenerateTrialParticle
   * from the ParticleSwarm class for testing purposes.
   *
   * \param member The population member for which we will create a new trial
   *        particle
   * \return The new trial particle
   */
  std::vector<double> AccessGenerateTrialParticle(std::size_t member)
  {
    return GenerateTrialParticle(member);
  }

  /**
   * This method provides access to the private method UpdatePopulationMember
   * from the ParticleSwarm class for testing purposes.
   *
   * \param CostFunction The function with which the cost of the new member will
   *        be calculated
   * \param member The member to be updated
   */
  void AccessUpdatePopulationMember(GenericCostFunction &CostFunction,
    std::size_t member)
  {
    UpdatePopulationMember(CostFunction, member);
  }

  /**
   * This method allows the best particle to be set directly for testing
   * purposes. The best particle is private to the ParticleSwarm class.
   *
   * \param particle The particle (including cost) to use as the best particle
   */
  void SetBestParticle(std::vector<double> &particle)
  {
    best_particle_ = particle;
  }

  /**
   * This method allows the number of dimensions to be set directly for testing
   * purposes. The number of dimensions is private to the ParticleSwarm class.
   *
   * \param dimensions The required number of dimensions
   */
  void SetDimensions(std::size_t dimensions)
  {
    dimensions_ = dimensions;
  }
};

namespace UnitTests
{
SUITE(UnitTestParticleSwarm)
{
TEST_FIXTURE(TestParticleSwarm, GenerateTrialParticleZeroAlpha)
{
  GenerateRandomEngines();
  population_ = {{1.0, 2.0, 3.0, 0.0}, {4.0, 5.0, 6.0, 2.0},
      {7.0, 8.0, 9.0, 2.0}, {5.0, 4.0, 3.0, 2.0}, {0.0, 1.0, 2.0, 0.0}};
  options.SetAlpha(0.0);  // No random component
  options.SetBeta(0.5);
  SetDimensions(3);
  SetBestParticle(population_[4]);
  const double geom_tol = 1e-8;
  auto new_particle = AccessGenerateTrialParticle(0);
  CHECK_CLOSE(0.5, new_particle[0], geom_tol);
  CHECK_CLOSE(1.5, new_particle[1], geom_tol);
  CHECK_CLOSE(2.5, new_particle[2], geom_tol);

  new_particle = AccessGenerateTrialParticle(1);
  CHECK_CLOSE(2.0, new_particle[0], geom_tol);
  CHECK_CLOSE(3.0, new_particle[1], geom_tol);
  CHECK_CLOSE(4.0, new_particle[2], geom_tol);

  new_particle = AccessGenerateTrialParticle(2);
  CHECK_CLOSE(3.5, new_particle[0], geom_tol);
  CHECK_CLOSE(4.5, new_particle[1], geom_tol);
  CHECK_CLOSE(5.5, new_particle[2], geom_tol);

  new_particle = AccessGenerateTrialParticle(3);
  CHECK_CLOSE(2.5, new_particle[0], geom_tol);
  CHECK_CLOSE(2.5, new_particle[1], geom_tol);
  CHECK_CLOSE(2.5, new_particle[2], geom_tol);

  new_particle = AccessGenerateTrialParticle(4);
  CHECK_CLOSE(0.0, new_particle[0], geom_tol);
  CHECK_CLOSE(1.0, new_particle[1], geom_tol);
  CHECK_CLOSE(2.0, new_particle[2], geom_tol);
}

TEST_FIXTURE(TestParticleSwarm, GenerateTrialParticleWithAlpha)
{
  GenerateRandomEngines();
  population_ = {{1.0, 2.0, 3.0, 0.0}, {4.0, 5.0, 6.0, 2.0},
      {7.0, 8.0, 9.0, 2.0}, {5.0, 4.0, 3.0, 2.0}, {0.0, 1.0, 2.0, 0.0}};
  options.SetAlpha(0.1);
  options.SetBeta(0.5);
  SetDimensions(3);
  SetBestParticle(population_[4]);
  const double geom_tol = 1.5;
  auto new_particle = AccessGenerateTrialParticle(0);
  CHECK_CLOSE(0.5, new_particle[0], geom_tol);
  CHECK_CLOSE(1.5, new_particle[1], geom_tol);
  CHECK_CLOSE(2.5, new_particle[2], geom_tol);

  new_particle = AccessGenerateTrialParticle(1);
  CHECK_CLOSE(2.0, new_particle[0], geom_tol);
  CHECK_CLOSE(3.0, new_particle[1], geom_tol);
  CHECK_CLOSE(4.0, new_particle[2], geom_tol);

  new_particle = AccessGenerateTrialParticle(2);
  CHECK_CLOSE(3.5, new_particle[0], geom_tol);
  CHECK_CLOSE(4.5, new_particle[1], geom_tol);
  CHECK_CLOSE(5.5, new_particle[2], geom_tol);

  new_particle = AccessGenerateTrialParticle(3);
  CHECK_CLOSE(2.5, new_particle[0], geom_tol);
  CHECK_CLOSE(2.5, new_particle[1], geom_tol);
  CHECK_CLOSE(2.5, new_particle[2], geom_tol);

  new_particle = AccessGenerateTrialParticle(4);
  CHECK_CLOSE(0.0, new_particle[0], geom_tol);
  CHECK_CLOSE(1.0, new_particle[1], geom_tol);
  CHECK_CLOSE(2.0, new_particle[2], geom_tol);
}

TEST_FIXTURE(TestParticleSwarm, UpdatePopulationMemberSimpleCostFunction)
{
  SimpleCostFunction cost_func;
  GenerateRandomEngines();
  population_ = {{1.0, 2.0, 3.0, 0.0}, {4.0, 5.0, 6.0, 2.0},
      {7.0, 8.0, 9.0, 2.0}, {5.0, 4.0, 3.0, 2.0}, {0.0, 1.0, 2.0, 0.0}};
  options.SetAlpha(0.0);  // No random component
  options.SetBeta(0.5);
  SetDimensions(3);
  SetBestParticle(population_[4]);
  const double geom_tol = 1e-8;
  AccessUpdatePopulationMember(cost_func, 0);
  CHECK_CLOSE(0.5, population_[0][0], geom_tol);
  CHECK_CLOSE(1.5, population_[0][1], geom_tol);
  CHECK_CLOSE(2.5, population_[0][2], geom_tol);

  AccessUpdatePopulationMember(cost_func, 1);
  CHECK_CLOSE(2.0, population_[1][0], geom_tol);
  CHECK_CLOSE(3.0, population_[1][1], geom_tol);
  CHECK_CLOSE(4.0, population_[1][2], geom_tol);

  AccessUpdatePopulationMember(cost_func, 2);
  CHECK_CLOSE(3.5, population_[2][0], geom_tol);
  CHECK_CLOSE(4.5, population_[2][1], geom_tol);
  CHECK_CLOSE(5.5, population_[2][2], geom_tol);

  AccessUpdatePopulationMember(cost_func, 3);
  CHECK_CLOSE(2.5, population_[3][0], geom_tol);
  CHECK_CLOSE(2.5, population_[3][1], geom_tol);
  CHECK_CLOSE(2.5, population_[3][2], geom_tol);

  AccessUpdatePopulationMember(cost_func, 4);
  CHECK_CLOSE(0.0, population_[4][0], geom_tol);
  CHECK_CLOSE(1.0, population_[4][1], geom_tol);
  CHECK_CLOSE(2.0, population_[4][2], geom_tol);
}

TEST_FIXTURE(TestParticleSwarm, UpdatePopulationMemberInvalidCostFunction)
{
  FirstNanCostFunction cost_func;
  GenerateRandomEngines();
  population_ = {{1.0, 2.0, 3.0, 0.0}, {4.0, 5.0, 6.0, 2.0},
      {7.0, 8.0, 9.0, 2.0}, {5.0, 4.0, 3.0, 2.0}, {0.0, 1.0, 2.0, 0.0}};
  options.SetAlpha(0.0);  // No random component
  options.SetBeta(0.5);
  SetDimensions(3);
  SetBestParticle(population_[4]);
  const double geom_tol = 1e-8;
  AccessUpdatePopulationMember(cost_func, 0);
  CHECK_CLOSE(1.0, population_[0][0], geom_tol);
  CHECK_CLOSE(2.0, population_[0][1], geom_tol);
  CHECK_CLOSE(3.0, population_[0][2], geom_tol);

  AccessUpdatePopulationMember(cost_func, 1);
  CHECK_CLOSE(4.0, population_[1][0], geom_tol);
  CHECK_CLOSE(5.0, population_[1][1], geom_tol);
  CHECK_CLOSE(6.0, population_[1][2], geom_tol);

  AccessUpdatePopulationMember(cost_func, 2);
  CHECK_CLOSE(7.0, population_[2][0], geom_tol);
  CHECK_CLOSE(8.0, population_[2][1], geom_tol);
  CHECK_CLOSE(9.0, population_[2][2], geom_tol);

  AccessUpdatePopulationMember(cost_func, 3);
  CHECK_CLOSE(5.0, population_[3][0], geom_tol);
  CHECK_CLOSE(4.0, population_[3][1], geom_tol);
  CHECK_CLOSE(3.0, population_[3][2], geom_tol);

  AccessUpdatePopulationMember(cost_func, 4);
  CHECK_CLOSE(0.0, population_[4][0], geom_tol);
  CHECK_CLOSE(1.0, population_[4][1], geom_tol);
  CHECK_CLOSE(2.0, population_[4][2], geom_tol);
}

TEST_FIXTURE(TestParticleSwarm, ChaosEnhancementStrategyZeroBetaBounds)
{
  // Check beta is unchanged for ranges between -1.0 and 1.0 (inclusive)
  // and is moved within if it is outside these ranges
  const unsigned strategy = 0;
  const double beta_zero_tol = 1e-12;
  options.SetBeta(-1.0001);
  AccessChaosEnhancement(strategy);
  CHECK_CLOSE(0.9999, options.GetBeta(), beta_zero_tol);
  options.SetBeta(-0.0001);
  AccessChaosEnhancement(strategy);
  CHECK_CLOSE(0.9999, options.GetBeta(), beta_zero_tol);
  options.SetBeta(0.0);
  AccessChaosEnhancement(strategy);
  CHECK_CLOSE(0.0, options.GetBeta(), beta_zero_tol);
  options.SetBeta(0.0001);
  AccessChaosEnhancement(strategy);
  CHECK_CLOSE(0.0001, options.GetBeta(), beta_zero_tol);
  options.SetBeta(1.0);
  AccessChaosEnhancement(strategy);
  CHECK_CLOSE(1.0, options.GetBeta(), beta_zero_tol);
  options.SetBeta(1.0001);
  AccessChaosEnhancement(strategy);
  CHECK_CLOSE(0.0001, options.GetBeta(), beta_zero_tol);
  options.SetBeta(2.0001);
  AccessChaosEnhancement(strategy);
  CHECK_CLOSE(0.0001, options.GetBeta(), beta_zero_tol);
}

TEST_FIXTURE(TestParticleSwarm, ChaosEnhancementStrategyOneSequence)
{
  // Check the first strategy generates a unique sequence of numbers between
  // 0.0 and 1.0, and check the critical values (0, 0.5, 1.0) that
  // can cause this method to stall unless corrected
  const unsigned strategy = 1;

  iterations_ = 1;
  options.SetBeta(0.2);
  for (auto i = 0; i < 10; ++i) {
    ++iterations_;
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  iterations_ = 1;
  options.SetBeta(0.0);
  for (auto i = 0; i < 10; ++i) {
    ++iterations_;
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  iterations_ = 1;
  options.SetBeta(0.5);
  for (auto i = 0; i < 10; ++i) {
    ++iterations_;
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  iterations_ = 1;
  options.SetBeta(1.0);
  for (auto i = 0; i < 10; ++i) {
    ++iterations_;
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }
}

TEST_FIXTURE(TestParticleSwarm, ChaosEnhancementStrategyTwoSequence)
{
  const unsigned strategy = 2;

  options.SetBeta(0.2);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.0);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.5);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(1.0);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }
}

TEST_FIXTURE(TestParticleSwarm, ChaosEnhancementStrategyThreeSequence)
{
  // Beta = 0 is a critical value for this strategy
  const unsigned strategy = 3;

  options.SetBeta(0.2);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.0);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.5);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(1.0);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }
}

TEST_FIXTURE(TestParticleSwarm, ChaosEnhancementStrategyFourSequence)
{
  // Beta = 1 is a critical value for this strategy
  const unsigned strategy = 4;

  options.SetBeta(0.2);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.0);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.5);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(1.0);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }
}

TEST_FIXTURE(TestParticleSwarm, ChaosEnhancementStrategyFiveSequence)
{
  // Beta = 0 is a critical value for this strategy
  const unsigned strategy = 5;

  options.SetBeta(0.2);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.0);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.5);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(1.0);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }
}

TEST_FIXTURE(TestParticleSwarm, ChaosEnhancementStrategySixSequence)
{
  // Beta = 0, 0.5, 1 are critical values for this strategy
  const unsigned strategy = 6;

  options.SetBeta(0.2);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.0);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.5);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(1.0);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }
}

TEST_FIXTURE(TestParticleSwarm, ChaosEnhancementStrategySevenSequence)
{
  const unsigned strategy = 7;

  options.SetBeta(0.2);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.0);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.5);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(1.0);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }
}

TEST_FIXTURE(TestParticleSwarm, ChaosEnhancementStrategyEightSequence)
{
  // Beta = 0, 0.5, 1 are critical values for this strategy
  const unsigned strategy = 8;

  options.SetBeta(0.2);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.0);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.5);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(1.0);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }
}

TEST_FIXTURE(TestParticleSwarm, ChaosEnhancementStrategyNineSequence)
{
  // Beta = 0 is a critical value for this strategy
  const unsigned strategy = 9;

  options.SetBeta(0.2);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.0);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.5);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(1.0);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }
}

TEST_FIXTURE(TestParticleSwarm, ChaosEnhancementStrategyTenSequence)
{
  // Beta = 0 is a critical value for this strategy
  const unsigned strategy = 10;

  options.SetBeta(0.2);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.0);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.5);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(1.0);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }
}

TEST_FIXTURE(TestParticleSwarm, ChaosEnhancementStrategyElevenSequence)
{
  // This map really only produces something useful for beta around 0.5 - 0.7.
  // Away from this, you tend to get sequences of zero (or close to it).
  const unsigned strategy = 11;

  options.SetBeta(0.5);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.6);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.7);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }
}

TEST_FIXTURE(TestParticleSwarm, ChaosEnhancementStrategyTwelveSequence)
{
  // Beta = 0 is a critical value for this strategy
  const unsigned strategy = 12;

  options.SetBeta(0.2);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.0);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(0.5);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }

  options.SetBeta(1.0);
  for (auto i = 0; i < 10; ++i) {
    auto beta = options.GetBeta();
    AccessChaosEnhancement(strategy);
    auto beta_new = options.GetBeta();
    CHECK(fabs(beta - beta_new) > 1e-12);
    CHECK(beta_new >= 0.0);
    CHECK(beta_new <= 1.0);
  }
}

TEST(ParticleSwarm_Constructor)
{
  ParticleSwarm ps;
  CHECK_EQUAL(0u, ps.bounds.GetNumberOfBounds());
  CHECK_EQUAL(20u, ps.options.GetPopulationSize());
  CHECK_EQUAL(0u, ps.options.GetRandomSeed());
  CHECK_EQUAL(1u, ps.options.GetStrategy());
  CHECK_CLOSE(1.0, ps.options.GetAlpha(), 1e-8);
  CHECK_CLOSE(0.5, ps.options.GetBeta(), 1e-8);
  CHECK_CLOSE(0.99, ps.options.GetDelta(), 1e-8);
}

TEST(ParticleSwarm_Reset)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  CHECK_EQUAL(2u, ps.bounds.GetNumberOfBounds());
  ps.options.SetPopulationSize(10u);
  CHECK_EQUAL(10u, ps.options.GetPopulationSize());
  ps.options.SetRandomSeed(42u);
  CHECK_EQUAL(42u, ps.options.GetRandomSeed());
  ps.options.SetStrategy(6u);
  CHECK_EQUAL(6u, ps.options.GetStrategy());
  ps.options.SetAlpha(0.75);
  CHECK_CLOSE(0.75, ps.options.GetAlpha(), 1e-8);
  ps.options.SetBeta(0.85);
  CHECK_CLOSE(0.85, ps.options.GetBeta(), 1e-8);
  ps.options.SetDelta(0.5);
  CHECK_CLOSE(0.5, ps.options.GetDelta(), 1e-8);
  ps.options.SetAddInitialToPopulation(true);
  CHECK(ps.options.GetAddInitialToPopulation());
  ps.options.SetUseHardBounds(true);
  CHECK(ps.options.GetUseHardBounds());
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);

  ps.Reset();
  CHECK_EQUAL(2u, ps.bounds.GetNumberOfBounds());
  CHECK_EQUAL(20u, ps.options.GetPopulationSize());
  CHECK_EQUAL(0u, ps.options.GetRandomSeed());
  CHECK_EQUAL(1u, ps.options.GetStrategy());
  CHECK_CLOSE(1.0, ps.options.GetAlpha(), 1e-8);
  CHECK_CLOSE(0.5, ps.options.GetBeta(), 1e-8);
  CHECK_CLOSE(0.99, ps.options.GetDelta(), 1e-8);
  CHECK(!ps.options.GetAddInitialToPopulation());
  CHECK(!ps.options.GetUseHardBounds());
}

TEST(ParticleSwarm_FindMinDefaultParams)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(ParticleSwarm_FindMinUserSetPopulationSize)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetPopulationSize(15u);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
  CHECK_EQUAL(15u, ps.options.GetPopulationSize());
}

TEST(ParticleSwarm_FindMinUserSetInvalidPopulationSize)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetPopulationSize(1u);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_EQUAL(2u, ps.options.GetPopulationSize());
}

TEST(ParticleSwarm_FindMinEmptyCoordinates)
{
  ParticleSwarm ps;
  SimpleCostFunction cost_func;
  std::vector<double> min_point;
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(-1, rc);
}

TEST(ParticleSwarm_FindMinHitMaxIterations)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetMaxIterations(1);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(2, rc);
}

TEST(ParticleSwarm_FindMinHitMaxFunctionEvaluations)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetMaxFunctionEvaluations(10);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(1, rc);
}

TEST(ParticleSwarm_FindMinConvergedOnCost)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetGeometricTolerance(1e-16);
  ps.options.SetCostTolerance(1e-6);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {0.0, 0.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(ParticleSwarm_OutputLevel1)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetMaxIterations(5);
  ps.options.SetPopulationSize(10);
  ps.options.SetOutputLevel(1);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(2, rc);
}

TEST(ParticleSwarm_OutputLevel2)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetMaxIterations(5);
  ps.options.SetPopulationSize(10);
  ps.options.SetOutputLevel(2);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(2, rc);
}

TEST(ParticleSwarm_OutputLevel3)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetMaxIterations(5);
  ps.options.SetPopulationSize(10);
  ps.options.SetOutputLevel(3);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(2, rc);
}

TEST(ParticleSwarm_GetCostGetSolution)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);

  // Best cost
  CHECK_CLOSE(0.0, ps.GetCost(), 1e-4);
  CHECK_CLOSE(0.0, ps.GetCost(0), 1e-4);
  // Last cost
  CHECK_CLOSE(0.0, ps.GetCost(19), 1e-4);
  // Invalid index
  CHECK(!std::isfinite(ps.GetCost(20)));

  // First solution
  auto first = ps.GetSolution(0);
  CHECK_CLOSE(0.0, first[0], 1e-2);
  CHECK_CLOSE(0.0, first[1], 1e-2);
  // Last solution
  auto last = ps.GetSolution(19);
  CHECK_CLOSE(0.0, last[0], 1e-2);
  CHECK_CLOSE(0.0, last[1], 1e-2);
  // Invalid index
  auto invalid = ps.GetSolution(20);
  CHECK(invalid.empty());
}

TEST(ParticleSwarm_AddInitialGuess)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetAddInitialToPopulation(true);
  SimpleCostFunction cost_func;
  std::vector<double> minpt {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, minpt);
  CHECK_EQUAL(0, rc);
}

TEST(ParticleSwarm_AddInvalidInitialGuess)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetAddInitialToPopulation(true);
  SimpleCostFunction cost_func;
  std::vector<double> minpt {1.0, std::numeric_limits<double>::signaling_NaN()};
  auto rc = ps.FindMin(cost_func, minpt);
  CHECK_EQUAL(0, rc);
}

TEST(ParticleSwarm_FindMinHardBounds)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, 5.0, 10.0);
  ps.bounds.SetBounds(1, 5.0, 10.0);
  ps.options.SetUseHardBounds(true);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(5.0, min_point[0], 1e-2);
  CHECK_CLOSE(5.0, min_point[1], 1e-2);
}

TEST(ParticleSwarm_FindMinUserSetPopulation)
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
  ParticleSwarm ps;
  ps.SetPopulation(population);
  CHECK_EQUAL(8u, ps.options.GetPopulationSize());
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 0.1);
  CHECK_CLOSE(0.0, min_point[1], 0.1);
  CHECK_EQUAL(8u, ps.options.GetPopulationSize());
}

TEST(ParticleSwarm_FindMinUserSetPopulationInvalidInitialGuess)
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
  ParticleSwarm ps;
  ps.SetPopulation(population);
  std::vector<double> min_point {1.0, 1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(-2, rc);
}

TEST(ParticleSwarm_FindMinUserSetPopulationInvalidPopulation)
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
  ParticleSwarm ps;
  ps.SetPopulation(population);
  std::vector<double> min_point {1.0, 1.0};
  auto rc =ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(-2, rc);
}

TEST(ParticleSwarm_FindMinUserSetPopulationTooSmall)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> population = {{1.0, 1.0}};
  for (auto &member : population) {
    auto residuals = cost_func(member);
    member.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  ParticleSwarm ps;
  ps.SetPopulation(population);
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(-2, rc);
}

TEST(ParticleSwarm_FindMinUserSetPopulationEmptyPopulation)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> population = {};
  ParticleSwarm ps;
  ps.SetPopulation(population);
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(-2, rc);
}

TEST(ParticleSwarm_FindMinStrategy0)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetUseAdaptiveParameters(true);
  ps.options.SetStrategy(0);  // No enhancement
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(ParticleSwarm_FindMinStrategy1)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetUseAdaptiveParameters(true);
  ps.options.SetStrategy(1);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(ParticleSwarm_FindMinStrategy2)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetUseAdaptiveParameters(true);
  ps.options.SetStrategy(2);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(ParticleSwarm_FindMinStrategy3)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetUseAdaptiveParameters(true);
  ps.options.SetStrategy(3);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(ParticleSwarm_FindMinStrategy4)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetUseAdaptiveParameters(true);
  ps.options.SetStrategy(4);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(ParticleSwarm_FindMinStrategy5)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetUseAdaptiveParameters(true);
  ps.options.SetStrategy(5);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(ParticleSwarm_FindMinStrategy6)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetUseAdaptiveParameters(true);
  ps.options.SetStrategy(6);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(ParticleSwarm_FindMinStrategy7)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetUseAdaptiveParameters(true);
  ps.options.SetStrategy(7);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(ParticleSwarm_FindMinStrategy8)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetUseAdaptiveParameters(true);
  ps.options.SetStrategy(8);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(ParticleSwarm_FindMinStrategy9)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetUseAdaptiveParameters(true);
  ps.options.SetStrategy(9);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(ParticleSwarm_FindMinStrategy10)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetUseAdaptiveParameters(true);
  ps.options.SetStrategy(10);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(ParticleSwarm_FindMinStrategy11)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetUseAdaptiveParameters(true);
  ps.options.SetStrategy(11);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(ParticleSwarm_FindMinStrategy12)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetUseAdaptiveParameters(true);
  ps.options.SetStrategy(12);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(ParticleSwarm_FindMinStrategyInvalid)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetUseAdaptiveParameters(true);
  ps.options.SetStrategy(20);  // Should work but without enhancement
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}

TEST(ParticleSwarm_FindMinMultiThreadedExecution)
{
  ParticleSwarm ps;
  ps.bounds.SetBounds(0, -10.0, 10.0);
  ps.bounds.SetBounds(1, -10.0, 10.0);
  ps.options.SetUseMultiThreaded(true);
  SimpleCostFunction cost_func;
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ps.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 1e-2);
  CHECK_CLOSE(0.0, min_point[1], 1e-2);
}
}  // suite UnitTestParticleSwarm
}  // namespace UnitTests
}  // namespace Unfit

