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
#include <limits>
#include <numeric>
#include <utility>
#include <vector>
#include "GenericCostFunction.hpp"
#include "GeneticAlgorithm.hpp"
#include "TestFunctions.hpp"
#include "UnitTest++.h"

namespace Unfit
{
/**
 * This class is designed solely to provide the functionality to access the
 * private methods and member variables of the GeneticAlgorithm class.
 */
class TestGeneticAlgorithm : public GeneticAlgorithm
{
 public:
  /**
   * A function to access the private method (GeneratePopulation) in the
   * GeneticAlgorithm class.
   *
   * \param CostFunction The equation to be used to calculate the cost
   */
  void AccessGeneratePopulation(GenericCostFunction &CostFunction)
  {
    GeneratePopulation(CostFunction);
  }

  /**
   * A function to access the private method (GetMatingPairs) in the
   * GeneticAlgorithm class.
   *
   * \return A pair for mating
   */
  std::pair<unsigned, unsigned> AccessGetMatingPair()
  {
    return GetMatingPair();
  }

  /**
   * A function to access the private method (Reproduce) in the
   * GeneticAlgorithm class.
   *
   * \param CostFunction The equation to be used to calculate the cost
   */
  void AccessReproduce(GenericCostFunction &CostFunction)
  {
    Reproduce(CostFunction);
  }

  /**
   * A function to access the private method (CrossOver) in the
   * GeneticAlgorithm class.
   *
   * \param parent_1 The first parent of the mating pair
   * \param parent_2 The second parent of the mating pair
   * \param offspring_1 The first offspring from the mating pair
   * \param offspring_2 The second offspring from the mating pair
   */
  void AccessCrossOver(const std::vector<double> &parent_1,
    const std::vector<double> &parent_2, std::vector<double> &offspring_1,
    std::vector<double> &offspring_2)
  {
    CrossOver(parent_1, parent_2, offspring_1, offspring_2);
  }

  /**
   * A function to access the private method (InitialiseBounds) in the
   * GeneticAlgorithm class.
   */
  void AccessInitialiseBounds()
  {
    InitialiseBounds();
  }

  /**
   * A function to access the private random number generator method in the
   * GeneticAlgorithm class and change its seed.
   *
   * \param seed The new seed
   */
  void AccessGeneratorSeed(unsigned seed)
  {
    generator_.seed(seed);
  }

  /**
   * A function to access the private random number generation method in the
   * GeneticAlgorithm class.
   *
   * \param i The index of the cost function parameter
   * \return A bounded random number for the requested cost function parameter
   */
  double AccessDistributionsGenerator(unsigned i)
  {
    if (i > dimensions_) return std::numeric_limits<double>::infinity();
    return distributions_[i](generator_);
  }

  /**
   * A function to access the private vector (population) in the
   * GeneticAlgorithm class.
   *
   * \return The population of chromosomes (including the appended costs)
   */
  std::vector<std::vector<double>> AccessPopulation() const
  {
    return population_;
  }

  /**
   * A function to access the private method (MutateGenes) in the
   * GeneticAlgorithm class.
   *
   * \param CostFunction The function to calculate the cost of a chromosome
   */
  void AccessMutateGenes(GenericCostFunction &CostFunction)
  {
    MutateGenes(CostFunction);
  }

  /**
   * A function to set the private variable (dimensions_) in the
   * GeneticAlgorithm class.
   *
   * \param dimensions The requested number of dimensions
   */
  void SetDimensions(unsigned dimensions)
  {
    dimensions_ = dimensions;
  }

  /**
   * A function to access the private variable (ranks_) in the
   * GeneticAlgorithm class.
   *
   * \return A vector containing the rank coefficients
   */
  std::vector<double> AccessRanks()
  {
    return ranks_;
  }

  /**
   * A function to set a whole population; for testing purposes only.
   *
   * \param CostFunction Returns the residuals of the model
   * \param chromosomes The intended population
   * \return true if the provided population is useable, false otherwise
   */
  bool AssignPopulation(GenericCostFunction &CostFunction,
      const std::vector<std::vector<double>> &chromosomes)
  {
    if (chromosomes.empty()) return false;
    if (chromosomes[0].empty()) return false;
    options.SetPopulationSize(static_cast<unsigned>(chromosomes.size()));
    SetDimensions(static_cast<unsigned>(chromosomes[0].size()));
    InitialiseBounds();
    CalculateRanks();
    // Check the vertices all lie within the domain, return false if not
    for (auto i : chromosomes) if (!bounds.IsWithinBounds(i)) return false;
    population_ = chromosomes;
    // Calculate the cost of each vertex, return false if invalid
    for (auto &i : population_) {
      i.push_back(0.0);  // default cost
      if (!CalculateCost(CostFunction, i)) return false;
    }
    return true;
  }
};

namespace UnitTests
{
SUITE(UnitTestGeneticAlgorithm)
{
TEST_FIXTURE(TestGeneticAlgorithm, SortPopulationOrdered)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{1.0, 2.0}, {2.0, 3.0},
      {3.0, 4.0}, {4.0, 5.0}};
  AssignPopulation(cost_func, chromosomes);
  SortPopulation();
  auto population = AccessPopulation();
  CHECK_EQUAL(4u, population.size());
  CHECK_EQUAL(3u, population[0].size());
  CHECK_CLOSE(5.0, population[0].back(), 1e-8);
  CHECK_CLOSE(13.0, population[1].back(), 1e-8);
  CHECK_CLOSE(25.0, population[2].back(), 1e-8);
  CHECK_CLOSE(41.0, population[3].back(), 1e-8);
}

TEST_FIXTURE(TestGeneticAlgorithm, SortPopulationReverseOrder)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{4.0, 5.0}, {3.0, 4.0},
      {2.0, 3.0}, {1.0, 2.0}};
  AssignPopulation(cost_func, chromosomes);
  SortPopulation();
  auto population = AccessPopulation();
  CHECK_EQUAL(4u, population.size());
  CHECK_EQUAL(3u, population[0].size());
  CHECK_CLOSE(5.0, population[0].back(), 1e-8);
  CHECK_CLOSE(13.0, population[1].back(), 1e-8);
  CHECK_CLOSE(25.0, population[2].back(), 1e-8);
  CHECK_CLOSE(41.0, population[3].back(), 1e-8);
}

TEST_FIXTURE(TestGeneticAlgorithm, SortPopulationRandomOrder)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}};
  AssignPopulation(cost_func, chromosomes);
  SortPopulation();
  auto population = AccessPopulation();
  CHECK_EQUAL(4u, population.size());
  CHECK_EQUAL(3u, population[0].size());
  CHECK_CLOSE(5.0, population[0].back(), 1e-8);
  CHECK_CLOSE(13.0, population[1].back(), 1e-8);
  CHECK_CLOSE(25.0, population[2].back(), 1e-8);
  CHECK_CLOSE(41.0, population[3].back(), 1e-8);
}

TEST_FIXTURE(TestGeneticAlgorithm, InitialiseBoundsDefaultBounds)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}};
  // This method calls InitialiseBounds, so call this then test it did what
  // was expected
  AssignPopulation(cost_func, chromosomes);
  std::vector<double> lower_bounds;
  std::vector<double> upper_bounds;
  bounds.GetBounds(lower_bounds, upper_bounds);
  // Check the bounds have been correctly set for the two variables
  CHECK_EQUAL(2u, lower_bounds.size());
  CHECK_EQUAL(2u, upper_bounds.size());
  CHECK_CLOSE(-std::numeric_limits<double>::max(), lower_bounds[0], 1e-8);
  CHECK_CLOSE(-std::numeric_limits<double>::max(), lower_bounds[1], 1e-8);
  CHECK_CLOSE(std::numeric_limits<double>::max(), upper_bounds[0], 1e-8);
  CHECK_CLOSE(std::numeric_limits<double>::max(), upper_bounds[1], 1e-8);
}

TEST_FIXTURE(TestGeneticAlgorithm, InitialiseBoundsTightBounds)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}};
  bounds.SetBounds(0, 0.0, 6.0);
  bounds.SetBounds(1, -1.0, 10.0);
  // This method calls InitialiseBounds, so call this then test it did what
  // was expected
  AssignPopulation(cost_func, chromosomes);
  std::vector<double> lower_bounds;
  std::vector<double> upper_bounds;
  bounds.GetBounds(lower_bounds, upper_bounds);
  // Check the bounds have been correctly set for the two variables
  CHECK_EQUAL(2u, lower_bounds.size());
  CHECK_EQUAL(2u, upper_bounds.size());
  CHECK_CLOSE(0.0, lower_bounds[0], 1e-8);
  CHECK_CLOSE(-1.0, lower_bounds[1], 1e-8);
  CHECK_CLOSE(6.0, upper_bounds[0], 1e-8);
  CHECK_CLOSE(10.0, upper_bounds[1], 1e-8);
  // Check that the random number generators for each dimension operate within
  // the stated bounds
  for (auto i = 0u; i < 100; ++i) {
    auto first_num = AccessDistributionsGenerator(0);
    auto second_num = AccessDistributionsGenerator(1);
    CHECK(first_num >= 0.0);
    CHECK(first_num <= 6.0);
    CHECK(second_num >= -1.0);
    CHECK(second_num <= 10.0);
  }
}

TEST_FIXTURE(TestGeneticAlgorithm, MutateGenes)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}, {4.0, 4.0}, {5.0, 3.0}, {6.0, 2.0}, {1.0, 7.0}};
  bounds.SetBounds(0, -1000.0, 1000.0);
  bounds.SetBounds(1, -1000.0, 1000.0);
  AssignPopulation(cost_func, chromosomes);
  options.SetElitism(0);  // All chromosomes can be mutated
  AccessGeneratorSeed(42);
  options.SetGamma(0.2);  // 20% mutations
  AccessMutateGenes(cost_func);

  // The expected number of mutations is Gamma*Dimensions*Population
  // which is 0.2*2*8 = 3 (int). Check we have three mutations.
  auto population = AccessPopulation();
  for (auto i : population) CHECK(std::isfinite(i.back()));  // Check the costs
  for (auto &i : population) i.pop_back();  // Remove the costs
  unsigned number_of_mutations = 0;
  for (auto i = 0u; i < population.size(); ++i) {
    if (fabs(population[i][0]-chromosomes[i][0]) > 1e-8) ++number_of_mutations;
    if (fabs(population[i][1]-chromosomes[i][1]) > 1e-8) ++number_of_mutations;
  }
  CHECK_EQUAL(3u, number_of_mutations);
}

TEST_FIXTURE(TestGeneticAlgorithm, MutateGenesOddNumberChromosomes)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}, {4.0, 4.0}, {5.0, 3.0}, {6.0, 2.0}};
  bounds.SetBounds(0, -1000.0, 1000.0);
  bounds.SetBounds(1, -1000.0, 1000.0);
  AssignPopulation(cost_func, chromosomes);
  options.SetElitism(0);  // All chromosomes can be mutated
  options.SetGamma(0.2);  // 20% mutations
  AccessMutateGenes(cost_func);

  // The expected number of mutations is Gamma*Dimensions*Population
  // which is 0.2*2*7 = 2 (int). Check we have two mutations.
  auto population = AccessPopulation();
  for (auto i : population) CHECK(std::isfinite(i.back()));  // Check the costs
  for (auto &i : population) i.pop_back();  // Remove the costs
  unsigned number_of_mutations = 0;
  for (auto i = 0u; i < population.size(); ++i) {
    if (fabs(population[i][0]-chromosomes[i][0]) > 1e-8) ++number_of_mutations;
    if (fabs(population[i][1]-chromosomes[i][1]) > 1e-8) ++number_of_mutations;
  }
  CHECK_EQUAL(2u, number_of_mutations);
}

TEST_FIXTURE(TestGeneticAlgorithm, MutateGenesWithOneElitism)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}, {4.0, 4.0}, {5.0, 3.0}, {6.0, 2.0}, {1.0, 7.0}};
  bounds.SetBounds(0, -1000.0, 1000.0);
  bounds.SetBounds(1, -1000.0, 1000.0);
  AssignPopulation(cost_func, chromosomes);
  options.SetElitism(1);  // Only 7/8 chromosomes can be mutated
  options.SetGamma(0.2);  // 20% mutations
  AccessMutateGenes(cost_func);

  auto population = AccessPopulation();
  for (auto i : population) CHECK(std::isfinite(i.back()));  // Check the costs
  for (auto &i : population) i.pop_back();  // Remove the costs
  // Check the elite chromosomes are unchanged
  CHECK_CLOSE(population[0][0], chromosomes[0][0], 1e-8);
  CHECK_CLOSE(population[0][1], chromosomes[0][1], 1e-8);
  // The expected number of mutations is Gamma*Dimensions*Population
  // which is 0.2*2*8 = 3 (int). Check we have three mutations.
  unsigned number_of_mutations = 0;
  for (auto i = 1u; i < population.size(); ++i) {
    if (fabs(population[i][0]-chromosomes[i][0]) > 1e-8) ++number_of_mutations;
    if (fabs(population[i][1]-chromosomes[i][1]) > 1e-8) ++number_of_mutations;
  }
  CHECK_EQUAL(3u, number_of_mutations);
}

TEST_FIXTURE(TestGeneticAlgorithm, MutateGenesWithHalfElitism)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}, {4.0, 4.0}, {5.0, 3.0}, {6.0, 2.0}, {1.0, 7.0}};
  bounds.SetBounds(0, -1000.0, 1000.0);
  bounds.SetBounds(1, -1000.0, 1000.0);
  AssignPopulation(cost_func, chromosomes);
  options.SetElitism(4);  // Only 4/8 chromosomes can be mutated
  AccessGeneratorSeed(42);
  options.SetGamma(0.2);  // 20% mutations
  AccessMutateGenes(cost_func);

  auto population = AccessPopulation();
  for (auto i : population) CHECK(std::isfinite(i.back()));  // Check the costs
  for (auto &i : population) i.pop_back();  // Remove the costs
  // Check the elite chromosomes are unchanged
  for (auto i = 0u; i < 4; ++i) {
    CHECK_CLOSE(population[i][0], chromosomes[i][0], 1e-8);
    CHECK_CLOSE(population[i][1], chromosomes[i][1], 1e-8);
  }
  // The expected number of mutations is Gamma*Dimensions*Population
  // which is 0.2*2*8 = 3 (int). Check we have three mutations.
  unsigned number_of_mutations = 0;
  for (auto i = 4u; i < population.size(); ++i) {
    if (fabs(population[i][0]-chromosomes[i][0]) > 1e-8) ++number_of_mutations;
    if (fabs(population[i][1]-chromosomes[i][1]) > 1e-8) ++number_of_mutations;
  }
  CHECK_EQUAL(3u, number_of_mutations);
}

TEST_FIXTURE(TestGeneticAlgorithm, MutateGenesWithAllElite)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}, {4.0, 4.0}, {5.0, 3.0}, {6.0, 2.0}, {1.0, 7.0}};
  bounds.SetBounds(0, -1000.0, 1000.0);
  bounds.SetBounds(1, -1000.0, 1000.0);
  AssignPopulation(cost_func, chromosomes);
  options.SetElitism(8);  // None of chromosomes can be mutated
  options.SetGamma(0.2);  // 20% mutations
  AccessMutateGenes(cost_func);

  auto population = AccessPopulation();
  for (auto i : population) CHECK(std::isfinite(i.back()));  // Check the costs
  for (auto &i : population) i.pop_back();  // Remove the costs
  // Check the elite chromosomes are unchanged
  for (auto i = 0u; i < population.size(); ++i) {
    CHECK_CLOSE(population[i][0], chromosomes[i][0], 1e-8);
    CHECK_CLOSE(population[i][1], chromosomes[i][1], 1e-8);
  }
}

TEST_FIXTURE(TestGeneticAlgorithm, MutateGenesWithAllButOneElite)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}, {4.0, 4.0}, {5.0, 3.0}, {6.0, 2.0}, {1.0, 7.0}};
  bounds.SetBounds(0, -1000.0, 1000.0);
  bounds.SetBounds(1, -1000.0, 1000.0);
  AssignPopulation(cost_func, chromosomes);
  options.SetElitism(7);  // Only 1/8 chromosomes can be mutated
  AccessGeneratorSeed(42);
  options.SetGamma(0.2);  // 20% mutations
  AccessMutateGenes(cost_func);

  auto population = AccessPopulation();
  for (auto i : population) CHECK(std::isfinite(i.back()));  // Check the costs
  for (auto &i : population) i.pop_back();  // Remove the costs
  // Check the elite chromosomes are unchanged
  for (auto i = 0u; i < 7; ++i) {
    CHECK_CLOSE(population[i][0], chromosomes[i][0], 1e-8);
    CHECK_CLOSE(population[i][1], chromosomes[i][1], 1e-8);
  }
  // The expected number of mutations is Gamma*Dimensions*Population
  // which is 0.2*2*8 = 3 (int). However, here we have only two dimensions, so
  // can expect only two mutations
  unsigned number_of_mutations = 0;
  for (auto i = 7u; i < population.size(); ++i) {
    if (fabs(population[i][0]-chromosomes[i][0]) > 1e-8) ++number_of_mutations;
    if (fabs(population[i][1]-chromosomes[i][1]) > 1e-8) ++number_of_mutations;
  }
  CHECK_EQUAL(2u, number_of_mutations);
}

TEST_FIXTURE(TestGeneticAlgorithm, GetMatingPair)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}, {4.0, 4.0}, {5.0, 3.0}, {6.0, 2.0}, {1.0, 7.0}};
  AssignPopulation(cost_func, chromosomes);
  for (auto i = 0u; i < 8; ++i) {
    const auto mating_pair = AccessGetMatingPair();
    CHECK(mating_pair.first < 4);  // 50% survival rate
    CHECK(mating_pair.second < 4);
    CHECK(mating_pair.first != mating_pair.second);  // unigue
  }
}

TEST_FIXTURE(TestGeneticAlgorithm, GetMatingPairKeepMinimum)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}, {4.0, 4.0}, {5.0, 3.0}, {6.0, 2.0}, {1.0, 7.0}};
  options.SetSurvivalRate(0.25);
  AssignPopulation(cost_func, chromosomes);
  for (auto i = 0u; i < 12; ++i) {
    const auto mating_pair = AccessGetMatingPair();
    CHECK(mating_pair.first < 2);  // 0% survival rate, but keep 2
    CHECK(mating_pair.second < 2);
    CHECK(mating_pair.first != mating_pair.second);  // unigue
  }
}

TEST_FIXTURE(TestGeneticAlgorithm, GetMatingPairReplaceNone)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}, {4.0, 4.0}, {5.0, 3.0}, {6.0, 2.0}, {1.0, 7.0}};
  options.SetSurvivalRate(1.0);
  AssignPopulation(cost_func, chromosomes);
  const auto mating_pair = AccessGetMatingPair();
  CHECK_EQUAL(0u, mating_pair.first);  // Default pair generated,
  CHECK_EQUAL(0u, mating_pair.second);  // Default pair generated
}

TEST_FIXTURE(TestGeneticAlgorithm, GetMatingPairReplaceOne)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}, {4.0, 4.0}, {5.0, 3.0}, {6.0, 2.0}, {1.0, 7.0}};
  options.SetSurvivalRate(7.0/8.0);
  AssignPopulation(cost_func, chromosomes);
  for (auto i = 0u; i < 8; ++i) {
    const auto mating_pair = AccessGetMatingPair();
    CHECK(mating_pair.first < 7);  // 87.5% survival rate
    CHECK(mating_pair.second < 7);
    CHECK(mating_pair.first != mating_pair.second);  // unigue
  }
}

TEST_FIXTURE(TestGeneticAlgorithm, GetMatingPairsReplaceAllButOne)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}, {4.0, 4.0}, {5.0, 3.0}, {6.0, 2.0}, {1.0, 7.0}};
  options.SetSurvivalRate(1.0/8.0);
  AssignPopulation(cost_func, chromosomes);
  for (auto i = 0u; i < 8; ++i) {
    const auto mating_pair = AccessGetMatingPair();
    CHECK(mating_pair.first < 2);  // 12.5% survival rate, but keep 2
    CHECK(mating_pair.second < 2);
    CHECK(mating_pair.first != mating_pair.second);  // unigue
  }
}

TEST_FIXTURE(TestGeneticAlgorithm, GetMatingPairsReplaceAll)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}, {4.0, 4.0}, {5.0, 3.0}, {6.0, 2.0}, {1.0, 7.0}};
  options.SetSurvivalRate(0.0/8.0);
  AssignPopulation(cost_func, chromosomes);
  for (auto i = 0u; i < 8; ++i) {
    const auto mating_pair = AccessGetMatingPair();
    CHECK(mating_pair.first < 2);  // 0% survival rate, but keep 2
    CHECK(mating_pair.second < 2);
    CHECK(mating_pair.first != mating_pair.second);  // unigue
  }
}

TEST_FIXTURE(TestGeneticAlgorithm, CrossOver)
{
  std::vector<double> parent_1 {1.0, 2.0, 3.0, 4.0, 5.0, 0.0};
  std::vector<double> parent_2 {6.0, 5.0, 4.0, 3.0, 2.0, 0.0};
  std::vector<double> offspring_1;
  std::vector<double> offspring_2;
  SetDimensions(5u);
  AccessCrossOver(parent_1, parent_2, offspring_1, offspring_2);
  CHECK(offspring_1.size() == parent_1.size());
  CHECK(offspring_2.size() == parent_1.size());
  CHECK_CLOSE(1.0, offspring_1[0], 1e-4);
  CHECK_CLOSE(2.0, offspring_1[1], 1e-4);
  CHECK_CLOSE(3.84427, offspring_1[2], 1e-4);
  CHECK_CLOSE(3.0, offspring_1[3], 1e-4);
  CHECK_CLOSE(2.0, offspring_1[4], 1e-4);
  CHECK_CLOSE(6.0, offspring_2[0], 1e-4);
  CHECK_CLOSE(5.0, offspring_2[1], 1e-4);
  CHECK_CLOSE(3.15573, offspring_2[2], 1e-4);
  CHECK_CLOSE(4.0, offspring_2[3], 1e-4);
  CHECK_CLOSE(5.0, offspring_2[4], 1e-4);

  AccessCrossOver(parent_1, parent_2, offspring_1, offspring_2);
  CHECK(offspring_1.size() == parent_1.size());
  CHECK(offspring_2.size() == parent_1.size());
  CHECK_CLOSE(1.0, offspring_1[0], 1e-4);
  CHECK_CLOSE(2.0, offspring_1[1], 1e-4);
  CHECK_CLOSE(3.0, offspring_1[2], 1e-4);
  CHECK_CLOSE(3.15275, offspring_1[3], 1e-4);
  CHECK_CLOSE(2.0, offspring_1[4], 1e-4);
  CHECK_CLOSE(6.0, offspring_2[0], 1e-4);
  CHECK_CLOSE(5.0, offspring_2[1], 1e-4);
  CHECK_CLOSE(4.0, offspring_2[2], 1e-4);
  CHECK_CLOSE(3.84725, offspring_2[3], 1e-4);
  CHECK_CLOSE(5.0, offspring_2[4], 1e-4);

  AccessCrossOver(parent_1, parent_2, offspring_1, offspring_2);
  CHECK(offspring_1.size() == parent_1.size());
  CHECK(offspring_2.size() == parent_1.size());
  CHECK_CLOSE(1.0, offspring_1[0], 1e-4);
  CHECK_CLOSE(2.0, offspring_1[1], 1e-4);
  CHECK_CLOSE(3.38438, offspring_1[2], 1e-4);
  CHECK_CLOSE(3.0, offspring_1[3], 1e-4);
  CHECK_CLOSE(2.0, offspring_1[4], 1e-4);
  CHECK_CLOSE(6.0, offspring_2[0], 1e-4);
  CHECK_CLOSE(5.0, offspring_2[1], 1e-4);
  CHECK_CLOSE(3.61562, offspring_2[2], 1e-4);
  CHECK_CLOSE(4.0, offspring_2[3], 1e-4);
  CHECK_CLOSE(5.0, offspring_2[4], 1e-4);

  AccessCrossOver(parent_1, parent_2, offspring_1, offspring_2);
  CHECK(offspring_1.size() == parent_1.size());
  CHECK(offspring_2.size() == parent_1.size());
  CHECK_CLOSE(1.0, offspring_1[0], 1e-4);
  CHECK_CLOSE(2.17014, offspring_1[1], 1e-4);
  CHECK_CLOSE(4.0, offspring_1[2], 1e-4);
  CHECK_CLOSE(3.0, offspring_1[3], 1e-4);
  CHECK_CLOSE(2.0, offspring_1[4], 1e-4);
  CHECK_CLOSE(6.0, offspring_2[0], 1e-4);
  CHECK_CLOSE(4.82986, offspring_2[1], 1e-4);
  CHECK_CLOSE(3.0, offspring_2[2], 1e-4);
  CHECK_CLOSE(4.0, offspring_2[3], 1e-4);
  CHECK_CLOSE(5.0, offspring_2[4], 1e-4);
}

TEST_FIXTURE(TestGeneticAlgorithm, GeneratePopulation)
{
  SimpleCostFunction cost_func;
  SetDimensions(3u);
  options.SetPopulationSize(5u);
  bounds.SetBounds(0, -1000.0, 1000.0);
  bounds.SetBounds(1, -1000.0, 1000.0);
  bounds.SetBounds(2, -1000.0, 1000.0);
  AccessInitialiseBounds();
  AccessGeneratePopulation(cost_func);
  auto population = AccessPopulation();
  CHECK_EQUAL(5u, population.size());
  CHECK_EQUAL(4u, population[0].size());
  for (auto i : population) CHECK(std::isfinite(i.back()));  // Check the costs
  for (auto &i : population) i.pop_back();  // Remove the costs
  for (auto i : population) CHECK(bounds.IsWithinBounds(i));  // Check bounds
}

TEST_FIXTURE(TestGeneticAlgorithm, GeneratePopulationCloseBounds)
{
  SimpleCostFunction cost_func;
  options.SetRandomSeed(0u);
  SetDimensions(2u);
  options.SetPopulationSize(3u);
  bounds.SetBounds(0, 0.0, 1.0);
  bounds.SetBounds(1, 1.0, 2.0);
  AccessInitialiseBounds();
  AccessGeneratePopulation(cost_func);
  SortPopulation();
  auto population = AccessPopulation();
  CHECK_EQUAL(3u, population.size());
  CHECK_EQUAL(3u, population[0].size());
  for (auto i : population) CHECK(std::isfinite(i.back()));  // Check the costs
  for (auto &i : population) i.pop_back();  // Remove the costs
  for (auto i : population) CHECK(bounds.IsWithinBounds(i));  // Check bounds
  CHECK_CLOSE(0.623564, population[0][0], 1e-4);
  CHECK_CLOSE(1.38438, population[0][1], 1e-4);
  CHECK_CLOSE(0.592845, population[1][0], 1e-4);
  CHECK_CLOSE(1.84427, population[1][1], 1e-4);
  CHECK_CLOSE(0.857946, population[2][0], 1e-4);
  CHECK_CLOSE(1.84725, population[2][1], 1e-4);
}

TEST_FIXTURE(TestGeneticAlgorithm, GeneratePopulationZeroDimensions)
{
  SimpleCostFunction cost_func;
  SetDimensions(0u);
  options.SetPopulationSize(5u);
  AccessInitialiseBounds();
  AccessGeneratePopulation(cost_func);
  auto population = AccessPopulation();
  CHECK(population.empty());
}

TEST_FIXTURE(TestGeneticAlgorithm, GeneratePopulationZeroPopulationSize)
{
  SimpleCostFunction cost_func;
  SetDimensions(3u);
  options.SetPopulationSize(0u);
  options.SetMaxIterations(1u);
  bounds.SetBounds(0, -1000.0, 1000.0);
  bounds.SetBounds(1, -1000.0, 1000.0);
  bounds.SetBounds(2, -1000.0, 1000.0);
  std::vector<double> min_point {1.0, 2.0, 3.0};
  FindMin(cost_func, min_point);
  auto population = AccessPopulation();
  CHECK_EQUAL(3u, population.size());
}

TEST_FIXTURE(TestGeneticAlgorithm, Reproduce)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{1.0, 2.0}, {2.0, 3.0},
      {3.0, 4.0}, {4.0, 5.0}, {5.0, 6.0}, {6.0, 7.0}, {7.0, 8.0}, {8.0, 9.0}};
  AssignPopulation(cost_func, chromosomes);
  AccessReproduce(cost_func);
  auto population = AccessPopulation();
  CHECK_EQUAL(8u, population.size());
  for (auto &i : population) i.pop_back();  // Remove the costs
  // Check the first half of the population is unchanged
  for (auto i = 0u; i < 4; ++i) {
    CHECK(fabs(population[i][0]-chromosomes[i][0]) < 1e-8);
    CHECK(fabs(population[i][1]-chromosomes[i][1]) < 1e-8);
  }
  // Check the second half of the population is changed
  for (auto i = 4u; i < population.size(); ++i) {
    CHECK(fabs(population[i][0]-chromosomes[i][0]) > 1e-8);
    CHECK(fabs(population[i][1]-chromosomes[i][1]) > 1e-8);
  }
}

TEST_FIXTURE(TestGeneticAlgorithm, ReproduceOddNumberOfChromosomes)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {5.0, 6.0}, {6.0, 7.0}, {7.0, 8.0}, {8.0, 9.0}};
  AssignPopulation(cost_func, chromosomes);
  AccessReproduce(cost_func);
  auto population = AccessPopulation();
  CHECK_EQUAL(7u, population.size());
  for (auto &i : population) i.pop_back();  // Remove the costs
  // Check the first half of the population is unchanged
  for (auto i = 0u; i < 3; ++i) {
    CHECK(fabs(population[i][0]-chromosomes[i][0]) < 1e-8);
    CHECK(fabs(population[i][1]-chromosomes[i][1]) < 1e-8);
  }
  // Check the second half of the population is changed
  for (auto i = 3u; i < population.size(); ++i) {
    CHECK(fabs(population[i][0]-chromosomes[i][0]) > 1e-8);
    CHECK(fabs(population[i][1]-chromosomes[i][1]) > 1e-8);
  }
}

TEST_FIXTURE(TestGeneticAlgorithm, ReproduceSurivialRateTooHigh)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}, {4.0, 4.0}, {5.0, 3.0}, {6.0, 2.0}, {1.0, 7.0}};
  options.SetSurvivalRate(1.0);
  AssignPopulation(cost_func, chromosomes);
  AccessReproduce(cost_func);
  auto population = AccessPopulation();
  CHECK_EQUAL(8u, population.size());
  for (auto &i : population) i.pop_back();  // Remove the costs
  // Check the population is unchanged
  for (auto i = 0u; i < population.size(); ++i) {
    CHECK(fabs(population[i][0]-chromosomes[i][0]) < 1e-8);
    CHECK(fabs(population[i][1]-chromosomes[i][1]) < 1e-8);
  }
}

TEST_FIXTURE(TestGeneticAlgorithm, ReproduceSurivialRateTooLow)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{1.0, 2.0}, {2.0, 3.0},
      {3.0, 4.0}, {4.0, 5.0}, {5.0, 6.0}, {6.0, 7.0}, {7.0, 8.0}, {8.0, 9.0}};
  options.SetSurvivalRate(0.0);
  AssignPopulation(cost_func, chromosomes);
  AccessReproduce(cost_func);
  auto population = AccessPopulation();
  CHECK_EQUAL(8u, population.size());
  for (auto &i : population) i.pop_back();  // Remove the costs
  // Check the first two members of the population is unchanged as there
  // must be two parents for reproduction
  for (auto i = 0u; i < 2; ++i) {
    CHECK(fabs(population[i][0]-chromosomes[i][0]) < 1e-8);
    CHECK(fabs(population[i][1]-chromosomes[i][1]) < 1e-8);
  }
  // Check the remainder of the population is changed
  for (auto i = 2u; i < population.size(); ++i) {
    CHECK(fabs(population[i][0]-chromosomes[i][0]) > 1e-8);
    CHECK(fabs(population[i][1]-chromosomes[i][1]) > 1e-8);
  }
}

TEST_FIXTURE(TestGeneticAlgorithm, GetCost)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}, {4.0, 4.0}, {5.0, 3.0}, {6.0, 2.0}, {1.0, 7.0}};
  AssignPopulation(cost_func, chromosomes);
  // Chromosomes within the defined range
  CHECK_CLOSE(13.0, GetCost(0), 1e-8);
  CHECK_CLOSE(25.0, GetCost(1), 1e-8);
  CHECK_CLOSE(50.0, GetCost(7), 1e-8);
  // Chromosomes outside the defined range
  CHECK(!std::isfinite(GetCost(8)));
  CHECK(!std::isfinite(GetCost(9)));
}

TEST_FIXTURE(TestGeneticAlgorithm, GetSolution)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}, {4.0, 4.0}, {5.0, 3.0}, {6.0, 2.0}, {1.0, 7.0}};
  AssignPopulation(cost_func, chromosomes);
  std::vector<double> individual;
  // Chromosomes within the defined range
  individual = GetSolution(0);
  CHECK_CLOSE(2.0, individual[0], 1e-8);
  CHECK_CLOSE(3.0, individual[1], 1e-8);
  individual = GetSolution(1);
  CHECK_CLOSE(3.0, individual[0], 1e-8);
  CHECK_CLOSE(4.0, individual[1], 1e-8);
  individual = GetSolution(7);
  CHECK_CLOSE(1.0, individual[0], 1e-8);
  CHECK_CLOSE(7.0, individual[1], 1e-8);
  // Chromosomes outside the defined range
  individual = GetSolution(8);
  CHECK(individual.empty());
  individual = GetSolution(9);
  CHECK(individual.empty());
}

TEST_FIXTURE(TestGeneticAlgorithm, IsConvergedCost)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}, {4.0, 4.0}, {5.0, 3.0}, {6.0, 2.0}, {1.0, 7.0}};
  AssignPopulation(cost_func, chromosomes);
  CHECK(!IsConverged(chromosomes[3]));  // Not converged @ default tol
  options.SetCostTolerance(1e6);
  CHECK(IsConverged(chromosomes[3]));   // Converved @ a high tolerance
}

TEST_FIXTURE(TestGeneticAlgorithm, IsConvergedGeometryNotConverged)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}, {4.0, 4.0}, {5.0, 3.0}, {6.0, 2.0}, {1.0, 7.0}};
  AssignPopulation(cost_func, chromosomes);
  CHECK(!IsConverged(chromosomes[3]));
}

TEST_FIXTURE(TestGeneticAlgorithm, CalculateRanks)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}, {4.0, 4.0}, {5.0, 3.0}, {6.0, 2.0}, {1.0, 7.0}};
  AssignPopulation(cost_func, chromosomes);
  const auto ranks = AccessRanks();
  CHECK_EQUAL(4u, ranks.size());
  CHECK_CLOSE(0.4, ranks[0], 1e-8);
  CHECK_CLOSE(0.7, ranks[1], 1e-8);
  CHECK_CLOSE(0.9, ranks[2], 1e-8);
  CHECK_CLOSE(1.0, ranks[3], 1e-8);
}

TEST_FIXTURE(TestGeneticAlgorithm, CalculateRanksAllSurvive)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}, {4.0, 4.0}, {5.0, 3.0}, {6.0, 2.0}, {1.0, 7.0}};
  options.SetSurvivalRate(1.0);
  AssignPopulation(cost_func, chromosomes);
  const auto ranks = AccessRanks();
  CHECK(ranks.empty());
}

TEST_FIXTURE(TestGeneticAlgorithm, CalculateRanksSurvivalRateTooLow)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> chromosomes {{2.0, 3.0}, {3.0, 4.0},
      {4.0, 5.0}, {1.0, 2.0}, {4.0, 4.0}, {5.0, 3.0}, {6.0, 2.0}, {1.0, 7.0}};
  options.SetSurvivalRate(0.0);
  AssignPopulation(cost_func, chromosomes);
  const auto ranks = AccessRanks();
  CHECK_EQUAL(2u, ranks.size());
  CHECK_CLOSE(2.0/3.0, ranks[0], 1e-8);
  CHECK_CLOSE(1.0, ranks[1], 1e-8);
}

TEST(GeneticAlgorithm_FindMin)
{
  SimpleCostFunction cost_func;
  GeneticAlgorithm ga;
  ga.options.SetPopulationSize(20u);
  ga.options.SetElitism(1u);
  ga.options.SetRandomSeed(0u);
  ga.options.SetCostTolerance(1e-4);
  ga.options.SetGamma(0.5);
  ga.bounds.SetBounds(0, -10.0, 10.0);
  ga.bounds.SetBounds(1, -10.0, 10.0);
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ga.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 0.01);
  CHECK_CLOSE(0.0, min_point[1], 0.01);
}

TEST(GeneticAlgorithm_FindMinDefaultParameters)
{
  SimpleCostFunction cost_func;
  GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, -0.25, 0.25);
  ga.bounds.SetBounds(1, -0.25, 0.25);
  std::vector<double> min_point {1.0, 1.0};
  ga.FindMin(cost_func, min_point);
  auto rc = ga.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 0.01);
  CHECK_CLOSE(0.0, min_point[1], 0.01);
}

TEST(GeneticAlgorithm_FindMinMoreDimensions)
{
  SimpleCostFunction cost_func;
  GeneticAlgorithm ga;
  ga.options.SetPopulationSize(30u);
  ga.options.SetElitism(1u);
  ga.options.SetGamma(0.5);
  ga.options.SetRandomSeed(0u);
  ga.options.SetCostTolerance(1e-4);
  ga.bounds.SetBounds(0, -10.0, 10.0);
  ga.bounds.SetBounds(1, -10.0, 10.0);
  ga.bounds.SetBounds(2, -10.0, 10.0);
  std::vector<double> min_point {1.0, 1.0, 1.0};
  auto rc = ga.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 0.01);
  CHECK_CLOSE(0.0, min_point[1], 0.01);
  CHECK_CLOSE(0.0, min_point[2], 0.01);
}

TEST(GeneticAlgorithm_FindMinEmtpyDimensions)
{
  SimpleCostFunction cost_func;
  GeneticAlgorithm ga;
  ga.options.SetPopulationSize(20u);
  ga.options.SetRandomSeed(0u);
  std::vector<double> min_point {};
  auto rc = ga.FindMin(cost_func, min_point);
  CHECK_EQUAL(-1, rc);
}

TEST(GeneticAlgorithm_FindMinExceedsMaximumIterations)
{
  SimpleCostFunction cost_func;
  GeneticAlgorithm ga;
  ga.options.SetPopulationSize(20u);
  ga.options.SetRandomSeed(0u);
  ga.bounds.SetBounds(0, -10.0, 10.0);
  ga.bounds.SetBounds(1, -10.0, 10.0);
  ga.options.SetMaxIterations(10u);
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ga.FindMin(cost_func, min_point);
  CHECK_EQUAL(2, rc);
}

TEST(GeneticAlgorithm_FindMinExceedsMaximumFunctionEvaluations)
{
  SimpleCostFunction cost_func;
  GeneticAlgorithm ga;
  ga.options.SetPopulationSize(20u);
  ga.options.SetRandomSeed(0u);
  ga.bounds.SetBounds(0, -10.0, 10.0);
  ga.bounds.SetBounds(1, -10.0, 10.0);
  ga.options.SetMaxFunctionEvaluations(10u);
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ga.FindMin(cost_func, min_point);
  CHECK_EQUAL(1, rc);
}

TEST(GeneticAlgorithm_OutputLevelOne)
{
  SimpleCostFunction cost_func;
  GeneticAlgorithm ga;
  ga.options.SetPopulationSize(10u);
  ga.options.SetRandomSeed(0u);
  ga.bounds.SetBounds(0, -10.0, 10.0);
  ga.bounds.SetBounds(1, -10.0, 10.0);
  ga.options.SetMaxIterations(5u);
  ga.options.SetOutputLevel(1);
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ga.FindMin(cost_func, min_point);
  CHECK_EQUAL(2, rc);
}

TEST(GeneticAlgorithm_OutputLevelTwo)
{
  SimpleCostFunction cost_func;
  GeneticAlgorithm ga;
  ga.options.SetPopulationSize(10u);
  ga.options.SetRandomSeed(0u);
  ga.bounds.SetBounds(0, -10.0, 10.0);
  ga.bounds.SetBounds(1, -10.0, 10.0);
  ga.options.SetMaxIterations(5u);
  ga.options.SetOutputLevel(2);
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ga.FindMin(cost_func, min_point);
  CHECK_EQUAL(2, rc);
}

TEST(GeneticAlgorithm_OutputLevelThree)
{
  SimpleCostFunction cost_func;
  GeneticAlgorithm ga;
  ga.options.SetPopulationSize(10u);
  ga.options.SetRandomSeed(0u);
  ga.bounds.SetBounds(0, -10.0, 10.0);
  ga.bounds.SetBounds(1, -10.0, 10.0);
  ga.options.SetMaxIterations(5u);
  ga.options.SetOutputLevel(3);
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ga.FindMin(cost_func, min_point);
  CHECK_EQUAL(2, rc);
}

TEST(GeneticAlgorithm_Reset)
{
  // A 2D problem with default parameters
  SimpleCostFunction cost_func;
  GeneticAlgorithm ga;
  ga.options.SetRandomSeed(0u);
  ga.bounds.SetBounds(0, -0.25, 0.25);
  ga.bounds.SetBounds(1, -0.25, 0.25);
  std::vector<double> min_point {1.0, 1.0};
  ga.FindMin(cost_func, min_point);
  auto rc = ga.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 0.01);
  CHECK_CLOSE(0.0, min_point[1], 0.01);

  // Reset, then 3D problem with non-default parameters
  ga.Reset();
  ga.options.SetPopulationSize(30u);
  ga.options.SetElitism(1u);
  ga.options.SetRandomSeed(0u);
  ga.options.SetCostTolerance(1e-4);
  ga.options.SetGamma(0.5);
  ga.bounds.SetBounds(0, -10.0, 10.0);
  ga.bounds.SetBounds(1, -10.0, 10.0);
  ga.bounds.SetBounds(2, -10.0, 10.0);
  min_point = {1.0, 1.0, 1.0};
  rc = ga.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 0.01);
  CHECK_CLOSE(0.0, min_point[1], 0.01);
  CHECK_CLOSE(0.0, min_point[2], 0.01);

  // Reset, then re-solve the original problem
  ga.Reset();
  ga.options.SetRandomSeed(0u);
  ga.bounds.SetBounds(0, -0.25, 0.25);
  ga.bounds.SetBounds(1, -0.25, 0.25);
  min_point = {1.0, 1.0};
  ga.FindMin(cost_func, min_point);
  rc = ga.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 0.01);
  CHECK_CLOSE(0.0, min_point[1], 0.01);
}

TEST(GeneticAlgorithm_GetSetParameters)
{
  GeneticAlgorithm ga;
  // Defaults
  CHECK_EQUAL(20u, ga.options.GetPopulationSize());
  CHECK_EQUAL(1u, ga.options.GetElitism());
  CHECK_CLOSE(0.5, ga.options.GetSurvivalRate(), 1e-8);
  // Set non-default values
  ga.options.SetRandomSeed(10);
  ga.options.SetPopulationSize(50u);
  CHECK_EQUAL(50u, ga.options.GetPopulationSize());
  ga.options.SetElitism(2u);
  CHECK_EQUAL(2u, ga.options.GetElitism());
  ga.options.SetSurvivalRate(0.25);
  CHECK_CLOSE(0.25, ga.options.GetSurvivalRate(), 1e-8);
  // Reset and re-check defaults
  ga.Reset();
  CHECK_EQUAL(20u, ga.options.GetPopulationSize());
  CHECK_EQUAL(1u, ga.options.GetElitism());
  CHECK_CLOSE(0.5, ga.options.GetSurvivalRate(), 1e-8);
}

TEST(GeneticAlgorithm_AddInitialGuess)
{
  GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, -10.0, 10.0);
  ga.bounds.SetBounds(1, -10.0, 10.0);
  ga.options.SetAddInitialToPopulation(true);
  SimpleCostFunction cost_func;
  std::vector<double> minpt {1.0, 1.0};
  auto rc = ga.FindMin(cost_func, minpt);
  CHECK_EQUAL(0, rc);
}

TEST(GeneticAlgorithm_AddInvalidInitialGuess)
{
  GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, -10.0, 10.0);
  ga.bounds.SetBounds(1, -10.0, 10.0);
  ga.options.SetAddInitialToPopulation(true);
  SimpleCostFunction cost_func;
  std::vector<double> minpt {1.0, std::numeric_limits<double>::signaling_NaN()};
  auto rc = ga.FindMin(cost_func, minpt);
  CHECK_EQUAL(0, rc);
}

TEST(GeneticAlgorithm_FindMinUserSetPopulation)
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
  GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, -1.0, 1.0);
  ga.bounds.SetBounds(1, -1.0, 1.0);
  ga.SetPopulation(population);
  CHECK_EQUAL(8u, ga.options.GetPopulationSize());
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ga.FindMin(cost_func, min_point);
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0, min_point[0], 0.1);
  CHECK_CLOSE(0.0, min_point[1], 0.1);
  CHECK_EQUAL(8u, ga.options.GetPopulationSize());
}

TEST(GeneticAlgorithm_FindMinUserSetPopulationInvalidInitialGuess)
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
  GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, -1.0, 1.0);
  ga.bounds.SetBounds(1, -1.0, 1.0);
  ga.SetPopulation(population);
  std::vector<double> min_point {1.0, 1.0, 1.0};
  auto rc = ga.FindMin(cost_func, min_point);
  CHECK_EQUAL(-2, rc);
}

TEST(GeneticAlgorithm_FindMinUserSetPopulationInvalidPopulation)
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
  GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, -1.0, 1.0);
  ga.bounds.SetBounds(1, -1.0, 1.0);
  ga.SetPopulation(population);
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ga.FindMin(cost_func, min_point);
  CHECK_EQUAL(-2, rc);
}

TEST(GeneticAlgorithm_FindMinUserSetPopulationTooSmall)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> population = {{1.0, 1.0}, {2.0, 2.0}};
  for (auto &member : population) {
    auto residuals = cost_func(member);
    member.emplace_back(std::inner_product(begin(residuals), end(residuals),
      begin(residuals), 0.0));
  }
  GeneticAlgorithm ga;
  ga.bounds.SetBounds(0, -1.0, 1.0);
  ga.bounds.SetBounds(1, -1.0, 1.0);
  ga.SetPopulation(population);
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ga.FindMin(cost_func, min_point);
  CHECK_EQUAL(-2, rc);
}

TEST(GeneticAlgorithm_FindMinUserSetPopulationEmptyPopulation)
{
  SimpleCostFunction cost_func;
  std::vector<std::vector<double>> population = {};
  GeneticAlgorithm ga;
  ga.SetPopulation(population);
  std::vector<double> min_point {1.0, 1.0};
  auto rc = ga.FindMin(cost_func, min_point);
  CHECK_EQUAL(-2, rc);
}
}  // suite UnitTestGeneticAlgorithm
}  // namespace UnitTests
}  // namespace Unfit
