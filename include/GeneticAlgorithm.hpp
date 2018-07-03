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
#ifndef UNFIT_INCLUDE_GENETICALGORITHM_HPP_
#define UNFIT_INCLUDE_GENETICALGORITHM_HPP_

#include <random>
#include <utility>
#include <vector>
#include "GenericCostFunction.hpp"
#include "GenericOptimizer.hpp"

namespace Unfit
{
/**
 * \brief A class to implement the Genetic Algorithm optimization method
 *
 * This class implements the Genetic Algorithm optimization method. It
 * requires a cost function (written by the user) and an initial guess of the
 * unknown parameters. The initial guess is only needed to obtain the number
 * of unknowns, so make sure it is of the correct length. The algorithm will
 * proceeds to try to find a set of parameters that minimize the cost function
 *
 * The inspiration for this implementation was derived from the book:
 *   Practical Genetic Algorithms, by Randy L. Haupt, Sue Ellen Haupt.
 * In particular, look at Chapter 3 on continuous GA as this implementation
 * largely follows their recommendations.
 */
class GeneticAlgorithm : public GenericOptimizer
{
  friend class TestGeneticAlgorithm;
 public:
  /**
   * The constructor sets all of the parameters to their default values.
   */
  GeneticAlgorithm();

  /**
   * As we could derive from this class, the destructor should be virtual. In
   * this case the default destructor is just fine as we are putting everything
   * into std::vectors which take care of themselves.
   */
  virtual ~GeneticAlgorithm() {}

  /**
   * \brief A method to find a minimum point of a function using a Genetic
   * Algorithm approach.
   *
   * This is the main method to call to perform a minimization of the provided
   * cost function. Here a Genetic Algorithm optimization technique is adopted.
   * There are a few points worth noting before starting the optimization.
   * First, make sure that the length of the coordinates vector that is passed
   * in to this function matches the number of unknowns you have in your cost
   * function. As the cost function is user supplied there is no way of checking
   * these match up. The values you supply are not important as the initial
   * population is randomly generated. The final result (coordinates with best
   * fit) is returned in the vector. Second, it is advisable to add bounds to
   * each of the variables in the fit. The algorithm generates possible answers
   * between the provided bounds, and if no bounds are provided the bounds are
   * +/- the maximum double the computer can represent. Chances are you know
   * enough about the problem to make more sensible choices. Convergence with
   * +/- 1e8, for example, is still much faster than the default bounds.
   *
   * Intended use:
   *   int rc = object.FindMin(CostFunction, coordinates)
   *
   * \param CostFunction The user-supplied function to minimise
   * \param coordinates On input must be the same size as the number of
   *        parameters in the cost function. On exit returns the unknowns
   *        that provided the minimum cost.
   * \return Zero upon success, non-zero upon failure. \n
   *        -1 = Empty coordinate vector \n
   *         1 = Maximum number of function evaluations was reached \n
   *         2 = Maximum number of iterations was reached
   */
  int FindMin(GenericCostFunction &CostFunction,
      std::vector<double> &coordinates) override;

  /**
   * Resets the GeneticAlgorithm back to its default state. This is equivalent
   * to destroying the object and creating a new one.
   */
  void Reset() override;

 private:
  /**
   * This method generates a population of random guesses of the unknown
   * parameters for the cost function between the bounds.
   * Any guess that does not compute a valid cost is automatically discarded.
   * Updates the population_ member variable with a sorted (based on cost)
   * initial population.
   *
   * \param CostFunction Returns the residuals of the model
   */
  void GeneratePopulation(GenericCostFunction &CostFunction);

  /**
   * Creates the default bounds for all variables (+/- maximum double) that have
   * not previously been set by the user. This method also initialises a random
   * number generator for each variable that generates uniform real numbers
   * between the specified bounds.
   */
  void InitialiseBounds();

  /**
   * CrossOver is called by Reproduce to make a new chromosome (offspring) from
   * a pair of parents. There are many approaches on how to do this. Here we
   * follow the Practical Genetic Algorithms book, swapping where needed and
   * blending at the cross over point.
   *
   * \param parent_1 First parent of the mating pair
   * \param parent_2 Second parent of the mating pair
   * \param offspring_1 First offspring of the parents
   * \param offspring_2 Second offspring of the parents
   */
  void CrossOver(const std::vector<double> &parent_1,
    const std::vector<double> &parent_2, std::vector<double> &offspring_1,
    std::vector<double> &offspring_2);

  /**
   * Reproduce takes the population that survived in the current generation,
   * generates random pairings, and repopulates the population via offspring
   * from these pairs such that the overall population size remains constant
   * across generations. Set/GetSurvivalRate determines how many survive and
   * how many are replaced.
   *
   * \param CostFunction Returns the residuals of the model
   */
  void Reproduce(GenericCostFunction &CostFunction);

  /**
   * GetMatingPair is called by Reproduce. This method randomly selects a pair
   * of surviving chromosomes to be parents for the next generation. At present
   * rank weighting is used, so the chromosomes with the best cost are more
   * likely to be selected. The pair is unique, in that the parent chromosomes
   * are forced to be different.
   *
   * \return The mating pair
   */
  std::pair<unsigned, unsigned> GetMatingPair();

  /**
   * As the name suggests, this method replaces randomly selected genes with
   * random mutations. The number of mutations to be inserted is gamma * the
   * population size. The elite chromosomes (with the lowest costs) are immune
   * from the mutation process (see Set/GetElitism). If a mutation causes an
   * invalid response from the cost function it is discarded and another
   * mutation is attempted.
   *
   * \param CostFunction Returns the residuals of the model
   */
  void MutateGenes(GenericCostFunction &CostFunction);

  /**
   * This implementation of GeneticAlgorithm uses rank weighting to determine
   * which surviving chromosomes make up a pair of parents. This information is
   * static for the life of the problem, and the ranking coefficients are
   * calculated using this method.
   */
  void CalculateRanks();

  /**
   * This method is called by FindMin and contains the main iteration loop. Here
   * is where the actual optimisation takes place, the population reproduces and
   * is mutated, and convergence is determined.
   *
   * \param CostFunction Returns the residuals of the model
   * \return Zero upon success, non-zero upon failure \n
   *         1 = Maximum number of function evaluations was reached \n
   *         2 = Maximum number of iterations was reached
   */
  int ProcessFindMin(GenericCostFunction &CostFunction);

  /** The number of chromosomes retained at each iteration*/
  unsigned n_keep_;
  /** The number of dimensions in the problem*/
  unsigned dimensions_;
  /** Stores the rank of each retained chromosome for rank weighting*/
  std::vector<double> ranks_;
  /** A random number engine for all random needs*/
  std::mt19937 generator_;
  /** Store the distributions for each variable*/
  std::vector<std::uniform_real_distribution<double>> distributions_;
};

}  // namespace Unfit

#endif
