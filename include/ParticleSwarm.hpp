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
#ifndef UNFIT_INCLUDE_PARTICLESWARM_HPP_
#define UNFIT_INCLUDE_PARTICLESWARM_HPP_

#include <vector>
#include "GenericCostFunction.hpp"
#include "GenericOptimizer.hpp"

namespace Unfit
{
/**
 * \brief A class to implement the Particle Swarm optimization method
 *
 * This class implements the Particle Swarm optimization method. It
 * requires a cost function (written by the user) and an initial guess of the
 * unknown parameters. The initial guess is only needed to obtain the number
 * of unknowns, so make sure it is of the correct length. The algorithm will
 * proceed to try to find a set of parameters that minimize the cost function.
 *
 * The idea behind Particle Swarm comes from Kennedy, Eberhart and Shi.
 * Inspiration for this implementation was derived from the paper entitled
 * "Chaos-enhanced accelerated particle swarm optimization" written by
 * Amir Hossein Gandomi, Gun Jin Yun, Xin-She Yang and Siamak Talatahari.
 *
 * This class also includes the chaotic maps from the same paper in an effort
 * to further enhance convergence. By default the chaotic maps are not enabled.
 * To use them you need to set ps.options.UseAdaptiveParameters(true). Then you
 * can select which map via ps.options.SetStrategy(x). In total, there are 12
 * chaotic maps included, selected via indices ranging from 1 to 12.
 */
class ParticleSwarm : public GenericOptimizer
{
  friend class TestParticleSwarm;
 public:
  /**
   * The constructor sets all of the parameters to their default values. Here
   * there are three parameters that are set to be different to the default
   * values provided in Options.hpp due to their different usage here. They are:
   *
   *   alpha = 1.0 \n
   *   beta  = 0.5 \n
   *   delta = 0.99
   */
  ParticleSwarm();

  /**
   * As we could derive from this class, the destructor should be virtual. In
   * this case the default destructor is just fine as we are putting everything
   * into std::vectors which take care of themselves.
   */
  virtual ~ParticleSwarm();

  /**
   * \brief A method to find the minimum of a model/function using a Particle
   * Swarm approach.
   *
   * This is the main method to call to perform a minimization of the provided
   * cost function. Here a Particle Swarm optimization technique is adopted.
   * There are a few points worth noting before starting the optimization.
   * First, make sure that the length of the coordinates vector that is passed
   * in to this method matches the number of unknowns you expect in your cost
   * function. As the cost function is user supplied there is no way of checking
   * these match. The values you supply are not important as the initial
   * population is randomly generated (unless you invoke the
   * SetAddInitialToPopulation option). The final result (coordinates
   * with best fit) is returned in the coordinates vector. Second, it is
   * advisable to add bounds to each of the variables in the fit. The algorithm
   * generates possible answers between the provided bounds, and if no bounds
   * are provided the bounds are +/- the maximum double the computer can
   * represent. Chances are you know enough about the problem to make more
   * sensible choices. Convergence with +/- 1e8, for example, is still much
   * faster than the default bounds. By default the solution can move outside
   * the bounds that are set. To prevent this, invoke the SetUseHardBounds
   * option. Finally, if you want to use an external population generator, you
   * can use that as the initial population instead of generating one here. You
   * can invoke the SetPopulation method (from the GenericOptimizer class) to
   * achieve this.
   *
   * Intended use:
   *   ParticleSwarm ps; \n
   *   auto rc = ps.FindMin(CostFunction, coordinates);
   *
   * \param CostFunction The user-supplied function to minimise
   * \param coordinates On input must be the same size as the number of
   *        parameters required by the cost function. On exit returns the
   *        coordinates that provided the minimum cost.
   * \return Zero upon success, non-zero upon failure. \n
   *        -2 = User supplied a population that was invalid \n
   *        -1 = Empty coordinate vector \n
   *         1 = Maximum number of function evaluations was reached \n
   *         2 = Maximum number of iterations was reached
   */
  int FindMin(GenericCostFunction &CostFunction, std::vector<double>
      &coordinates) override;

  /**
   * Resets Particle Swarm back to its default state. This is equivalent
   * to destroying the object and creating a new one. The values of alpha, beta
   * and delta are set as per the constructor.
   */
  void Reset() override;

 private:
  /**
   * This method is called by FindMin and contains the main iteration loop. Here
   * is where the actual optimisation takes place. At each iteration a new trial
   * particle is generated for each member of the existing population. If a
   * new particle has a finite cost it replaces the existing particle, otherwise
   * the existing particle is retained. The best particle found so far is also
   * tracked. See the IsConverged documentation from the GenericOptimizer class
   * for convergence information. The iterations will also cease if the maximum
   * number of iterations or cost function evaluations is exceeded. Note that
   * the amplitude of the random component of each trial particle is reduced
   * in each iteration. The alpha parameter scales the random component and at
   * each iteration alpha is updated as alpha*delta (default delta = 0.99).
   *
   * \param CostFunction Returns the residuals (cost) of the model given a set
   *        of parameters
   * \return Zero upon success, non-zero upon failure \n
   *         1 = Maximum number of function evaluations was reached \n
   *         2 = Maximum number of iterations was reached
   */
  int ProcessFindMin(GenericCostFunction &CostFunction);

  /**
   * This method updates the beta parameter via one of 12 different chaotic
   * maps, which the user can select via ps.options.SetStrategy(x). Each map
   * generates a random (chaotic) sequence of numbers within the range 0 to 1.
   * The implemented maps and their strategy number are listed below.
   *
   *   1  - Chebyshev Map \n
   *   2  - Circle Map \n
   *   3  - Gauss Map \n
   *   4  - Intermittency Map \n
   *   5  - Iterative Map \n
   *   6  - Liebovitch Map \n
   *   7  - Logistic Map \n
   *   8  - Piecewise Map \n
   *   9  - Sine Map \n
   *   10 - Singer Map \n
   *   11 - Sinusoidal Map \n
   *   12 - Tent Map
   *
   * \param enhancement_strategy The index of the desired chaotic map (1 to 12)
   */
  void ChaosEnhancement(unsigned enhancement_strategy);

  /**
   * This method generates a trial particle from an existing population member.
   * Given the current particle, P, and the best particle, B, the trial particle
   * is generated via:
   *
   *   trial = (1-beta)*P + beta*B + alpha*P*rand
   *
   * where rand is a random number that follows a normal distribution with mean
   * 0 and standard deviation 1. You can see that beta scales how much the
   * trial particle moves from P towards the best particle B, and alpha scales
   * the size of the random component.
   *
   * \param member The population member for which we will create a new trial
   *        particle
   * \return The new trial particle
   */
  std::vector<double> GenerateTrialParticle(std::size_t member);

  /**
   * This method calls GenerateTrialParticle to create a new candidate particle.
   * It then calculates the cost of the trial particle and if the cost is
   * finite then the new particle replaces the existing particle at position
   * 'member' in the population. If the cost of the trial particle is not finite
   * then the existing particle will remain as part of the population.
   *
   * \param CostFunction The function with which the cost of the new member will
   *        be calculated
   * \param member The member to be updated
   */
  void UpdatePopulationMember(GenericCostFunction &CostFunction,
    std::size_t member);

  /** A vector that contains the best solution found so far*/
  std::vector<double> best_particle_;
  /** A variable to store the number of dimensions*/
  std::size_t dimensions_;
  /** The index of the cost of each vertex*/
  std::size_t cost_;
};

}  // namespace Unfit

#endif

