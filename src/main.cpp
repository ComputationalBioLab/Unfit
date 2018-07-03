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
#include <iostream>
#include <string>
#include <vector>
#include "Unfit.hpp"
#include "Parabolic.hpp"

/**
 * This file provides an example of how to use Unfit.
 *
 * In this problem, we have some experimental data, y, at observations, x.
 * We wish to find the parameters, a, b, and c, of a parabolic equation:
 *
 *   y(x) = ax^2+bx+c
 *
 * such that y(x) best approximates the experimental data.
 *
 * We will solve the problem using six algorithms implemented in Unfit:
 * Differential Evolution, Genetic Algorithm, Levenberg-Marquardt,
 * Nelder-Mead simplex, Particle Swarm and Simulated Annealing.
 * Both Nelder-Mead and Levenberg-Marquardt are considered
 * as local minimizers and can quickly converge if given a nice cost
 * function and/or a good initial guess. The remaining algorithms are
 * global minimizers and thus, they can be substantially slower. Try
 * Differential Evolution first as it is less sensitive to the input parameters.
 * It is also possible to use the global minimizers to find a reasonable
 * starting point and then use either Nelder-Mead or Levenberg-Marquardt to
 * find the nearest minimum.
 *
 * The correct solution to the problem is approximately:
 *
 *   a = -21.15
 *   b = 61.04
 *   c = -9.61
 */
int main()
{
  // First we set/get our data points
  std::vector<double> x {0.498531, 0.622145, 0.746551, 0.899687, 0.995019,
    1.24803, 1.49695, 1.7464, 1.86737, 1.92478, 2.07206, 2.12789, 2.23212};
  std::vector<double> y {17.0676, 20.0356, 24.0914, 27.5963, 28.9598, 31.9736,
    34.6866, 33.7931, 31.9415, 30.6897, 28.3853, 23.9687, 18.5146};

  // Now that x and y contains the data, we declare the objective function
  Unfit::Examples::Parabolic parabolic_func(x, y);

  // First, try the Nelder-Mead simplex algorithm
  Unfit::NelderMead nm_opt;
  // Choose your initial guess for (a,b,c)
  std::vector<double> init_guess_nm {1.0, 1.0, 1.0};
  //  Call the Nelder-Mead minimization algorithm...
  nm_opt.FindMin(parabolic_func, init_guess_nm);
  //  ... and voila', your optimized parameter set
  std::cout << "NM Result: " << init_guess_nm[0] << " " << init_guess_nm[1]
  << " " << init_guess_nm[2] << std::endl;
  std::cout << "NM Cost:   " << nm_opt.GetCost() << std::endl;
  std::cout << std::endl;

  // Now try the Levenberg-Marquardt algorithm
  Unfit::LevenbergMarquardt lm_opt;
  // Choose the same initial guess
  std::vector<double> init_guess_lm {1.0, 1.0, 1.0};
  //  Call the Levenberg-Marquardt minimization algorithm...
  lm_opt.FindMin(parabolic_func, init_guess_lm);
  //  ... and voila' again, your optimized parameter set
  std::cout << "LM Result: " << init_guess_lm[0] << " " << init_guess_lm[1]
  << " " << init_guess_lm[2] << std::endl;
  std::cout << "LM Cost:   " << lm_opt.GetCost() << std::endl;
  std::cout << std::endl;

  // Now try the Genetic Algorithm approach
  Unfit::GeneticAlgorithm ga_opt;
  // Here the initial guess just needs to have the correct number of entries
  std::vector<double> init_guess_ga {1.0, 1.0, 1.0};
  // For the Genetic Algorithm, we need to set some sensible bounds
  ga_opt.bounds.SetBounds(0, -100.0, 100.0);  // First parameter
  ga_opt.bounds.SetBounds(1, -100.0, 100.0);  // Second parameter
  ga_opt.bounds.SetBounds(2, -100.0, 100.0);  // Third parameter
  //  Call the Genetic Algorithm minimizer...
  ga_opt.FindMin(parabolic_func, init_guess_ga);
  //  ... and voila' again, your optimized parameter set
  std::cout << "GA Result: " << init_guess_ga[0] << " " << init_guess_ga[1]
  << " " << init_guess_ga[2] << std::endl;
  std::cout << "GA Cost:   " << ga_opt.GetCost() << std::endl;
  std::cout << std::endl;

  // Now try the Differential Evolution algorithm
  Unfit::DifferentialEvolution de_opt;
  // Here the initial guess just needs to have the correct number of entries
  std::vector<double> init_guess_de {1.0, 1.0, 1.0};
  // For Differential Evolution, we need to set some sensible bounds
  de_opt.bounds.SetBounds(0, -100.0, 100.0);  // First parameter
  de_opt.bounds.SetBounds(1, -100.0, 100.0);  // Second parameter
  de_opt.bounds.SetBounds(2, -100.0, 100.0);  // Third parameter
  //  Call the Genetic Algorithm minimizer...
  de_opt.FindMin(parabolic_func, init_guess_de);
  //  ... and voila' again, your optimized parameter set
  std::cout << "DE Result: " << init_guess_de[0] << " " << init_guess_de[1]
  << " " << init_guess_de[2] << std::endl;
  std::cout << "DE Cost:   " << de_opt.GetCost() << std::endl;
  std::cout << std::endl;

  // Now try the Particle Swarm algorithm
  Unfit::ParticleSwarm ps_opt;
  // Here the initial guess just needs to have the correct number of entries
  std::vector<double> init_guess_ps {1.0, 1.0, 1.0};
  // For Particle Swarm, we need to set some sensible bounds
  ps_opt.bounds.SetBounds(0, -100.0, 0.0);  // First parameter
  ps_opt.bounds.SetBounds(1, 0.0, 100.0);   // Second parameter
  ps_opt.bounds.SetBounds(2, -10.0, 10.0);  // Third parameter
  //  Call the Genetic Algorithm minimizer...
  ps_opt.FindMin(parabolic_func, init_guess_ps);
  //  ... and voila' again, your optimized parameter set
  std::cout << "PS Result: " << init_guess_ps[0] << " " << init_guess_ps[1]
  << " " << init_guess_ps[2] << std::endl;
  std::cout << "PS Cost:   " << ps_opt.GetCost() << std::endl;
  std::cout << std::endl;

  // Now try the Simulated Annealing algorithm
  Unfit::SimulatedAnnealing sa_opt;
  // Here the initial guess just needs to have the correct number of entries
  std::vector<double> init_guess_sa {1.0, 1.0, 1.0};
  // For Simulated Annealing, we need to set some sensible bounds
  sa_opt.bounds.SetBounds(0, -50.0, 0.0);  // First parameter
  sa_opt.bounds.SetBounds(1, 0.0, 100.0);   // Second parameter
  sa_opt.bounds.SetBounds(2, -10.0, 10.0);  // Third parameter
  //  Call the Genetic Algorithm minimizer...
  sa_opt.FindMin(parabolic_func, init_guess_sa);
  //  ... and voila' again, your optimized parameter set
  std::cout << "SA Result: " << init_guess_sa[0] << " " << init_guess_sa[1]
  << " " << init_guess_sa[2] << std::endl;
  std::cout << "SA Cost:   " << sa_opt.GetCost() << std::endl;
  std::cout << std::endl;
  return 0;
}
