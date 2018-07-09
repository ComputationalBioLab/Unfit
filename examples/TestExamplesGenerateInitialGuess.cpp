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
#include <chrono>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "CardiacAlphaN.hpp"
#include "DataFileReader.hpp"
#include "DifferentialEvolution.hpp"
#include "Exponential.hpp"
#include "GaussianEquation.hpp"
#include "GenerateInitialGuess.hpp"
#include "GeneticAlgorithm.hpp"
#include "GeneticSwitch.hpp"
#include "LevenbergMarquardt.hpp"
#include "NelderMead.hpp"
#include "ODE3DVariant.hpp"
#include "Osborne.hpp"
#include "Parabolic.hpp"
#include "ParticleSwarm.hpp"
#include "SimulatedAnnealing.hpp"
#include "ThreeNaMarkov.hpp"
#include "Woods.hpp"
#include "UnitTest++.h"

namespace  // file scope
{
typedef std::chrono::high_resolution_clock hrclock_t;

// Needed because at present mingw has not implemented this (C++11).
// Replace to_string with std::to_string when possible
template <typename T> std::string to_string(T duration)
{
  std::stringstream converter;
  converter << duration;
  return converter.str();
}

std::string TestTime(hrclock_t::time_point t1, hrclock_t::time_point t2)
{
  using namespace std::chrono;
  auto elapsed = duration_cast<nanoseconds>(t2 - t1).count();
  if (elapsed < 10000) return {to_string(elapsed) + "\t nsec for : "};
  elapsed = duration_cast<microseconds>(t2 - t1).count();
  if (elapsed < 10000) return {to_string(elapsed) + "\t usec for : "};
  elapsed = duration_cast<milliseconds>(t2 - t1).count();
  if (elapsed < 10000) return {to_string(elapsed) + "\t msec for : "};
  elapsed = duration_cast<seconds>(t2 - t1).count();
  return {to_string(elapsed) + "\t sec  for : "};
}
}  // namespace

SUITE(ExamplesGenerateInitialGuess)
{
//
// The first set of tests here illustrate how to use the interface to the
// initial guess generation functions
//
TEST(GenerateInitialGuessSimpleInterface)
{
  const int number_of_unknowns {3};
  std::vector<double> lb {0.0, -10.0, -10.0};  // Lower bounds
  std::vector<double> ub {100.0, 10.0, 10.0};  // Upper bounds

  // Create the cost function
  std::vector<double> x;
  std::vector<double> y;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/exponential_exp_data.dat"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, x));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, y));
  CHECK_EQUAL(x.size(), y.size());
  Unfit::Examples::Exponential cost_func(x, y);

  Unfit::NelderMead nm;
  auto t1 = hrclock_t::now();  // Start time
  auto coordinates = Unfit::GenerateInitialGuess<Unfit::ParticleSwarm>(
    cost_func, lb, ub, number_of_unknowns);
  auto rc = nm.FindMin(cost_func, coordinates);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2);
  std::cout << "GenerateInitialGuess, Simple Interface" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(3.66696, coordinates[0], 1e-2);
  CHECK_CLOSE(0.0885905, coordinates[1], 1e-2);
  CHECK_CLOSE(3.06235, coordinates[2], 1e-2);
}

TEST(GenerateInitialGuessDetailedInterface)
{
  const int number_of_unknowns {3};
  std::vector<double> lb {0.0, -10.0, -10.0};  // Lower bounds
  std::vector<double> ub {100.0, 10.0, 10.0};  // Upper bounds
  const int number_of_trials {16};
  const int population_per_trial {16};
  const double cost_tol {1e-2};
  const double geom_tol {1e-4};

  // Create the cost function
  std::vector<double> x;
  std::vector<double> y;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/exponential_exp_data.dat"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, x));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, y));
  CHECK_EQUAL(x.size(), y.size());
  Unfit::Examples::Exponential cost_func(x, y);

  Unfit::LevenbergMarquardt lm;
  auto t1 = hrclock_t::now();  // Start time
  auto coordinates = Unfit::GenerateInitialGuess<Unfit::ParticleSwarm>(
    cost_func, lb, ub, number_of_unknowns, number_of_trials,
    population_per_trial, cost_tol, geom_tol);
  auto rc = lm.FindMin(cost_func, coordinates);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2);
  std::cout << "GenerateInitialGuess, Detailed Interface" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(3.66696, coordinates[0], 1e-2);
  CHECK_CLOSE(0.0885905, coordinates[1], 1e-2);
  CHECK_CLOSE(3.06235, coordinates[2], 1e-2);
}

TEST(GenerateInitialPopulationSimpleInterface)
{
  const int number_of_unknowns {3};
  std::vector<double> lb {0.0, -10.0, -10.0};  // Lower bounds
  std::vector<double> ub {100.0, 10.0, 10.0};  // Upper bounds
  std::vector<double> coordinates(number_of_unknowns, 0.0);

  // Create the cost function
  std::vector<double> x;
  std::vector<double> y;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/exponential_exp_data.dat"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, x));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, y));
  CHECK_EQUAL(x.size(), y.size());
  Unfit::Examples::Exponential cost_func(x, y);

  Unfit::DifferentialEvolution de;
  auto t1 = hrclock_t::now();  // Start time
  de.SetPopulation(Unfit::GenerateInitialPopulation<Unfit::ParticleSwarm>(
    cost_func, lb, ub, number_of_unknowns));
  auto rc = de.FindMin(cost_func, coordinates);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2);
  std::cout << "GenerateInitialPopulation, Simple Interface" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(3.66696, coordinates[0], 1e-2);
  CHECK_CLOSE(0.0885905, coordinates[1], 1e-2);
  CHECK_CLOSE(3.06235, coordinates[2], 1e-2);
}

TEST(GenerateInitialPopulationDetailedInterface)
{
  const int number_of_unknowns {3};
  std::vector<double> lb {0.0, -10.0, -10.0};  // Lower bounds
  std::vector<double> ub {100.0, 10.0, 10.0};  // Upper bounds
  const int number_of_trials {16};
  const int population_per_trial {16};
  const double cost_tol {1e-2};
  const double geom_tol {1e-4};
  std::vector<double> coordinates(number_of_unknowns, 0.0);

  // Create the cost function
  std::vector<double> x;
  std::vector<double> y;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/exponential_exp_data.dat"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, x));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, y));
  CHECK_EQUAL(x.size(), y.size());
  Unfit::Examples::Exponential cost_func(x, y);

  Unfit::DifferentialEvolution de;
  auto t1 = hrclock_t::now();  // Start time
  de.SetPopulation(Unfit::GenerateInitialPopulation<Unfit::ParticleSwarm>(
    cost_func, lb, ub, number_of_unknowns, number_of_trials,
    population_per_trial, cost_tol, geom_tol));
  auto rc = de.FindMin(cost_func, coordinates);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2);
  std::cout << "GenerateInitialPopulation, Detailed Interface" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(3.66696, coordinates[0], 1e-2);
  CHECK_CLOSE(0.0885905, coordinates[1], 1e-2);
  CHECK_CLOSE(3.06235, coordinates[2], 1e-2);
}

//
// For Differential Evolution, there are three examples that cannot find the
// correct answer with the default parameters. Here we use the initial guess
// framework to generate a good starting population then run Differential
// Evolution. All three examples then get the correct answer with the default
// parameters.
//
TEST(DifferentialEvolution_x2_CardiacAlphaN)  // From Unfit 1
{
  const int number_of_trials {16};
  const int population_per_trial {16};
  const int number_of_unknowns {3};
  std::vector<double> coordinates(number_of_unknowns, 0.0);
  std::vector<double> lb {-100.0, -100.0, -10.0};  // Lower bounds
  std::vector<double> ub {100.0, 100.0, 10.0};     // Upper bounds

  // Create the cost function
  std::vector<double> vm;
  std::vector<double> alpha_n;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/alphadata.csv"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, vm));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, alpha_n));
  CHECK_EQUAL(vm.size(), alpha_n.size());
  Unfit::Examples::CardiacAlphaN cost_func(vm, alpha_n);

  Unfit::DifferentialEvolution de;
  auto t1 = hrclock_t::now();  // Start time
  de.SetPopulation(
      Unfit::GenerateInitialPopulation<Unfit::DifferentialEvolution>(cost_func,
      lb, ub, number_of_unknowns, number_of_trials, population_per_trial));
  auto rc = de.FindMin(cost_func, coordinates);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "CardiacAlphaN ";
  std::cout << "(DifferentialEvolution/DifferentialEvolution)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0731685, coordinates[0], 1e-3);
  CHECK_CLOSE(-0.033983, coordinates[1], 1e-3);
  CHECK_CLOSE(3.10956, coordinates[2], 1e-3);
}

TEST(DifferentialEvolution_x2_Gaussian)
{
  const int number_of_trials {16};
  const int population_per_trial {16};
  const int number_of_unknowns {4};
  std::vector<double> coordinates(number_of_unknowns, 0.0);
  std::vector<double> lb {0.0, -100.0, 0.0, -100.0};    // Lower bounds
  std::vector<double> ub {100.0, 100.0, 100.0, 100.0};  // Upper bounds

  // Create the cost function
  std::vector<double> x;
  std::vector<double> y;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/gaussian_exp_data.dat"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, x));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, y));
  CHECK_EQUAL(x.size(), y.size());
  Unfit::Examples::GaussianEquation cost_func(x, y);

  Unfit::DifferentialEvolution de;
  auto t1 = hrclock_t::now();  // Start time
  de.SetPopulation(
      Unfit::GenerateInitialPopulation<Unfit::DifferentialEvolution>(cost_func,
      lb, ub, number_of_unknowns, number_of_trials, population_per_trial));
  auto rc = de.FindMin(cost_func, coordinates);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Gaussian ";
  std::cout << "(DifferentialEvolution/DifferentialEvolution)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(12.4144, coordinates[0], 1e-4);
  CHECK_CLOSE(22.7355, coordinates[1], 1e-4);
  CHECK_CLOSE(4.0337, fabs(coordinates[2]), 1e-4);
  CHECK_CLOSE(12.5486, coordinates[3], 1e-4);
}

// MB: This test example does run, but it takes a long time and so has
//     been omitted from the regular tests for now.
//
// TEST(DifferentialEvolution_LevenbergMarquardt_GeneticSwitch)
// {
//  const int number_of_trials {16};
//  const int population_per_trial {16};
//  const int number_of_unknowns {4};
//  std::vector<double> lb {0.0, 0.0, -10.0, -1.0};    // Lower bounds
//  std::vector<double> ub {1000.0, 10.0, 10.0, 1.0};  // Upper bounds
//
//  // Create the cost function
//  std::vector<double> t;
//  std::vector<double> z_a;
//  Unfit::DataFileReader<double> dfr;
//  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/training_data_zA.txt"));
//  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, t));
//  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, z_a));
//  auto dt = t[1] - t[0];
//  Unfit::Examples::GeneticSwitch cost_func(z_a, dt);
//
//  Unfit::LevenbergMarquardt lm;
//  auto t1 = hrclock_t::now();  // Start time
//  auto coordinates = Unfit::GenerateInitialGuess<Unfit::DifferentialEvolution>
//      (cost_func, lb, ub, number_of_unknowns, number_of_trials,
//      population_per_trial);
//  auto rc = lm.FindMin(cost_func, coordinates);
//  auto t2 = hrclock_t::now();  // End time
//  std::cout << TestTime(t1, t2) << "GeneticSwitch ";
//  std::cout << "(DifferentialEvolution/LevenbergMarquardt)" << std::endl;
//
//  // Check the result matches what we expect
//  CHECK_EQUAL(0, rc);
//  CHECK_CLOSE(216.433109901728, coordinates[0], 1e-2);
//  CHECK_CLOSE(2.00476590272927, coordinates[1], 1e-3);
//  CHECK_CLOSE(-2.17034564257086, coordinates[2], 1e-3);
//  CHECK_CLOSE(-0.226712906915202, coordinates[3], 1e-3);
// }

//
// For Particle Swarm, there are nine examples that cannot find the
// correct answer with the default parameters. Here we use the initial guess
// framework to generate a good starting population with Particle Swarm. For
// the second round of optimization three different strategies have been
// adoppted:
//   1. Another round of Particle Swarm
//   2. Take the best starting location and then run Levenberg Marquardt
//   3. Another round of Particle Swarm followed by Levenberg Marquardt
// All nine of the examples then get the correct answer with the default
// parameters.
//
TEST(ParticleSwarm_x2_CardiacAlphaN)  // From Unfit 1
{
  const int number_of_trials {16};
  const int population_per_trial {16};
  const int number_of_unknowns {3};
  std::vector<double> coordinates(number_of_unknowns, 0.0);
  std::vector<double> lb {-100.0, -100.0, -10.0};  // Lower bounds
  std::vector<double> ub {100.0, 100.0, 10.0};     // Upper bounds

  // Create the cost function
  std::vector<double> vm;
  std::vector<double> alpha_n;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/alphadata.csv"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, vm));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, alpha_n));
  CHECK_EQUAL(vm.size(), alpha_n.size());
  Unfit::Examples::CardiacAlphaN cost_func(vm, alpha_n);

  Unfit::ParticleSwarm ps;
  auto t1 = hrclock_t::now();  // Start time
  ps.SetPopulation(Unfit::GenerateInitialPopulation<Unfit::ParticleSwarm>(
      cost_func, lb, ub, number_of_unknowns, number_of_trials,
      population_per_trial));
  auto rc = ps.FindMin(cost_func, coordinates);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "CardiacAlphaN ";
  std::cout << "(ParticleSwarm/ParticleSwarm)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0731685, coordinates[0], 1e-3);
  CHECK_CLOSE(-0.033983, coordinates[1], 1e-3);
  CHECK_CLOSE(3.10956, coordinates[2], 1e-3);
}

TEST(ParticleSwarm_LevenbergMarquardt_Exponential)
{
  const int number_of_trials {16};
  const int population_per_trial {16};
  const int number_of_unknowns {3};
  std::vector<double> lb {0.0, -10.0, -10.0};  // Lower bounds
  std::vector<double> ub {100.0, 10.0, 10.0};  // Upper bounds

  // Create the cost function
  std::vector<double> x;
  std::vector<double> y;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/exponential_exp_data.dat"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, x));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, y));
  CHECK_EQUAL(x.size(), y.size());
  Unfit::Examples::Exponential cost_func(x, y);

  Unfit::LevenbergMarquardt lm;
  auto t1 = hrclock_t::now();  // Start time
  auto coordinates = Unfit::GenerateInitialGuess<Unfit::ParticleSwarm>(
      cost_func, lb, ub, number_of_unknowns, number_of_trials,
      population_per_trial);
  auto rc = lm.FindMin(cost_func, coordinates);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Exponential ";
  std::cout << "(ParticleSwarm/LevenbergMarquardt)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(3.66696, coordinates[0], 1e-2);
  CHECK_CLOSE(0.0885905, coordinates[1], 1e-2);
  CHECK_CLOSE(3.06235, coordinates[2], 1e-2);
}

TEST(ParticleSwarm_x2_Gaussian)
{
  const int number_of_trials {16};
  const int population_per_trial {16};
  const int number_of_unknowns {4};
  std::vector<double> coordinates(number_of_unknowns, 0.0);
  std::vector<double> lb {0.0, -100.0, 0.0, -100.0};    // Lower bounds
  std::vector<double> ub {100.0, 100.0, 100.0, 100.0};  // Upper bounds

  // Create the cost function
  std::vector<double> x;
  std::vector<double> y;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/gaussian_exp_data.dat"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, x));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, y));
  CHECK_EQUAL(x.size(), y.size());
  Unfit::Examples::GaussianEquation cost_func(x, y);

  Unfit::ParticleSwarm ps;
  auto t1 = hrclock_t::now();  // Start time
  ps.SetPopulation(Unfit::GenerateInitialPopulation<Unfit::ParticleSwarm>(
      cost_func, lb, ub, number_of_unknowns, number_of_trials,
      population_per_trial));
  auto rc = ps.FindMin(cost_func, coordinates);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Gaussian ";
  std::cout << "(ParticleSwarm/ParticleSwarm)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(12.4144, coordinates[0], 1e-4);
  CHECK_CLOSE(22.7355, coordinates[1], 1e-4);
  CHECK_CLOSE(4.0337, fabs(coordinates[2]), 1e-4);
  CHECK_CLOSE(12.5486, coordinates[3], 1e-4);
}

// MB: This test example does run, but it takes a long time and so has
//     been omitted from the regular tests for now.
//
// TEST(ParticleSwarm_LevenbergMarquardt_GeneticSwitch)
// {
//  const int number_of_trials {16};
//  const int population_per_trial {16};
//  const int number_of_unknowns {4};
//  std::vector<double> lb {0.0, 0.0, -10.0, -1.0};    // Lower bounds
//  std::vector<double> ub {1000.0, 10.0, 10.0, 1.0};  // Upper bounds
//
//  // Create the cost function
//  std::vector<double> t;
//  std::vector<double> z_a;
//  Unfit::DataFileReader<double> dfr;
//  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/training_data_zA.txt"));
//  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, t));
//  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, z_a));
//  auto dt = t[1] - t[0];
//  Unfit::Examples::GeneticSwitch cost_func(z_a, dt);
//
//  Unfit::LevenbergMarquardt lm;
//  auto t1 = hrclock_t::now();  // Start time
//  auto coordinates = Unfit::GenerateInitialGuess<Unfit::ParticleSwarm>(
//      cost_func, lb, ub, number_of_unknowns, number_of_trials,
//      population_per_trial);
//  auto rc = lm.FindMin(cost_func, coordinates);
//  auto t2 = hrclock_t::now();  // End time
//  std::cout << TestTime(t1, t2) << "GeneticSwitch ";
//  std::cout << "(ParticleSwarm/LevenbergMarquardt)" << std::endl;
//
//  // Check the result matches what we expect
//  CHECK_EQUAL(0, rc);
//  CHECK_CLOSE(216.433109901728, coordinates[0], 1e-2);
//  CHECK_CLOSE(2.00476590272927, coordinates[1], 1e-3);
//  CHECK_CLOSE(-2.17034564257086, coordinates[2], 1e-3);
//  CHECK_CLOSE(-0.226712906915202, coordinates[3], 1e-3);
// }

TEST(ParticleSwarm_LevenbergMarquardt_ODE3DVariant)  // From Unfit 1
{
  const int number_of_trials {16};
  const int population_per_trial {16};
  const int number_of_unknowns {3};
  std::vector<double> lb {-100.0, -100.0, -10.0};  // Lower bounds
  std::vector<double> ub {100.0, 100.0, 10.0};     // Upper bounds

  // Create the cost function
  std::vector<double> t;
  std::vector<double> x;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/test.txt"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, t));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, x));
  auto dt = t[1] - t[0];
  Unfit::Examples::ODE3DVariant cost_func(x, dt);

  Unfit::LevenbergMarquardt lm;
  auto t1 = hrclock_t::now();  // Start time
  auto coordinates = Unfit::GenerateInitialGuess<Unfit::ParticleSwarm>(
      cost_func, lb, ub, number_of_unknowns, number_of_trials,
      population_per_trial);
  auto rc = lm.FindMin(cost_func, coordinates);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "ODE3DVariant ";
  std::cout << "(ParticleSwarm/LevenbergMarquardt)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(2.21403, coordinates[0], 1e-2);
  CHECK_CLOSE(3.32104, coordinates[1], 1e-2);
  CHECK_CLOSE(0.0, coordinates[2], 1e-2);
}

TEST(ParticleSwarm_x2_LevenbergMarquardt_Osborne)  // From LEVMAR
{
  const int number_of_trials {32};
  const int population_per_trial {16};
  const int number_of_unknowns {5};
  std::vector<double> coordinates(number_of_unknowns, 0.0);
  std::vector<double> lb {-10.0, -10.0, -10.0, -10.0, -10.0};  // Lower bounds
  std::vector<double> ub {10.0, 10.0, 10.0, 10.0, 10.0};       // Upper bounds

  // Create the cost function
  std::vector<double> x {8.44e-1, 9.08e-1, 9.32e-1, 9.36e-1, 9.25e-1, 9.08e-1,
      8.81e-1, 8.5e-1, 8.18e-1, 7.84e-1, 7.51e-1, 7.18e-1, 6.85e-1, 6.58e-1,
      6.28e-1, 6.03e-1, 5.8e-1, 5.58e-1, 5.38e-1, 5.22e-1, 5.06e-1, 4.9e-1,
      4.78e-1, 4.67e-1, 4.57e-1, 4.48e-1, 4.38e-1, 4.31e-1, 4.24e-1, 4.2e-1,
      4.14e-1, 4.11e-1, 4.06e-1};
  Unfit::Examples::Osborne cost_func(x);

  Unfit::ParticleSwarm ps;
  auto t1 = hrclock_t::now();  // Start time
  ps.SetPopulation(Unfit::GenerateInitialPopulation<Unfit::ParticleSwarm>(
      cost_func, lb, ub, number_of_unknowns, number_of_trials,
      population_per_trial));
  auto rc = ps.FindMin(cost_func, coordinates);
  CHECK_EQUAL(0, rc);

  Unfit::LevenbergMarquardt lm;
  rc = lm.FindMin(cost_func, coordinates);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Osborne ";
  std::cout << "(ParticleSwarm_x2/LevenbergMarquardt)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.3754, coordinates[0], 1e-3);
  // Note there is the possibility of [1] & [3] swapping with [2] & [4] and the
  // solution is still valid as we are fitting two exponentials
  if (coordinates[1] > coordinates[2]) {
    CHECK_CLOSE(1.9358, coordinates[1], 1e-2);
    CHECK_CLOSE(-1.4647, coordinates[2], 1e-2);
    CHECK_CLOSE(0.0129, coordinates[3], 1e-4);
    CHECK_CLOSE(0.0221, coordinates[4], 1e-4);
  }
  else {
    CHECK_CLOSE(1.9358, coordinates[2], 1e-2);
    CHECK_CLOSE(-1.4647, coordinates[1], 1e-2);
    CHECK_CLOSE(0.0129, coordinates[4], 1e-4);
    CHECK_CLOSE(0.0221, coordinates[3], 1e-4);
  }
}

TEST(ParticleSwarm_x2_Parabolic)
{
  const int number_of_trials {16};
  const int population_per_trial {16};
  const int number_of_unknowns {3};
  std::vector<double> coordinates(number_of_unknowns, 0.0);
  std::vector<double> lb {-100.0, -100.0, -100.0};  // Lower bounds
  std::vector<double> ub {100.0, 100.0, 100.0};     // Upper bounds

  // Create the cost function
  std::vector<double> x;
  std::vector<double> y;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/parabolic_exp_data.dat"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, x));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, y));
  CHECK_EQUAL(x.size(), y.size());
  Unfit::Examples::Parabolic cost_func(x, y);

  Unfit::ParticleSwarm ps;
  auto t1 = hrclock_t::now();  // Start time
  ps.SetPopulation(Unfit::GenerateInitialPopulation<Unfit::ParticleSwarm>(
      cost_func, lb, ub, number_of_unknowns, number_of_trials,
      population_per_trial));
  auto rc = ps.FindMin(cost_func, coordinates);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Parabolic ";
  std::cout << "(ParticleSwarm/ParticleSwarm)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(-21.1144, coordinates[0], 1e-2);
  CHECK_CLOSE(61.0351, coordinates[1], 1e-2);
  CHECK_CLOSE(-9.61397, coordinates[2], 1e-2);
}

TEST(ParticleSwarm_LevenbergMarquardt_ThreeNaMarkov)  // From Unfit 1
{
  const int number_of_trials {16};
  const int population_per_trial {16};
  const int number_of_unknowns {4};
  std::vector<double> lb {0.0, 0.0, 0.0, 0.0};      // Lower bounds
  std::vector<double> ub {10.0, 10.0, 10.0, 10.0};  // Upper bounds

  // Create the cost function
  std::vector<double> t;
  std::vector<double> po;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/ThreeNaMarkov.csv"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, t));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, po));
  auto dt = t[1] - t[0];
  Unfit::Examples::ThreeNaMarkov cost_func(po, dt);

  Unfit::LevenbergMarquardt lm;
  auto t1 = hrclock_t::now();  // Start time
  auto coordinates = Unfit::GenerateInitialGuess<Unfit::ParticleSwarm>(
      cost_func, lb, ub, number_of_unknowns, number_of_trials,
      population_per_trial);
  auto rc = lm.FindMin(cost_func, coordinates);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "ThreeNaMarkov ";
  std::cout << "(ParticleSwarm/LevenbergMarquardt)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(0.0494253, coordinates[0], 1e-4);
  CHECK_CLOSE(0.0283517, coordinates[1], 1e-4);
  CHECK_CLOSE(0.321102, coordinates[2], 1e-4);
  CHECK_CLOSE(0.00174363, coordinates[3], 1e-4);
}

TEST(ParticleSwarm_LevenbergMarquardt_Woods)  // From LEVMAR
{
  const int number_of_trials {16};
  const int population_per_trial {16};
  const int number_of_unknowns {4};
  std::vector<double> lb {-10.0, -10.0, -10.0, -10.0};  // Lower bounds
  std::vector<double> ub {10.0, 10.0, 10.0, 10.0};      // Upper bounds

  // Create the cost function
  Unfit::Examples::Woods cost_func;

  Unfit::LevenbergMarquardt lm;
  auto t1 = hrclock_t::now();  // Start time
  auto coordinates = Unfit::GenerateInitialGuess<Unfit::ParticleSwarm>(
      cost_func, lb, ub, number_of_unknowns, number_of_trials,
      population_per_trial);
  auto rc = lm.FindMin(cost_func, coordinates);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Woods ";
  std::cout << "(ParticleSwarm/LevenbergMarquardt)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(1, coordinates[0], 1e-3);
  CHECK_CLOSE(1, coordinates[1], 1e-3);
  CHECK_CLOSE(1, coordinates[2], 1e-3);
  CHECK_CLOSE(1, coordinates[3], 1e-3);
}

//
// These last examples show other possible combinations in a multi-level
// optimization scheme
//
TEST(NelderMead_DifferentialEvolution_Exponential)  // From Unfit 1
{
  const int number_of_trials {16};
  const int population_per_trial {16};
  const int number_of_unknowns {3};
  std::vector<double> lb {0.0, -10.0, -10.0};  // Lower bounds
  std::vector<double> ub {100.0, 10.0, 10.0};  // Upper bounds
  std::vector<double> coordinates(number_of_unknowns, 0.0);

  // Create the cost function
  std::vector<double> x;
  std::vector<double> y;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/exponential_exp_data.dat"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, x));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, y));
  CHECK_EQUAL(x.size(), y.size());
  Unfit::Examples::Exponential cost_func(x, y);

  Unfit::DifferentialEvolution de;
  auto t1 = hrclock_t::now();  // Start time
  de.SetPopulation(Unfit::GenerateInitialPopulation<Unfit::NelderMead>(
      cost_func, lb, ub, number_of_unknowns, number_of_trials,
      population_per_trial));
  auto rc = de.FindMin(cost_func, coordinates);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Exponential ";
  std::cout << "(NelderMead/DifferentialEvolution)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(3.66696, coordinates[0], 1e-2);
  CHECK_CLOSE(0.0885905, coordinates[1], 1e-2);
  CHECK_CLOSE(3.06235, coordinates[2], 1e-2);
}

TEST(LevenbergMarquardt_DifferentialEvolution_Exponential)  // From Unfit 1
{
  const int number_of_trials {16};
  const int population_per_trial {16};
  const int number_of_unknowns {3};
  std::vector<double> lb {0.0, -10.0, -10.0};  // Lower bounds
  std::vector<double> ub {100.0, 10.0, 10.0};  // Upper bounds
  std::vector<double> coordinates(number_of_unknowns, 0.0);

  // Create the cost function
  std::vector<double> x;
  std::vector<double> y;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/exponential_exp_data.dat"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, x));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, y));
  CHECK_EQUAL(x.size(), y.size());
  Unfit::Examples::Exponential cost_func(x, y);

  Unfit::DifferentialEvolution de;
  auto t1 = hrclock_t::now();  // Start time
  de.SetPopulation(Unfit::GenerateInitialPopulation<Unfit::LevenbergMarquardt>(
      cost_func, lb, ub, number_of_unknowns, number_of_trials,
      population_per_trial));
  auto rc = de.FindMin(cost_func, coordinates);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "Exponential ";
  std::cout << "(LevenbergMarquardt/DifferentialEvolution)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(3.66696, coordinates[0], 1e-2);
  CHECK_CLOSE(0.0885905, coordinates[1], 1e-2);
  CHECK_CLOSE(3.06235, coordinates[2], 1e-2);
}

TEST(GeneticAlgorithm_LevenbergMarquardt_ODE3DVariant)  // From Unfit 1
{
  const int number_of_trials {16};
  const int population_per_trial {16};
  const int number_of_unknowns {3};
  std::vector<double> lb {-100.0, -100.0, -10.0};  // Lower bounds
  std::vector<double> ub {100.0, 100.0, 10.0};     // Upper bounds

  // Create the cost function
  std::vector<double> t;
  std::vector<double> x;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/test.txt"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, t));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, x));
  auto dt = t[1] - t[0];
  Unfit::Examples::ODE3DVariant cost_func(x, dt);

  Unfit::LevenbergMarquardt lm;
  auto t1 = hrclock_t::now();  // Start time
  auto coordinates = Unfit::GenerateInitialGuess<Unfit::GeneticAlgorithm>(
      cost_func, lb, ub, number_of_unknowns, number_of_trials,
      population_per_trial);
  auto rc = lm.FindMin(cost_func, coordinates);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "ODE3DVariant ";
  std::cout << "(GeneticAlgorithm/LevenbergMarquardt)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(2.21403, coordinates[0], 1e-2);
  CHECK_CLOSE(3.32104, coordinates[1], 1e-2);
  CHECK_CLOSE(0.0, coordinates[2], 1e-2);
}

TEST(DifferentialEvolution_SimulatedAnnealing_ODE3DVariant)  // From Unfit 1
{
  const int number_of_trials {16};
  const int population_per_trial {16};
  const int number_of_unknowns {3};
  std::vector<double> lb {-100.0, -100.0, -10.0};  // Lower bounds
  std::vector<double> ub {100.0, 100.0, 10.0};     // Upper bounds

  // Create the cost function
  std::vector<double> t;
  std::vector<double> x;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/test.txt"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, t));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, x));
  auto dt = t[1] - t[0];
  Unfit::Examples::ODE3DVariant cost_func(x, dt);

  Unfit::SimulatedAnnealing sa;
  sa.bounds.SetBounds(lb, ub);
  sa.options.SetAddInitialToPopulation(true);
  sa.options.SetMaxFunctionEvaluations(200000);
  auto t1 = hrclock_t::now();  // Start time
  auto coordinates = Unfit::GenerateInitialGuess<Unfit::DifferentialEvolution>(
      cost_func, lb, ub, number_of_unknowns, number_of_trials,
      population_per_trial);
  auto rc = sa.FindMin(cost_func, coordinates);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "ODE3DVariant ";
  std::cout << "(DifferentialEvolution/SimulatedAnnealing)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(2.21403, coordinates[0], 1e-2);
  CHECK_CLOSE(3.32104, coordinates[1], 1e-2);
  CHECK_CLOSE(0.0, coordinates[2], 1e-2);
}

TEST(SimulatedAnnealing_DifferentialEvolution_ODE3DVariant)  // From Unfit 1
{
  const int number_of_unknowns {3};
  std::vector<double> lb {-100.0, -100.0, -10.0};  // Lower bounds
  std::vector<double> ub {100.0, 100.0, 10.0};     // Upper bounds
  std::vector<double> coordinates(number_of_unknowns, 0.0);

  // Create the cost function
  std::vector<double> t;
  std::vector<double> x;
  Unfit::DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("examples/data/test.txt"));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(0, t));
  CHECK_EQUAL(0u, dfr.RetrieveColumn(1, x));
  auto dt = t[1] - t[0];
  Unfit::Examples::ODE3DVariant cost_func(x, dt);

  Unfit::DifferentialEvolution de;
  auto t1 = hrclock_t::now();  // Start time
  de.SetPopulation(Unfit::GenerateInitialPopulation<Unfit::SimulatedAnnealing>(
    cost_func, lb, ub, number_of_unknowns));
  auto rc = de.FindMin(cost_func, coordinates);
  auto t2 = hrclock_t::now();  // End time
  std::cout << TestTime(t1, t2) << "ODE3DVariant ";
  std::cout << "(SimulatedAnnealing/DifferentialEvolution)" << std::endl;

  // Check the result matches what we expect
  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(2.21403, coordinates[0], 1e-2);
  CHECK_CLOSE(3.32104, coordinates[1], 1e-2);
  CHECK_CLOSE(0.0, coordinates[2], 1e-2);
}
}  // suite ExamplesMultiLevelParticleSwarm
