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
#ifndef UNFIT_HPP_
#define UNFIT_HPP_

#include "DataFileReader.hpp"
#include "DifferentialEvolution.hpp"
#include "GeneticAlgorithm.hpp"
#include "LevenbergMarquardt.hpp"
#include "NelderMead.hpp"
#include "ParticleSwarm.hpp"
#include "SimulatedAnnealing.hpp"

/**
 * \mainpage Welcome to Unfit 3.0: Data fitting and optimization software
 *
 * Does Matlab's fminsearch solve your problems but run too slow? Does Scilab
 * handle your problems but runs even slower? Unfit may be the solution you are
 * looking for. Unfit provides access to six popular optimization algorithms
 * through a common and highly customizable interface. The algorithms are:
 *
 *   - Differential Evolution (global minimizer, non-derivative)
 *   - Genetic Algorithm (global minimizer, non-derivative)
 *   - Levenberg Marquardt (local minimizer, derivative based)
 *   - Nelder Mead Simplex (local minimizer, non-derivative)
 *   - Particle Swarm (global minimizer, non-derivative)
 *   - Simulated Annealing (global minimizer, non-derivative)
 *
 * Want to try, for example, a Hybrid Genetic Algorithm approach? No problem.
 * Just pass the solution from your Genetic Algorithm to one of the local
 * minimizers and you are done. All it takes is two extra lines of code.
 *
 * To solve your own optimization problems, all you need to do is to write your
 * own cost function. It can be very simple (see our Parabolic example) or it
 * can be as complex as you like.
 *
 * We have designed Unfit with the following key features:
 *
 *   - Speed
 *   - Bounds (box) constraints
 *   - Handles discontinuous functions
 *   - Pure C++ implementation (not a wrapper for C)
 *   - Modern C++ style (uses C++11)
 *   - Highly extensible cost functions
 *   - Flexibility to set your own options or use our defaults
 *   - Well tested (100% line and function coverage)
 *
 * Check out our web page for a tutorial on how to write your own cost function,
 * and (if interested) our guide on how to write your own optimizer using the
 * Unfit framework.
 *
 * Contact:
 *   - Dr Martin Buist martin.buist _at_ nus.edu.sg
 *   - Dr Alberto Corrias alberto _at_ nus.edu.sg
 */

#endif
