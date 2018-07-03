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
#ifndef UNFIT_INCLUDE_CONJUGATEGRADIENT_HPP_
#define UNFIT_INCLUDE_CONJUGATEGRADIENT_HPP_

#include <cmath>
#include <vector>
#include "Matrix.hpp"

namespace Unfit
{
/**
 * This method solves a linear system of equations by the conjugate gradient
 * method. The conjugate gradient method is an iterative method to solve for
 * a vector x in the equation : Ax = b, where A is a symmetric square matrix
 * and b is a vector.
 *
 * This algorithm is from Conjugate Gradients Without the Agonizing Pain
 * by Jonathan Shewchuk, page 50.
 *
 * Implemetation notes:
 *  - If vector_x (solution) is the correct size, whatever is passed in will
 *    be used as the initial guess
 *  - If vector_x (solution) is empty, or the wrong size, a vector of zeros will
 *    be used as the initial guess
 *  - The convergence tolerance is optional, with the default being 1e-8
 *  - For n equations, the maximum number of iterations is also n
 *
 * Intended use:
 *   auto rc = ConjugateGradient(x, a, b);
 *
 * \param    x Contains the initial guess for x (optional) on input. Contains
 *           the resulting solution vector on output.
 * \param    a A square matrix of size (n x n). The matrix should be
 *           symmetric and positive definite.
 * \param    b The right-hand-side vector.
 * \param    epsilon (optional) is the tolerance for the convergence.
 * \return   an error code, 0 means success, 1 signifies invalid input,
 *           2 signifies failure to converge.
 */
int ConjugateGradient(std::vector<double> &x, const Matrix &a,
    const std::vector<double> &b, const double epsilon = 1.0e-8);

}  // namespace Unfit

#endif
