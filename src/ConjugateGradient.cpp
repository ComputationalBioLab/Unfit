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
#include <vector>
#include "ConjugateGradient.hpp"
#include "Matrix.hpp"

namespace Unfit
{
int ConjugateGradient(std::vector<double> &x, const Matrix &a,
    const std::vector<double> &b, const double epsilon)
{
  // Check to see if a solution is feasible based on the inputs provided
  if (a.Size() == 0 || b.size() == 0) return 1;
  if (a.GetNumberOfColumns() % b.size() != 0) return 1;
  // Check if vector_x is the right size or initialise it
  if (x.size() != b.size()) x.assign(b.size(), 0.0);
  auto residual = AddSubtractVectors(b, InnerProduct(a, x), false);
  auto d = residual;
  auto delta_new = InnerProduct(residual, residual);
  auto delta_zero = delta_new;
  // Check to see if we have convergence before we start iterating
  if (delta_new < epsilon*epsilon*epsilon) return 0;

  for (auto i = 1u; i <= b.size(); ++i) {
    auto q = InnerProduct(a, d);
    auto alpha = delta_new/InnerProduct(d, q);
    x = AddScaleVectors(x, d, alpha);
    if (i % 50 == 0) {
      residual = AddSubtractVectors(b, InnerProduct(a, x), false);
    } else {
      residual = AddScaleVectors(residual, q, -alpha);
    }
    auto delta_old = delta_new;
    delta_new = InnerProduct(residual, residual);
    auto beta = delta_new/delta_old;
    d = AddScaleVectors(residual, d, beta);
    if (delta_new < epsilon*epsilon*delta_zero) return 0;
  }
  return 2;  // failure to converge at end of for loop
}

}  // namespace Unfit
