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
#include "DataFileReader.hpp"
#include "UnitTest++.h"

static const double tol {1.0e-8};

namespace Unfit
{
namespace UnitTests
{
SUITE(UnitTesetConjugateGradient)
{
TEST(ConjugateGradient_SimpleSolve)
{
  Matrix hessian_a {2, 2, {4, 1, 1, 3}};
  std::vector<double> gradient_b {-1.0, -2.0};
  std::vector<double> delta_x {100.0, 100.0};
  int rc {1};
  rc = ConjugateGradient(delta_x, hessian_a, gradient_b);

  CHECK_EQUAL(0, rc);
  CHECK_CLOSE(-1.0/11, delta_x[0], tol);
  CHECK_CLOSE(-7.0/11, delta_x[1], tol);
}

TEST(ConjugateGradient_MatrixAZeroSize)
{
  Matrix hessian_a;
  std::vector<double> gradient_b {-1.0, -2.0};
  std::vector<double> delta_x {100.0, 100.0};
  int rc {1};
  rc = ConjugateGradient(delta_x, hessian_a, gradient_b);

  CHECK_EQUAL(1, rc);
}

TEST(ConjugateGradient_MatrixAInvalidSize)
{
  Matrix hessian_a {1, 3, {4, 1, 1}};
  std::vector<double> gradient_b {-1.0, -2.0};
  std::vector<double> delta_x {100.0, 100.0};
  int rc {1};
  rc = ConjugateGradient(delta_x, hessian_a, gradient_b);

  CHECK_EQUAL(1, rc);
}

TEST(ConjugateGradient_VectorBInvalidSize)
{
  Matrix hessian_a {2, 2, {4, 1, 1, 3}};
  std::vector<double> gradient_b {-1.0, -2.0, -3.0};
  std::vector<double> delta_x {100.0, 100.0};
  int rc {1};
  rc = ConjugateGradient(delta_x, hessian_a, gradient_b);

  CHECK_EQUAL(1, rc);
}

TEST(ConjugateGradient_VectorXInvalidSize)
{
  Matrix hessian_a {2, 2, {4, 1, 1, 3}};
  std::vector<double> gradient_b {-1.0, -2.0};
  std::vector<double> delta_x {100.0, 100.0, 100.0};
  int rc {1};
  rc = ConjugateGradient(delta_x, hessian_a, gradient_b);

  CHECK_EQUAL(0, rc);
}

TEST(ConjugateGradient_AllInvalidSize)
{
  Matrix hessian_a {0, 0, {}};
  std::vector<double> gradient_b {};
  std::vector<double> delta_x {};
  int rc =ConjugateGradient(delta_x, hessian_a, gradient_b);
  CHECK_EQUAL(1, rc);
}

TEST(ConjugateGradient_a0check)
{
  Matrix hessian_a {2, 2, {0, 0, 0, 0}};
  std::vector<double> gradient_b {-1.0, -2.0};
  std::vector<double> delta_x {10.0, 10.0};
  int rc {1};

  rc = ConjugateGradient(delta_x, hessian_a, gradient_b);

  CHECK_EQUAL(2, rc);
}

TEST(ConjugateGradient_LargerMatrix)
{
  DataFileReader<double> dfr;
  CHECK_EQUAL(0u, dfr.ReadFile("./unittests/data/matrix.csv"));
  auto matrix_data = dfr.RetrieveDataRowWiseAsVector();
  Matrix a(101, 101, matrix_data);
  CHECK_EQUAL(10201u, a.Size());
  a = InnerProduct(a, true);       // Make it symmetric
  a.AddConstantToDiagonal(100.0);  // Make it more diagonally dominant
  CHECK_EQUAL(0u, dfr.ReadFile("./unittests/data/vector.csv"));
  auto b = dfr.RetrieveDataRowWiseAsVector();
  CHECK_EQUAL(101u, b.size());
  std::vector<double> x;
  CHECK_EQUAL(0, ConjugateGradient(x, a, b));
  CHECK_EQUAL(101u, x.size());
  CHECK_CLOSE(-0.00959955, x[0], tol);
  CHECK_CLOSE(-0.00269212, x[10], tol);
  CHECK_CLOSE(0.00305225, x[100], tol);
}
}  // suite UnitTestConjugateGradient
}  // namespace UnitTests
}  // namespace Unfit
