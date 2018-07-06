// Unfit: Data fitting and optimization software
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
#include <vector>
#include "Matrix.hpp"
#include "UnitTest++.h"

static const double tol {1.0e-8};

namespace Unfit
{
namespace UnitTests
{
SUITE(UnitTestMatrix)
{
TEST(Matrix_Constructor1)
{
  Matrix matrix_a;
  CHECK_EQUAL(0u, matrix_a.GetNumberOfColumns());
  CHECK_EQUAL(0u, matrix_a.GetNumberOfRows());
  CHECK_EQUAL(0u, matrix_a.Size());
}

TEST(Matrix_Constructor2)
{
  std::vector<double> dummy;
  dummy.assign(20u, 0);
  for (unsigned i = 0; i < 20; ++i)    dummy[i] = i*22+i;
  Matrix matrix3(5, 4, dummy);
  CHECK_EQUAL(4u, matrix3.GetNumberOfColumns());
  CHECK_EQUAL(5u, matrix3.GetNumberOfRows());
  for (unsigned i = 0; i < 4*5; ++i) {
    CHECK_CLOSE(dummy[i], matrix3.values_[i], tol);
  }
}

TEST(Matrix_Constructor3)
{  // error msg WILL come out
  std::vector<double> dummy;
  dummy.assign(20u, 0);
  for (unsigned i = 0; i < 20; ++i)    dummy[i] = i*22+i;
  Matrix matrix4(4, 4, dummy);
  CHECK_EQUAL(0u, matrix4.GetNumberOfColumns());
  CHECK_EQUAL(0u, matrix4.GetNumberOfRows());
  CHECK_EQUAL(0u, matrix4.Size());
}

TEST(Matrix_AssignValue)
{
  Matrix matrix_test;
  matrix_test.Assign(2, 3, 0);
  CHECK_EQUAL(2u, matrix_test.GetNumberOfRows());
  CHECK_EQUAL(6u, matrix_test.Size());
}

TEST(Matrix_AssignVector)
{
  Matrix matrix_test;
  std::vector<double> vector_a {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  matrix_test.Assign(2, 3, vector_a);
  CHECK_EQUAL(2u, matrix_test.GetNumberOfRows());
  CHECK_EQUAL(6u, matrix_test.Size());
}

TEST(Matrix_AssignIncompatibleVector)
{
  Matrix matrix_test;
  std::vector<double> vector_a {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  matrix_test.Assign(2, 2, vector_a);
  CHECK(matrix_test.Empty());
}

TEST(Matrix_GetValue1)
{
  Matrix matrix_a;
  CHECK_CLOSE(0.0, matrix_a.GetValue(1, 1), 1e-8);
}

TEST(Matrix_GetValue2)
{
  Matrix matrix_a {2, 2, {1, 2,
                          3, 4}};
  CHECK_CLOSE(4.0, matrix_a.GetValue(1, 1), 1e-8);
}

TEST(Matrix_GetValue3RowsExceed)
{
  Matrix matrix_a {2, 2, {1, 2,
                          3, 4}};
  CHECK_CLOSE(0.0, matrix_a.GetValue(3, 1), 1e-8);
}

TEST(Matrix_GetValue4ColumnsExceed)
{
  Matrix matrix_a {2, 2, {1, 2,
                          3, 4}};
  CHECK_CLOSE(0.0, matrix_a.GetValue(1, 3), 1e-8);
}

TEST(Matrix_GetNumberOfRowsNullMatrix)
{
  Matrix matrix_a;
  CHECK_EQUAL(0u, matrix_a.GetNumberOfRows());
}

TEST(Matrix_GetNumberOfColumnsNullMatrix)
{
  Matrix matrix_a;
  CHECK_EQUAL(0u, matrix_a.GetNumberOfColumns());
}

TEST(Matrix_GetNumberOfRows1x1)
{
  Matrix matrix_a {1, 1, 0};
  CHECK_EQUAL(1u, matrix_a.GetNumberOfRows());
}

TEST(Matrix_GetNumberOfColumns1x1)
{
  Matrix matrix_a {1, 1, 0};
  CHECK_EQUAL(1u, matrix_a.GetNumberOfColumns());
}

TEST(Matrix_GetNumberOfRows2x3)
{
  Matrix matrix_a {2, 3, {1, 2,
                          3, 4,
                          5, 6}};
  CHECK_EQUAL(2u, matrix_a.GetNumberOfRows());
}

TEST(Matrix_GetNumberOfColumns2x3)
{
  Matrix matrix_a {2, 3, {1, 2,
                          3, 4,
                          5, 6}};
  CHECK_EQUAL(3u, matrix_a.GetNumberOfColumns());
}

TEST(Matrix_SetValue)
{
  Matrix matrix_a(2, 2, 1.0);
  CHECK_CLOSE(1.0, matrix_a.GetValue(1, 1), 1e-8);
  matrix_a.SetValue(1, 1, 3.0);
  CHECK_CLOSE(3.0, matrix_a.GetValue(1, 1), 1e-8);
  CHECK_CLOSE(1.0, matrix_a.GetValue(0, 0), 1e-8);
  CHECK_CLOSE(1.0, matrix_a.GetValue(0, 1), 1e-8);
  CHECK_CLOSE(1.0, matrix_a.GetValue(1, 0), 1e-8);
}

TEST(Matrix_SetValueOutOfBounds)
{
  Matrix matrix_a(2, 2, 1.0);
  CHECK_CLOSE(1.0, matrix_a.GetValue(1, 1), 1e-8);
  matrix_a.SetValue(2, 1, 3.0);
  CHECK_CLOSE(1.0, matrix_a.GetValue(1, 1), 1e-8);
  CHECK_CLOSE(1.0, matrix_a.GetValue(0, 0), 1e-8);
  CHECK_CLOSE(1.0, matrix_a.GetValue(0, 1), 1e-8);
  CHECK_CLOSE(1.0, matrix_a.GetValue(1, 0), 1e-8);
  matrix_a.SetValue(1, 2, 3.0);
  CHECK_CLOSE(1.0, matrix_a.GetValue(1, 1), 1e-8);
  CHECK_CLOSE(1.0, matrix_a.GetValue(0, 0), 1e-8);
  CHECK_CLOSE(1.0, matrix_a.GetValue(0, 1), 1e-8);
  CHECK_CLOSE(1.0, matrix_a.GetValue(1, 0), 1e-8);
  matrix_a.SetValue(2, 2, 3.0);
  CHECK_CLOSE(1.0, matrix_a.GetValue(1, 1), 1e-8);
  CHECK_CLOSE(1.0, matrix_a.GetValue(0, 0), 1e-8);
  CHECK_CLOSE(1.0, matrix_a.GetValue(0, 1), 1e-8);
  CHECK_CLOSE(1.0, matrix_a.GetValue(1, 0), 1e-8);
}

TEST(Matrix_SingleMatrix)
{
  Matrix matrix_a {1, 1, 1};
  CHECK_CLOSE(1.0, matrix_a.GetValue(0, 0), 1e-8);
}

TEST(Matrix_TwoxTwoMatrix)
{
  Matrix matrix_a {2, 2, {1, 2,
                          3, 4}};
  CHECK_CLOSE(3.0, matrix_a.GetValue(1, 0), 1e-8);
}

TEST(Matrix_ThreexTwoMatrix)
{
  Matrix matrix_a {3, 2, {1, 2, 3,
                          4, 5, 6}};
  CHECK_CLOSE(4.0, matrix_a.GetValue(1, 1), 1e-8);
}

TEST(Matrix_TwoxFiveMatrix)
{
  Matrix matrix_a {2, 5, {1, 2,
                          3, 4,
                          5, 6,
                          7, 8,
                          9, 10}};
  CHECK_CLOSE(9.0, matrix_a.GetValue(1, 3), 1e-8);
}

TEST(Matrix_MatrixVectorProductPositiveElements)
{
  Matrix matrix_a {2, 2, {1, 0, 0, 1}};
  std::vector<double> vector_x {1.0, 2.0};

  CHECK_CLOSE(1, InnerProduct(matrix_a, vector_x)[0], tol);
  CHECK_CLOSE(2, InnerProduct(matrix_a, vector_x)[1], tol);
}

TEST(Matrix_MatrixVectorProductNullMatrixA)
{
  Matrix matrix_a {};
  std::vector<double> vector_x {1.0, 2.0};

  CHECK_EQUAL(0u, InnerProduct(matrix_a, vector_x).size());
}

TEST(Matrix_MatrixVectorProductNullVectorX)
{
  Matrix matrix_a {1, 2, {1, 2}};
  std::vector<double> vector_x {};

  CHECK_EQUAL(0u, InnerProduct(matrix_a, vector_x).size());
}

TEST(Matrix_MatrixVectorProductNullMatrixAndVector)
{
  Matrix matrix_a {};
  std::vector<double> vector_x {};

  CHECK_EQUAL(0u, InnerProduct(matrix_a, vector_x).size());
}

TEST(Matrix_MatrixVectorProductIncompatibleSize)
{
  Matrix matrix_a {1, 3, {1, 2, 3}};
  std::vector<double> vector_x {1.0, 2.0};

  CHECK_EQUAL(0u, InnerProduct(matrix_a, vector_x).size());  // 1
}

TEST(Matrix_MatrixVectorProductIncompatibleSizeTwo)
{
  Matrix matrix_a {1, 3, {1, 2, 3}};
  std::vector<double> vector_x {1.0, 2.0, 3.0, 4.0, 5.0};

  CHECK_EQUAL(0u, InnerProduct(matrix_a, vector_x).size());  // 1
}

TEST(Matrix_MatrixVectorProductNotSquareMatrixFail)
{
  Matrix matrix_a {2, 3, {1, 2, 3,
                          4, 5, 6}};
  std::vector<double> vector_x {1.0, 2.0};

  CHECK_EQUAL(0u, InnerProduct(matrix_a, vector_x).size());  // 2
}

TEST(Matrix_MatrixVectorProductNotSquareMatrixPass)
{
  Matrix matrix_a {3, 2, {1, 2,
                          3, 4,
                          5, 6}};
  std::vector<double> vector_x {1.0, 2.0};

  CHECK_CLOSE(5u, InnerProduct(matrix_a, vector_x)[0], tol);
  CHECK_CLOSE(11u, InnerProduct(matrix_a, vector_x)[1], tol);
  CHECK_CLOSE(17u, InnerProduct(matrix_a, vector_x)[2], tol);
  CHECK_EQUAL(3u, InnerProduct(matrix_a, vector_x).size());
}

TEST(Matrix_MatrixVectorProductNotSquareMatrixFailTwo)
{
  Matrix matrix_a {3, 2, {1, 2,
                          3, 4,
                          5, 6}};
  std::vector<double> vector_x {1.0, 2.0, 3.0};
  CHECK_EQUAL(0u, InnerProduct(matrix_a, vector_x).size());  // 3
}

TEST(Matrix_MatrixVectorProductCalculationCheckOne)
{
  Matrix matrix_a {8, 8, {1, 2, 3, 4, 5, 6, 7, 8,
                          1, 2, 3, 4, 5, 6, 7, 8,
                          1, 2, 3, 4, 5, 6, 7, 8,
                          1, 2, 3, 4, 5, 6, 7, 8,
                          1, 2, 3, 4, 5, 6, 7, 8,
                          1, 2, 3, 4, 5, 6, 7, 8,
                          1, 2, 3, 4, 5, 6, 7, 8,
                          1, 2, 3, 4, 5, 6, 7, 8}};
  std::vector<double> vector_x {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
  std::vector<double> result = InnerProduct(matrix_a, vector_x);
  for (unsigned i = 0; i < 8; ++i) CHECK_CLOSE(204, result[i], tol);
}

TEST(Matrix_InnerProductPositiveElements)
{
  std::vector<double> vector_a {1.0, 2.0, 3.0, 4.0};
  std::vector<double> vector_b {1.0, 2.0, 3.0, 4.0};
  CHECK_CLOSE(30.0, InnerProduct(vector_a, vector_b), tol);
}

TEST(Matrix_InnerProductZeroProduct)
{
  std::vector<double> vector_a {2.0, 1.0, 1.0};
  std::vector<double> vector_b {-1.0, 1.0, 1.0};
  CHECK_CLOSE(0.0, InnerProduct(vector_a, vector_b), tol);
}

TEST(Matrix_InnerProductNullVectorA)
{
  std::vector<double> vector_a {};
  std::vector<double> vector_b {1.0, 2.0, 3.0, 4.0};
  CHECK(std::isnan(InnerProduct(vector_a, vector_b)));
}

TEST(Matrix_InnerProductNullVectorB)
{
  std::vector<double> vector_a {1.0, 2.0};
  std::vector<double> vector_b {};
  CHECK(std::isnan(InnerProduct(vector_a, vector_b)));
}

TEST(Matrix_InnerProductImproperSizes)
{
  std::vector<double> vector_a {1.0, 2.0, 3.0};
  std::vector<double> vector_b {1.0, 2.0, 3.0, 4.0};
  CHECK(std::isnan(InnerProduct(vector_a, vector_b)));
}

TEST(Matrix_ScaleVector)
{
  double beta = 2.0;
  std::vector<double> vector_a {1.0, 2.0, 3.0, 4.0};
  CHECK_CLOSE(2.0, Scale(vector_a, beta)[0], tol);
  CHECK_CLOSE(4.0, Scale(vector_a, beta)[1], tol);
  CHECK_CLOSE(6.0, Scale(vector_a, beta)[2], tol);
  CHECK_CLOSE(8.0, Scale(vector_a, beta)[3], tol);
}

TEST(Matrix_AddVectorPositiveElements)
{
  std::vector<double> vector_a {1.0, 2.0, 3.0};
  std::vector<double> vector_b {4.0, 5.0, 6.0};
  CHECK_CLOSE(5, AddSubtractVectors(vector_a, vector_b)[0], tol);
  CHECK_CLOSE(7, AddSubtractVectors(vector_a, vector_b)[1], tol);
  CHECK_CLOSE(9, AddSubtractVectors(vector_a, vector_b)[2], tol);
  CHECK_EQUAL(3u, AddSubtractVectors(vector_a, vector_b).size());
}

TEST(Matrix_AddVectorNullVectorA)
{
  std::vector<double> vector_a {};
  std::vector<double> vector_b {1.0, 2.0, 3.0};
  CHECK_EQUAL(0u, AddSubtractVectors(vector_a, vector_b).size());
}

TEST(Matrix_AddVectorNullVectors)
{
  std::vector<double> vector_a {};
  std::vector<double> vector_b {};
  CHECK_EQUAL(0u, AddSubtractVectors(vector_a, vector_b).size());
}

TEST(Matrix_AddVectorZeroVectors)
{
  std::vector<double> vector_a {0.0, 0.0, 0.0};
  std::vector<double> vector_b {0.0, 0.0, 0.0};
  CHECK_CLOSE(0, AddSubtractVectors(vector_a, vector_b)[0], tol);
  CHECK_CLOSE(0, AddSubtractVectors(vector_a, vector_b)[1], tol);
  CHECK_CLOSE(0, AddSubtractVectors(vector_a, vector_b)[2], tol);
  CHECK_EQUAL(3u, AddSubtractVectors(vector_a, vector_b).size());
}

TEST(Matrix_SubtractVectorPositiveElements)
{
  std::vector<double> vector_a {1.0, 2.0, 3.0, 4.0};
  std::vector<double> vector_b {1.0, 2.0, 3.0, 4.0};
  // for (auto i= 0; i < 4; ++i) {

  CHECK_CLOSE(0u, AddSubtractVectors(vector_a, vector_b, false)[0], tol);
  CHECK_CLOSE(0u, AddSubtractVectors(vector_a, vector_b, false)[1], tol);
  CHECK_CLOSE(0u, AddSubtractVectors(vector_a, vector_b, false)[2], tol);
  CHECK_CLOSE(0u, AddSubtractVectors(vector_a, vector_b, false)[3], tol);
  CHECK_EQUAL(4u, AddSubtractVectors(vector_a, vector_b).size());
}

TEST(Matrix_SubtractVectorNegativeElements)
{
  std::vector<double> vector_a {1.0, 2.0, 3.0, 4.0};
  std::vector<double> vector_b {-1.0, -2.0, -3.0, -4.0};
  for (auto i = 0u; i < 4; ++i) {
    CHECK_CLOSE((i+1.0)*2.0,
        AddSubtractVectors(vector_a, vector_b, false)[i], 1e-8);
  }
  CHECK_EQUAL(4u, AddSubtractVectors(vector_a, vector_b).size());
}

TEST(Matrix_AddVectorImproperSizesVectorA)
{
  std::vector<double> vector_a {1.0, 2.0, 3.0, 4.0};
  std::vector<double> vector_b {4.0, 5.0, 6.0};
  CHECK_EQUAL(0u, AddSubtractVectors(vector_a, vector_b).size());
}

TEST(Matrix_AddVectorImproperSizesVectorB)
{
  std::vector<double> vector_a {1.0, 2.0};
  std::vector<double> vector_b {4.0, 5.0, 6.0};
  CHECK_EQUAL(0u, AddSubtractVectors(vector_a, vector_b).size());
}

TEST(Matrix_AddMuAlongDiagonal)
{
  Matrix a {3, 3, 0};
  a.AddConstantToDiagonal(1.0);
  CHECK_CLOSE(1.0, a.values_[0], 1e-8);
}

TEST(Matrix_TransposeMatrix1)
{
  Matrix a {2, 3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  Matrix b {Transpose(a)};
  CHECK_EQUAL(2u, b.GetNumberOfColumns());
  CHECK_EQUAL(3u, b.GetNumberOfRows());
  CHECK_CLOSE(1, b.values_[0], tol);
  CHECK_CLOSE(4, b.values_[1], tol);
  CHECK_CLOSE(2, b.values_[2], tol);
  CHECK_CLOSE(5, b.values_[3], tol);
  CHECK_CLOSE(3, b.values_[4], tol);
  CHECK_CLOSE(6, b.values_[5], tol);
}

TEST(Matrix_TransposeEmptyMatrix)
{
  Matrix a {0, 0, 0.0};
  Matrix b {Transpose(a)};
  CHECK_EQUAL(0u, b.GetNumberOfColumns());
  CHECK_EQUAL(0u, b.GetNumberOfRows());
  CHECK_EQUAL(0u, b.Size());
}

TEST(Matrix_FindL2NormOfVectorVectorAEmpty)
{
  std::vector<double> vector_a {};
  CHECK_CLOSE(0.0, FindL2NormOfVector(vector_a), 1e-15);
}

TEST(Matrix_MaxDiagonalMatrix1)
{
  Matrix a = {2, 2, {1.0, 2.0, 3.0, 4.0}};
  CHECK_CLOSE(4, MaxDiagonalMat(a), tol);
}

TEST(Matrix_MaxDiagonalMatrix2)
{
  Matrix a = {1, 5, {1.0, 2.0, 3.0, 4.0, 5.0}};
  CHECK_CLOSE(0, MaxDiagonalMat(a), tol);
}

TEST(Matrix_MatrixTransposeTimesMatrix)
{
  Matrix a {2, 3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  // Calculate aTa
  Matrix c = InnerProduct(a, false);
  CHECK_EQUAL(9u, c.Size());
  CHECK_CLOSE(17.0, c.values_[0], 1.0e-8);
  CHECK_CLOSE(22.0, c.values_[1], 1.0e-8);
  CHECK_CLOSE(27.0, c.values_[2], 1.0e-8);
  CHECK_CLOSE(22.0, c.values_[3], 1.0e-8);
  CHECK_CLOSE(29.0, c.values_[4], 1.0e-8);
  CHECK_CLOSE(36.0, c.values_[5], 1.0e-8);
  CHECK_CLOSE(27.0, c.values_[6], 1.0e-8);
  CHECK_CLOSE(36.0, c.values_[7], 1.0e-8);
  CHECK_CLOSE(45.0, c.values_[8], 1.0e-8);

  // Calculate aaT
  Matrix d = InnerProduct(a, true);
  CHECK_EQUAL(4u, d.Size());
  CHECK_CLOSE(14.0, d.values_[0], 1.0e-8);
  CHECK_CLOSE(32.0, d.values_[1], 1.0e-8);
  CHECK_CLOSE(32.0, d.values_[2], 1.0e-8);
  CHECK_CLOSE(77.0, d.values_[3], 1.0e-8);

  // Now try an empty matrix
  Matrix b;
  Matrix e = InnerProduct(b, true);
  CHECK(e.Empty());
}

TEST(Matrix_AddConstantToDiagonalSquareMatrix)
{
  Matrix a {3, 3, 1};
  a.AddConstantToDiagonal(1.0);
  CHECK_EQUAL(9u, a.Size());
  CHECK_CLOSE(2.0, a.values_[0], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[1], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[2], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[3], 1.0e-8);
  CHECK_CLOSE(2.0, a.values_[4], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[5], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[6], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[7], 1.0e-8);
  CHECK_CLOSE(2.0, a.values_[8], 1.0e-8);
}

TEST(Matrix_AddConstantToUpperDiagonalSquareMatrix)
{
  Matrix a {3, 3, 1};
  a.AddConstantToDiagonal(1.0, 1, true);
  CHECK_EQUAL(9u, a.Size());
  CHECK_CLOSE(1.0, a.values_[0], 1.0e-8);
  CHECK_CLOSE(2.0, a.values_[1], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[2], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[3], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[4], 1.0e-8);
  CHECK_CLOSE(2.0, a.values_[5], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[6], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[7], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[8], 1.0e-8);
}

TEST(Matrix_AddConstantToLowerDiagonalSquareMatrix)
{
  Matrix a {3, 3, 1};
  a.AddConstantToDiagonal(1.0, 1, false);
  CHECK_EQUAL(9u, a.Size());
  CHECK_CLOSE(1.0, a.values_[0], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[1], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[2], 1.0e-8);
  CHECK_CLOSE(2.0, a.values_[3], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[4], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[5], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[6], 1.0e-8);
  CHECK_CLOSE(2.0, a.values_[7], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[8], 1.0e-8);
}

TEST(Matrix_AddConstantToDiagonalNonSquareMatrix)
{
  Matrix a {4, 3, 1};
  a.AddConstantToDiagonal(1.0);
  CHECK_EQUAL(12u, a.Size());
  CHECK_CLOSE(2.0, a.values_[0], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[1], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[2], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[3], 1.0e-8);
  CHECK_CLOSE(2.0, a.values_[4], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[5], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[6], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[7], 1.0e-8);
  CHECK_CLOSE(2.0, a.values_[8], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[9], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[10], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[11], 1.0e-8);
}

TEST(Matrix_AddConstantToDiagonalNonSquareMatrix2)
{
  Matrix a {3, 4, 1};
  a.AddConstantToDiagonal(1.0);
  CHECK_EQUAL(12u, a.Size());
  CHECK_CLOSE(2.0, a.values_[0], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[1], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[2], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[3], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[4], 1.0e-8);
  CHECK_CLOSE(2.0, a.values_[5], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[6], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[7], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[8], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[9], 1.0e-8);
  CHECK_CLOSE(2.0, a.values_[10], 1.0e-8);
  CHECK_CLOSE(1.0, a.values_[11], 1.0e-8);
}

TEST(Matrix_AddScaleVectors)
{
  std::vector<double> v1(3, 2.0);
  std::vector<double> v2(3, 1.0);
  auto result = AddScaleVectors(v1, v2, 2.0);
  CHECK_EQUAL(3u, result.size());
  CHECK_CLOSE(4.0, result[0], 1.0e-8);
  CHECK_CLOSE(4.0, result[1], 1.0e-8);
  CHECK_CLOSE(4.0, result[2], 1.0e-8);
}

TEST(Matrix_AddScaleVectorsIncompatibleSizes)
{
  std::vector<double> v1(4, 2.0);
  std::vector<double> v2(3, 1.0);
  auto result = AddScaleVectors(v1, v2, 2.0);
  CHECK(result.empty());
}
}  // suite UnitTestMatrix
}  // namespace UnitTests
}  // namespace Unfit
