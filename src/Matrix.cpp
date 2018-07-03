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
#include <numeric>
#include <vector>
#include "Matrix.hpp"

namespace Unfit
{
Matrix::Matrix()
  : values_ {},
    number_of_rows_ {0},
    number_of_columns_ {0}
{}

Matrix::Matrix(std::size_t num_rows, std::size_t num_cols, double the_value)
  : values_(num_rows*num_cols, the_value),
    number_of_rows_ {num_rows},
    number_of_columns_ {num_cols}
{}

Matrix::Matrix(std::size_t num_rows, std::size_t num_cols,
    const std::vector<double> &the_values)
  : values_ {0},
    number_of_rows_ {num_rows},
    number_of_columns_ {num_cols}
{
  if (the_values.size() == num_rows*num_cols) {
    values_ = the_values;
  } else {
    // should change to throw catch!!
    number_of_rows_ = 0;
    number_of_columns_ = 0;
    values_.clear();
  }
}

void Matrix::Assign(std::size_t num_rows, std::size_t num_cols,
    double the_value)
{
  number_of_rows_ = num_rows;
  number_of_columns_ = num_cols;
  values_.assign(num_rows*num_cols, the_value);
}

void Matrix::Assign(std::size_t num_rows, std::size_t num_cols,
    const std::vector<double> &the_values)
{
  if (the_values.size() == num_rows*num_cols) {
    number_of_rows_ = num_rows;
    number_of_columns_ = num_cols;
    values_ = the_values;
  } else {
    // should change to throw catch!!
    number_of_rows_ = 0;
    number_of_columns_ = 0;
    values_.clear();
  }
}

double Matrix::GetValue(std::size_t row, std::size_t col) const noexcept
{
  if (row > number_of_rows_ || col > number_of_columns_) return 0.0;
  return values_[row * number_of_columns_ + col];
}

void Matrix::SetValue(std::size_t row, std::size_t column, double value)
{
  if ((row < number_of_rows_) && (column < number_of_columns_)) {
    values_[(row*number_of_columns_)+column] = value;
  }
}

std::size_t Matrix::GetNumberOfColumns() const noexcept
{
  return number_of_columns_;
}

std::size_t Matrix::GetNumberOfRows() const noexcept
{
  return number_of_rows_;
}

std::size_t Matrix::Size() const noexcept
{
  return values_.size();
}

bool Matrix::Empty() const noexcept
{
  return values_.empty();
}

void Matrix::AddConstantToDiagonal(double value, std::size_t offset,
    bool upper_triangle)
{
  std::size_t col = 0u;
  auto row = offset;
  if (upper_triangle) {
    col = offset;
    row = 0u;
  }
  while ((col < number_of_columns_) && (row < number_of_rows_)) {
    values_[(row*number_of_columns_)+col] += value;
    ++col;
    ++row;
  }
}

std::vector<double> InnerProduct(const Matrix &m, const std::vector<double> &v)
{
  auto num_col = m.GetNumberOfColumns();
  auto num_row = m.GetNumberOfRows();
  // Check whether we have a zero length matrix or x vector
  if ((num_col*num_row == 0) || (v.size() == 0)) return {};
  // Check whether a matrix-vector product is possible from the dimensions
  if (num_col != v.size()) return {};

  std::vector<double> result(num_row, 0.0);
  auto it_m = begin(m.values_);
  for (auto &r : result) {
    for (auto x : v) r += (*it_m++) * x;
  }
  return result;
}

double InnerProduct(const std::vector<double> &v1,
    const std::vector<double> &v2)
{
  if (v1.empty()) return std::nan("");
  if (v1.size() != v2.size()) return std::nan("");
  return std::inner_product(begin(v1), end(v1), begin(v2), 0.0);
}

std::vector<double> AddSubtractVectors(const std::vector<double> &v1,
    const std::vector<double> &v2, const bool isadd)
{
  auto vector_size = v1.size();
  std::vector<double> result {};
  if (vector_size != v2.size()) return result;
  result = v1;
  if (isadd) {
    for (auto i = 0u; i < vector_size; ++i) result[i] += v2[i];
  }
  else {
    for (auto i = 0u; i < vector_size; ++i) result[i] -= v2[i];
  }
  return result;
}

std::vector<double> AddScaleVectors(const std::vector<double> &v1,
    const std::vector<double> &v2, const double scale_factor)
{
  auto vector_size = v1.size();
  if (vector_size != v2.size()) return {};
  std::vector<double> result(vector_size);
  std::transform(begin(v1), end(v1), begin(v2), begin(result),
      [&](double e1, double e2) { return e1 + scale_factor * e2; } );
  return result;
}

Matrix InnerProduct(const Matrix &m, bool transpose)
{
  auto num_col = m.GetNumberOfColumns();
  auto num_row = m.GetNumberOfRows();
  // Check our matrix has a non-zero size
  if ((num_col*num_row == 0)) return {};

  if (transpose) {
    Matrix result(num_row, num_row, 0.0);
    for (auto row = 0u; row < num_row; ++row) {
      auto row_st = row*num_col;
      for (auto col = 0u; col < num_row; ++col) {
        auto pos = row*num_row + col;
        if (col < row) {  // Copy what we have calculated already (symmetry)
          result.values_[pos] = result.values_[col*num_row + row];
        }
        else {            // or Calculate via an inner product
          auto col_st = col*num_col;
          for (auto i = 0u; i < num_col; ++i) {
            result.values_[pos] += m.values_[row_st + i]*m.values_[col_st + i];
          }
        }
      }
    }
    return result;
  }
  else {
    Matrix result(num_col, num_col, 0.0);
    for (auto row = 0u; row < num_col; ++row) {
      for (auto col = 0u; col < num_col; ++col) {
        auto pos = row*num_col + col;
        if (col < row) {  // Copy what we have calculated already (symmetry)
          result.values_[pos] = result.values_[col*num_col + row];
        }
        else {            // or Calculate via an inner product
          for (auto i = 0u; i < num_row; ++i) {
            result.values_[pos] += m.values_[row + i*num_col] *
                m.values_[col + i*num_col];
          }
        }
      }
    }
    return result;
  }
}

Matrix Transpose(Matrix &m)
{
  auto a_row = m.GetNumberOfRows();
  auto a_col = m.GetNumberOfColumns();
  // row and col is swapped
  Matrix matrix_transpose(a_col, a_row, 0);
  if (m.Empty()) return matrix_transpose;

  for (auto i = 0u; i < a_row; ++i) {
    for (auto j = 0u; j < a_col; ++j) {
      matrix_transpose.values_[j* a_row+ i] = m.values_[i*a_col + j];
    }
  }
  return matrix_transpose;
}

double FindL2NormOfVector(const std::vector<double> &v)
{
  return sqrt(std::inner_product(begin(v), end(v), begin(v), 0.0));
}

double SumOfSquares(const std::vector<double> &v)
{
  return std::inner_product(begin(v), end(v), begin(v), 0.0);
}

double MaxDiagonalMat(const Matrix &m)
{
  double the_max {m.values_[0]};
  auto length = m.GetNumberOfRows();
  if (m.Size() != length*length) return 0;
  for (auto i = 0u; i < length; ++i) {
    the_max = std::max(m.values_[i*length+i], the_max);
  }
  return the_max;
}

std::vector<double> Scale(const std::vector<double> &v,
    const double scale_factor)
{
  auto result = v;
  for (auto &i : result) i *= scale_factor;
  return result;
}

}  // namespace Unfit
