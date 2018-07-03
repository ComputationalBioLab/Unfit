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
#ifndef UNFIT_INCLUDE_MATRIX_HPP_
#define UNFIT_INCLUDE_MATRIX_HPP_

#include <vector>

namespace Unfit
{
/**
 * A class to store and operate on a 2D matrix. The matrix does not have to be
 * square. The data is stored in a 1D vector row-wise.
 */
class Matrix
{
 public:
  /**
   * Constructor: (default) Override the default constructor to set the number
   * of rows and columns to zero explicitly so we can't access any entries by
   * accident.
   */
  Matrix();

  /**
   * Constructor: Creates a matrix with the specified dimensions
   * (num_rows x num_cols) and initialises each of the entries
   * to a constant value. This constant must be the same value type as
   * is specified in the matrix template.
   *
   * \param num_rows  the number of rows in the resulting matrix
   * \param num_cols  the number of columns in the resulting matrix
   * \param the_value a constant value to be assigned to all matrix entries
   */
  Matrix(std::size_t num_rows, std::size_t num_cols, double the_value);

  /**
   * Constructor: Creates a matrix with the specified dimensions
   * (number_of_rows x number_of_columns) from a std::vector. Note that this
   * does a copy construction such that the vector of values that was passed
   * in remains unchanged. This is not as fast as a move construction via
   * std::move, especially for large objects, but I suspect that this is the
   * behaviour most users would expect. The vector of values must have the
   * same value type as is specified in the matrix template.
   *
   * \param num_rows   the number of rows in the resulting matrix
   * \param num_cols the number of columns in the resulting matrix
   * \param the_values a constant value to be assigned to all matrix entries
   */
  Matrix(std::size_t num_rows, std::size_t num_cols,
      const std::vector<double> &the_values);

  /**
   * This method assigns and resizes the matrix to the corrected dimensions
   *
   * \param num_rows  the new number of rows in the resulting matrix
   * \param num_cols  the new number of columns in the resulting matrix
   * \param the_value  the value to be inserted into the matrix
   */
  void Assign(std::size_t num_rows, std::size_t num_cols, double the_value);

  /**
   * This method assigns and resizes the matrix to the corrected dimensions
   *
   * \param num_rows  the new number of rows in the resulting matrix
   * \param num_cols  the new number of columns in the resulting matrix
   * \param matrix_values the matrix vector to be inserted into the matrix
   */
  void Assign(std::size_t num_rows, std::size_t num_cols,
    const std::vector<double> &matrix_values);

  /**
   * Get the value stored in the matrix at the position (row,column).
   * Note: The matrix coordinates start from (0,0).
   *
   * \param row the row number of the requested value
   * \param col the column number of the requested value
   * \return the value of the matrix at position (row,column)
   */
  double GetValue(std::size_t row, std::size_t col) const noexcept;

  /**
   * Set the value stored in the matrix at the position (row,column).
   * Note: The matrix coordinates start from (0,0).
   *
   * \param row the row number of the requested value
   * \param col the column number of the requested value
   * \param the_value  the value to be placed at (row,column)
   */
  void SetValue(std::size_t row, std::size_t col, double the_value);

  /**
   * Get the number of columns in the stored matrix. The entries will be
   * stored from 0 to rows-1.
   *
   * \return the number columns in the matrix
   */
  std::size_t GetNumberOfColumns() const noexcept;

  /**
   * Get the number of rows in the stored matrix. The entries will be stored
   * from 0 to cols-1.
   *
   * \return the number rows in the matrix
   */
  std::size_t GetNumberOfRows() const noexcept;

  /**
   * Get the number of entries in the stored matrix.
   *
   * \return the number entries in the matrix
   */
  std::size_t Size() const noexcept;

  /**
   * Get the number of entries in the stored matrix.
   *
   * \return the number entries in the matrix
   */
  bool Empty() const noexcept;

  /**
   * Add a given value to a matrix diagonal. If offset is 0 (the default
   * value), this method will set the values along the main (leading) diagonal.
   * Otherwise, if an offset of "n" is specified, it will set the values along
   * the nth sub-diagonal. By default, the diagonal values will be set in the
   * upper triangular part of the matrix, meaning that the offset represents a
   * number of columns. If the diagonal of interest is in the lower triangular
   * part of the matrix then just set the upper_triangular boolean to false. In
   * this case the offset represents a number of rows. If an offset is
   * specified that is larger than the matrix dimensions no values in the
   * matrix are changed. If the offset is zero it does not matter if the
   * upper_triangular argument is true or not. The value type for the new
   * values must be the same as the value type in the matrix. Note that the
   * matrix coordinates start from (0,0).
   *
   * \param value  the value to be set on the diagonal
   * \param offset the distance from the main diagonal
   *     (optional, default = 0 = main diagonal)
   * \param upper_triangle elect the direction for the offset
   *     (optional, default = true = upper triangle)
   */
  void AddConstantToDiagonal(double value, std::size_t offset = 0,
      bool upper_triangle = true);

  /**
   * Matrix entries stored row-wise in a 1D vector. Note that these values
   * are public to allow direct (fast) access. This is fine for reading, but
   * if you are writing here directly it is up to you to check you are within
   * the bounds of the matrix.
   */
  std::vector<double> values_;

 private:
  /** Number of rows in the matrix */
  std::size_t number_of_rows_;
  /** Number of columns in the matrix */
  std::size_t number_of_columns_;
};

// Here are some functions that operate on matrices and vectors that we need.
// For simplicity these have all been collected here.

/**
 * Calculate the inner product of a matrix and a vector: r = M*v. This function
 * will check if the operation is possible by size(m)%size(v)=0.
 *
 * \param m input matrix
 * \param v input vector
 * \return a vector containing m*v
 */
std::vector<double> InnerProduct(const Matrix &m, const std::vector<double> &v);

/**
 * This method performs an inner product of two vectors and returns the (scalar)
 * result.
 *
 * \param v1 the first vector
 * \param v2 the first vector
 * \return the inner product v1 dot v2
 */
double InnerProduct(const std::vector<double> &v1,
    const std::vector<double> &v2);

/**
 * This function takes in a matrix and performs an inner product of the matrix
 * with itself. Given a matrix M, this means it will calculate (M^T)M, where M^T
 * is the transpose of M. The boolean variable sets whether or not you are
 * passing in M (transpose = false) or M^T (transpose = true). Please note that
 * because of how the matrix data is stored (row-wise) that this function will
 * be faster if you pass in M^T compared to passing in M.
 *
 * \param m The matrix to be multiplied by itself
 * \param transpose True if the transposed matrix is passed in, false otherwise
 * \return A matrix containing (M^T)M
 */
Matrix InnerProduct(const Matrix &m, bool transpose);

/**
 * This method simply multiplies every element in the vector by a constant and
 * returns the result in a new vector
 *
 * \param v the vector to be multiplied
 * \param scale_factor the value by which all vector elements are multiplied
 * \return v*scale_factor
 */
std::vector<double> Scale(const std::vector<double> &v,
    const double scale_factor);

/**
 * This method simply adds or subtracts one vector to or from another.
 *
 * \param v1 the vector by which all elements is added or substracted
 *        against another parameter addend_b (i.e vector_a -/+ vector_b)
 * \param v2 vector used as mentioned in addend_a
 * \param isadd If isadd is true, do addition else subtraction.
 *        Is set true as default.
 * \return returns a resulting vector
 */
std::vector<double> AddSubtractVectors(const std::vector<double> &v1,
    const std::vector<double> &v2, const bool isadd = true);

/**
 * This function performs the (vector) operation:
 *   result = v1 + v2*scale_factor
 *
 * If the vectors are of different sizes this function will return an empty
 * vector
 *
 * \param v1 The vector which is not scaled
 * \param v2 The vector that is scaled
 * \param scale_factor The scale factor for the second vector
 * \return A vector containing v1 + v2*scale_factor
 */
std::vector<double> AddScaleVectors(const std::vector<double> &v1,
    const std::vector<double> &v2, const double scale_factor);

/**
 * This method find the transpose of matrix m
 *
 * \param m the matrix that need to be transposed
 * \return the transposed matrix.
 */
Matrix Transpose(Matrix &m);

/**
 * This method simply finds the L2 norm of a vector, example :
 * vector_a {1, 2, 3} and the L2 norm is = sqrt(1^2 + 2^2 + 3^2)
 *
 * \param v the transposed vector J(x).
 * \return returns the value of the normal.
 */
double FindL2NormOfVector(const std::vector<double> &v);

/**
 * This method simply finds the sum of the squared elements of a vector. If the
 * vector is (v1,v2,v3), the result will be v1*v1 + v2*v2 +v3*v3.
 *
 * \param v The input vector
 * \return The sum of squares
 */
double SumOfSquares(const std::vector<double> &v);

/**
 * This method finds the maximum value along the diagonal. Note that this only
 * works for square matrices at present
 * example : vector_a {1, 4, 3, 2}
 *           maximum along diagonal = 2.
 *
 * \param m the matrix required to find the maximum.
 * \return returns the maximum value found.
 *
 */
double MaxDiagonalMat(const Matrix &m);

}  // namespace Unfit

#endif
