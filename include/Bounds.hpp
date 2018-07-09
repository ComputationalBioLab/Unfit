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
#ifndef UNFIT_INCLUDE_BOUNDS_HPP_
#define UNFIT_INCLUDE_BOUNDS_HPP_

#include <vector>

namespace Unfit
{
/**
 * This class is designed to handle all of the bounds work for Unfit. To ensure
 * safety, internally a pair of bounds should be created for all variables, even
 * if there is no desire to set them to specific values. This is as easy as
 * passing the number of variables to the constructor (if known), or later
 * calling SetNumberOfBounds which will also create the default bounds for you.
 * Remember that all indices start from zero in C++, so that is the location of
 * the first bound. Note that bounds are always set (or removed) in pairs, so
 * if there is an upper bound there will always be a lower bound.
 */
class Bounds
{
 public:
  /**
   * Create an empty Bounds object where no bounds are set.
   */
  Bounds();

  /**
   * Create a Bounds object and sets number_of_bounds variables to have a lower
   * bound at equal to the maximum negative float and an upper bound equal to
   * the maximum positive float (stored as doubles). The reason for this is that
   * on some platforms using the minimum and maximum doubles has caused issues.
   *
   * \param number_of_bounds The number of variables in the problem
   */
  Bounds(std::size_t number_of_bounds);

  /**
   * Checks if the coordinates of the point provided to the method are within
   * bounds of the domain. The check is INCLUSIVE. This means that if the
   * value of a variable is on the domain boundary it is considered to be within
   * the domain. If any of the variables are outside the bounds it will change
   * that coordinate such that it lies on the bound.
   *
   * \param point the coordinates of the point to check
   */
  void ClampWithinBounds(std::vector<double> &point);

  /**
   * Here, the current values for all of the upper and lower bounds are
   * written into the upper_bound and lower_bound arguments, respectively.
   *
   * \param lower_bound (input) vector containing the lower bounds;
   * \param upper_bound (input) vector containing the upper bounds;
   */
  void GetBounds(std::vector<double> &lower_bound,
      std::vector<double> &upper_bound) const;

  /**
   * Returns the number of variables that are currently bounded. This should
   * be the same as the number of variables in the problem. If not, use
   * SetNumberOfBounds to fix this.
   *
   * \return The number of bounded variables
   */
  std::size_t GetNumberOfBounds() const noexcept;

  /**
   * Checks if the point provided to the method is above the lower bound on that
   * variable. The check is INCLUSIVE. This means that if the value of a
   * variable is on the domain boundary it is considered to be within
   * the domain. An invalid index will return false.
   *
   * \param index the location of the point to check
   * \param point the value of the point to check
   *
   * \return true if the point is within the bounds \n
   *         false if the point is not within the bounds
   */
  bool IsAboveLowerBound(std::size_t index, double point) const noexcept;

  /**
   * Checks if the point provided to the method is below the upper bound on that
   * variable. The check is INCLUSIVE. This means that if the value of a
   * variable is on the domain boundary it is considered to be within
   * the domain. An invalid index will return false.
   *
   * \param index the location of the point to check
   * \param point the value of the point to check
   *
   * \return true if the point is within the bounds \n
   *         false if the point is not within the bounds
   */
  bool IsBelowUpperBound(std::size_t index, double point) const noexcept;

  /**
   * Checks if the point provided to the method is within the bounds of
   * the variable. The check is INCLUSIVE. This means that if the value of a
   * variable is on the domain boundary it is considered to be within
   * the domain. An invalid index will return false.
   *
   * \param index the location of the point to check
   * \param point the value of the point to check
   *
   * \return true if the point is within the bounds \n
   *         false if the point is not within the bounds
   */
  bool IsWithinBounds(std::size_t index, double point) const noexcept;

  /**
   * Checks if the coordinates of the point provided to the method are within
   * bounds of the domain. The check is INCLUSIVE. This means that if the
   * value of a variable is on the domain boundary it is considered to be within
   * the domain.
   *
   * \param point the coordinates of the point to check
   *
   * \return true if the coordinates are within the bounds \n
   *         false if the coordinates are not within the bounds
   */
  bool IsWithinBounds(const std::vector<double> &point) const noexcept;

  /**
   * This method keeps the number of bounds unchanged, but resets them such
   * that all of the lower bounds are the maximum negative float and all of the
   * upper bounds are the maximum positive float (stored as doubles).
   */
  void ResetBounds();

  /**
   * Sets the upper and lower bound at the position "index". If setting the
   * bound fails, the original state of the bounds will be preserved (including
   * remaining undefined if they have not been set previously).
   *
   * \param index The position of the variable on which the bound will be set
   * \param lower_bound The value of the lower bound
   * \param upper_bound The value of the upper bound
   *
   * \return true = Sucess \n
   *         false = Failure; invalid bound or bounds
   */
  bool SetBounds(std::size_t index, double lower_bound, double upper_bound);

  /**
   * Sets the upper and lower bounds at all positions. All existing bounds will
   * be overwritten if this method succeeds. This includes changing the number
   * of bounds, if necessary. The number of upper bounds must be equal to the
   * number of lower bounds, which should be equal to the number of variables
   * you have. If setting the bounds fails, the original state of the bounds
   * will be preserved (including remaining undefined if they have not been set
   * previously).
   *
   * \param lower_bound A vector containing the lower bounds
   * \param upper_bound A vector containing the upper bounds
   *
   * \return true = Sucess \n
   *         false = Failure; invalid bound or bounds
   */
  bool SetBounds(const std::vector<double> &lower_bound,
      const std::vector<double> &upper_bound);

  /**
   * Set the number of bounds in the problem. This should be the same as the
   * number of variables. If number_of_bounds is equal to the current number of
   * bounds no updates will occur (existing bounds will be preserved). If
   * number_of_bounds is less that the current number of bounds, the bounds
   * above number_of_bounds will be deleted (other bounds will be preserved). If
   * number_of_bounds is greater than the current number of bounds, new bounds
   * will be appended (lower = -max; upper = +max) to the existing set of bounds
   * (if any). NOTE: An upper and a lower bound on the same variable only count
   * as one bound.
   *
   * \param number_of_bounds The number of variables in the problem
   */
  void SetNumberOfBounds(std::size_t number_of_bounds);

  /**
   * Returns the lower bound value of given index. The index should
   * be the smaller than the size of the lower bound vector. If not,
   * method assumes that there is no contraints and returns the
   * largest negative double precision number.
   *
   * \param index The position of the variable on which the bound will be set
   *
   * \return The value of lower bound of given index
   */
  double GetLowerBound(std::size_t index) const noexcept;

  /**
   * Set the lower bound at the position "index". The corresponding upper bound
   * will be set to the largest positive double precision number. If setting the
   * bound fails, the original state of the bounds will be preserved (including
   * remaining undefined if they have not been set previously).
   *
   * \param index The position of the variable on which the bound will be set
   * \param lower_bound The value of the lower bound
   *
   * \return true = Sucess \n
   *         false = Failure; invalid upper bound
   */
  bool SetLowerBound(std::size_t index, double lower_bound);

  /**
   * Returns the upper bound value of given index. The index should
   * be the smaller than the size of the upper bound vector. If not,
   * method assumes that there is no contraints and returns the
   * largest positive double precision number.
   *
   * \param index The position of the variable on which the bound will be set
   *
   * \return The value of upper bound of given index
   */
  double GetUpperBound(std::size_t index) const noexcept;

  /**
   * Set the upper bound at the position "index". The corresponding lower bound
   * will be set to the largest negative double precision number. If setting the
   * bound fails, the original state of the bounds will be preserved (including
   * remaining undefined if they have not been set previously).
   *
   * \param index The position of the variable on which the bound will be set
   * \param upper_bound The value of the upper bound
   *
   * \return true = Sucess \n
   *         false = Failure; invalid upper bound
   */
  bool SetUpperBound(std::size_t index, double upper_bound);

 private:
  /**Stores the upper bounds of the variables*/
  std::vector<double> upper_bound_;
  /**Stores the lower bounds of the variables*/
  std::vector<double> lower_bound_;
};

}  // namespace Unfit

#endif
