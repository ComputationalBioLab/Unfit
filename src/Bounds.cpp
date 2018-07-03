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
#include <limits>
#include <vector>
#include "Bounds.hpp"

namespace Unfit
{
static const double max_positive_num {std::numeric_limits<double>::max()};
static const double max_negative_num {-std::numeric_limits<double>::max()};

Bounds::Bounds()
  : upper_bound_ {},
    lower_bound_ {}
{}

Bounds::Bounds(std::size_t number_of_bounds)
  : upper_bound_(number_of_bounds, max_positive_num),
    lower_bound_(number_of_bounds, max_negative_num)
{}

void Bounds::GetBounds(std::vector<double> &lower_bound,
    std::vector<double> &upper_bound) const
{
  lower_bound = lower_bound_;
  upper_bound = upper_bound_;
}

void Bounds::ClampWithinBounds(std::vector<double> &point)
{
  for (std::size_t i = 0; i < upper_bound_.size(); ++i) {
    if (point[i] < lower_bound_[i]) point[i] = lower_bound_[i];
    if (point[i] > upper_bound_[i]) point[i] = upper_bound_[i];
  }
}

std::size_t Bounds::GetNumberOfBounds() const noexcept
{
  return upper_bound_.size();
}

bool Bounds::IsAboveLowerBound(std::size_t index, double point) const noexcept
{
  if (index >= lower_bound_.size()) return false;
  if (point < lower_bound_[index]) return false;
  return true;
}

bool Bounds::IsBelowUpperBound(std::size_t index, double point) const noexcept
{
  if (index >= upper_bound_.size()) return false;
  if (point > upper_bound_[index]) return false;
  return true;
}

bool Bounds::IsWithinBounds(std::size_t index, double point) const noexcept
{
  if (index >= upper_bound_.size()) return false;
  if (point < lower_bound_[index] || point > upper_bound_[index]) return false;
  return true;
}

bool Bounds::IsWithinBounds(const std::vector<double> &point) const noexcept
{
  for (std::size_t i = 0; i < upper_bound_.size(); ++i) {
    if (point[i] < lower_bound_[i] || point[i] > upper_bound_[i]) return false;
  }
  return true;
}

void Bounds::ResetBounds()
{
  for (auto &bnd : upper_bound_) bnd = max_positive_num;
  for (auto &bnd : lower_bound_) bnd = max_negative_num;
}

bool Bounds::SetBounds(std::size_t index, double lower_bound,
    double upper_bound)
{
  // Check to see the bounds are reasonable, else do nothing
  if (!std::isfinite(lower_bound)) return false;
  if (!std::isfinite(upper_bound)) return false;
  if (lower_bound > upper_bound) return false;
  // Expand the bounds arrays, if needed
  if (index >= upper_bound_.size()) SetNumberOfBounds(index+1);
  // Set the requested bounds
  lower_bound_[index] = lower_bound;
  upper_bound_[index] = upper_bound;
  return true;
}

bool Bounds::SetBounds(const std::vector<double> &lower_bound,
  const std::vector<double> &upper_bound)
{
  // Check the bounds are all reasonable, else do nothing
  if (lower_bound.size() != upper_bound.size()) return false;
  for (std::size_t i = 0; i < lower_bound.size(); ++i) {
    if (!std::isfinite(lower_bound[i])) return false;
    if (!std::isfinite(upper_bound[i])) return false;
    if (lower_bound[i] >= upper_bound[i]) return false;
  }
  // Set the values of the bounds as given by the arguments
  lower_bound_ = lower_bound;
  upper_bound_ = upper_bound;
  return true;
}

void Bounds::SetNumberOfBounds(std::size_t number_of_bounds)
{
  upper_bound_.resize(number_of_bounds, max_positive_num);
  lower_bound_.resize(number_of_bounds, max_negative_num);
}

double Bounds::GetLowerBound(std::size_t index) const noexcept
{
  if (index < lower_bound_.size()) return lower_bound_[index];
  return max_negative_num;
}

bool Bounds::SetLowerBound(std::size_t index, double lower_bound)
{
  return SetBounds(index, lower_bound, max_positive_num);
}

double Bounds::GetUpperBound(std::size_t index) const noexcept
{
  if (index < upper_bound_.size()) return upper_bound_[index];
  return max_positive_num;
}

bool Bounds::SetUpperBound(std::size_t index, double upper_bound)
{
  return SetBounds(index, max_negative_num, upper_bound);
}

}  // namespace Unfit

