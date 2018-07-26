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
#include <iostream>
#include <vector>
#include "GenericNDModel.hpp"
#include "GenericNDCostFunction.hpp"

namespace Unfit
{
GenericNDCostFunction::GenericNDCostFunction(GenericNDModel &model)
  : model_ {model},
    x_ {},
    y_ {},
    has_data_ {false}
{}

GenericNDCostFunction::GenericNDCostFunction(GenericNDModel &model,
    const std::vector<std::vector<double>> &x, const std::vector<double> &y)
  : model_ {model},
    x_ {x},
    y_ {y},
    has_data_ {false}
{}

std::vector<double> GenericNDCostFunction::operator()(
    const std::vector<double> &c)
{
  // Need to check the data as we don't check during the constructor
  if (!has_data_) CheckData(x_, y_);
  // CheckData can change has_data_, so we must check again in case it did
  if (has_data_) {
    auto model_result = model_(c, x_);
    // Want r = y - f(x), and get this via r = y, then r -= f(x)
    auto residuals = y_;
    for (auto i = 0u; i < residuals.size(); ++i) {
      residuals[i] -= model_result[i];
    }
    return residuals;
  } else {
    std::cout << "CheckData failed. Model has no valid data set" << std::endl;
    return std::vector<double>(); // empty vector
  }
}

bool GenericNDCostFunction::CheckData(const std::vector<std::vector<double>> &x,
    const std::vector<double> &y)
{
  has_data_ = false;
  if (x.empty()) return false;
  if (y.empty()) return false;
  if (x[0].empty()) return false;
  auto data_size = y.size();
  for (auto x_i : x) {
    if (x_i.size() != data_size) return false;
  }
  has_data_ = true;
  return true;
}

void GenericNDCostFunction::GetData(std::vector<std::vector<double>> &x,
    std::vector<double> &y)
{
  x = x_;
  y = y_;
}

bool GenericNDCostFunction::SetData(const std::vector<std::vector<double>> &x,
    const std::vector<double> &y)
{
  auto data_is_okay = CheckData(x, y);
  if (data_is_okay) {
    x_ = x;
    y_ = y;
    has_data_ = true;
  }
  return data_is_okay;
}

}  // namespace Unfit
