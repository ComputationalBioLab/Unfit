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
#include "CardiacAlphaNModel.hpp"
#include "GenericNDCostFunction.hpp"
#include "UnitTest++.h"

static const double tolerance {1.0e-6};

namespace Unfit
{
namespace UnitTests
{
SUITE(UnitTestGenericNDCostFunction)
{
TEST(GenericNDCostFunction_2DResidualCalculation)
{
  std::vector<double> alpha_n {0.0, 0.05, 0.1, 0.15, 0.2};
  std::vector<std::vector<double>> vm {{-80.0, -40.0, 0.0, 40.0, 80.0}};
  Examples::CardiacAlphaNModel c_alpha_n;
  GenericNDCostFunction cost_func(c_alpha_n, vm, alpha_n);
  std::vector<double> c {1.0, -0.01, 2.0};
  auto residuals = cost_func(c);
  CHECK(residuals.size() == vm[0].size());
  CHECK_CLOSE(-0.057324, residuals[0], tolerance);
  CHECK_CLOSE(-0.033173, residuals[1], tolerance);
  CHECK_CLOSE(-0.019203, residuals[2], tolerance);
  CHECK_CLOSE(-0.017982, residuals[3], tolerance);
  CHECK_CLOSE(-0.031475, residuals[4], tolerance);
}

TEST(GenericNDCostFunction_SetData2D)
{
  std::vector<double> alpha_n {0.0, 0.05, 0.1, 0.15, 0.2};
  std::vector<std::vector<double>> vm {{-80.0, -40.0, 0.0, 40.0, 80.0}};
  Examples::CardiacAlphaNModel c_alpha_n;
  GenericNDCostFunction cost_func(c_alpha_n, vm, alpha_n);
  vm.clear();
  vm = {{-60.0, -40.0, -20.0, 0.0, 20.0}};
  auto rc = cost_func.SetData(vm, alpha_n);  // change vm
  CHECK(rc);
  std::vector<double> c {1.0, -0.01, 2.0};
  auto residuals = cost_func(c);
  CHECK(residuals.size() == vm[0].size());
  CHECK_CLOSE(-0.069138, residuals[0], tolerance);
  CHECK_CLOSE(-0.033173, residuals[1], tolerance);
  CHECK_CLOSE( 0.000250, residuals[2], tolerance);
  CHECK_CLOSE( 0.030797, residuals[3], tolerance);
  CHECK_CLOSE( 0.058149, residuals[4], tolerance);
  alpha_n.clear();
  alpha_n = {0.1, 0.12, 0.14, 0.16, 0.18};
  rc = cost_func.SetData(vm, alpha_n);  // then change alpha_n
  CHECK(rc);
  residuals = cost_func(c);
  CHECK(residuals.size() == vm[0].size());
  CHECK_CLOSE(0.030862, residuals[0], tolerance);
  CHECK_CLOSE(0.036827, residuals[1], tolerance);
  CHECK_CLOSE(0.040250, residuals[2], tolerance);
  CHECK_CLOSE(0.040797, residuals[3], tolerance);
  CHECK_CLOSE(0.038149, residuals[4], tolerance);
  vm[0].push_back(40.0);
  alpha_n.push_back(0.2);
  rc = cost_func.SetData(vm, alpha_n); // then change both
  CHECK(rc);
  residuals = cost_func(c);
  CHECK(residuals.size() == vm[0].size());
  CHECK_CLOSE(0.030862, residuals[0], tolerance);
  CHECK_CLOSE(0.036827, residuals[1], tolerance);
  CHECK_CLOSE(0.040250, residuals[2], tolerance);
  CHECK_CLOSE(0.040797, residuals[3], tolerance);
  CHECK_CLOSE(0.038149, residuals[4], tolerance);
  CHECK_CLOSE(0.032018, residuals[5], tolerance);
}

TEST(Generic2DCostFunction_SetDataEmptyData)
{
  std::vector<double> alpha_n {0.0, 0.05, 0.1, 0.15, 0.2};
  std::vector<std::vector<double>> vm {{-80.0, -40.0, 0.0, 40.0, 80.0}};
  Examples::CardiacAlphaNModel c_alpha_n;
  GenericNDCostFunction cost_func(c_alpha_n, vm, alpha_n);
  vm.clear();
  auto rc = cost_func.SetData(vm, alpha_n);  // vm is empty
  CHECK(!rc);
  vm = {{-80.0, -40.0, 0.0, 40.0, 80.0}};
  alpha_n.clear();
  rc = cost_func.SetData(vm, alpha_n);  // alpha is empty
  CHECK(!rc);
  vm[0].clear();
  rc = cost_func.SetData(vm, alpha_n);  // both are empty
  CHECK(!rc);
  alpha_n = {0.0, 0.05, 0.1, 0.15, 0.2};
  vm = {{-80.0, -40.0, 0.0, 40.0, 80.0}};
  rc = cost_func.SetData(vm, alpha_n);  // both are okay
  CHECK(rc);
}

TEST(Generic2DCostFunction_SetDataWrongSize)
{
  std::vector<double> alpha_n {0.0, 0.05, 0.1, 0.15, 0.2};
  std::vector<std::vector<double>> vm {{-80.0, -40.0, 0.0, 40.0, 80.0}};
  Examples::CardiacAlphaNModel c_alpha_n;
  GenericNDCostFunction cost_func(c_alpha_n, vm, alpha_n);
  vm[0].pop_back();
  auto rc = cost_func.SetData(vm, alpha_n);  // vm shorter
  CHECK(!rc);
  vm[0].push_back(0.2);
  vm[0].push_back(0.25);
  rc = cost_func.SetData(vm, alpha_n);  // vm is longer
  CHECK(!rc);
  vm[0].pop_back();
  rc = cost_func.SetData(vm, alpha_n);  // both are okay
  CHECK(rc);
}
}  // suite UnitTestGeneric2DCostFunction
}  // namespace UnitTests
}  // namespace Unfit
