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
#include <limits>
#include "UnitTest++.h"
#include "Bounds.hpp"

static const double max_positive_num
    {static_cast<double>(std::numeric_limits<float>::max())};
static const double max_negative_num
    {-static_cast<double>(std::numeric_limits<float>::max())};
static const double bound_test_tol {1e-16};

namespace Unfit
{
namespace UnitTests
{
SUITE(UnitTestBounds)
{
TEST(Bounds_Constructors)
{
  // Start with no bounds
  Unfit::Bounds no_bounds;
  CHECK_EQUAL(0u, no_bounds.GetNumberOfBounds());

  // Now set the number of bounds
  no_bounds.SetNumberOfBounds(5);
  std::vector<double> lower;
  std::vector<double> upper;
  no_bounds.GetBounds(lower, upper);
  for (auto bnd : lower) CHECK_CLOSE(max_negative_num, bnd, bound_test_tol);
  for (auto bnd : upper) CHECK_CLOSE(max_positive_num, bnd, bound_test_tol);

  // Start with some bounds
  Unfit::Bounds some_bounds(5);
  CHECK_EQUAL(5u, some_bounds.GetNumberOfBounds());
  some_bounds.GetBounds(lower, upper);
  for (auto bnd : lower) CHECK_CLOSE(max_negative_num, bnd, bound_test_tol);
  for (auto bnd : upper) CHECK_CLOSE(max_positive_num, bnd, bound_test_tol);
}

TEST(Bounds_SetAndReset)
{
  Unfit::Bounds some_bounds(5);
  CHECK_EQUAL(5u, some_bounds.GetNumberOfBounds());

  // Set bounds and check they are set
  some_bounds.SetBounds(2, 0.0, 1.0);
  std::vector<double> lower;
  std::vector<double> upper;
  some_bounds.GetBounds(lower, upper);
  CHECK_CLOSE(0.0, lower[2], bound_test_tol);
  CHECK_CLOSE(1.0, upper[2], bound_test_tol);

  // Reset bounds and check they have been reset
  some_bounds.ResetBounds();
  some_bounds.GetBounds(lower, upper);
  for (auto bnd : lower) CHECK_CLOSE(max_negative_num, bnd, bound_test_tol);
  for (auto bnd : upper) CHECK_CLOSE(max_positive_num, bnd, bound_test_tol);
}

TEST(Bounds_SetUpperAndLower)
{
  Unfit::Bounds some_bounds(5);
  CHECK_EQUAL(5u, some_bounds.GetNumberOfBounds());

  // Set bounds and check they are set
  some_bounds.SetLowerBound(2, 0.0);
  some_bounds.SetUpperBound(3, 1.0);
  std::vector<double> lower;
  std::vector<double> upper;
  some_bounds.GetBounds(lower, upper);
  CHECK_CLOSE(0.0, lower[2], bound_test_tol);
  CHECK_CLOSE(max_positive_num, upper[2], bound_test_tol);
  CHECK_CLOSE(max_negative_num, lower[3], bound_test_tol);
  CHECK_CLOSE(1.0, upper[3], bound_test_tol);
}

TEST(Bounds_OnePointWithinBounds)
{
  Unfit::Bounds some_bounds;
  std::vector<double> lower {0.0, 0.0, 0.0};
  std::vector<double> upper {1.0, 1.0, 1.0};
  some_bounds.SetBounds(lower, upper);

  CHECK(some_bounds.IsWithinBounds(1, 0.5));
  CHECK(!some_bounds.IsWithinBounds(1, -1.0));
  CHECK(!some_bounds.IsWithinBounds(1, 2.0));
  CHECK(!some_bounds.IsWithinBounds(5, 0.5));
}

TEST(Bounds_OnePointWithinOneSidedBounds)
{
  Unfit::Bounds some_bounds;
  std::vector<double> lower {0.0, 0.0, 0.0};
  std::vector<double> upper {1.0, 1.0, 1.0};
  some_bounds.SetBounds(lower, upper);

  CHECK(some_bounds.IsBelowUpperBound(1, 0.5));
  CHECK(some_bounds.IsAboveLowerBound(1, 0.5));
  CHECK(!some_bounds.IsBelowUpperBound(3, 0.5));
  CHECK(!some_bounds.IsAboveLowerBound(3, 0.5));
  CHECK(!some_bounds.IsBelowUpperBound(1, 2.0));
  CHECK(!some_bounds.IsAboveLowerBound(1, -1.0));
}

TEST(Bounds_PointsWithinBounds)
{
  Unfit::Bounds some_bounds;
  std::vector<double> lower {0.0, 0.0, 0.0};
  std::vector<double> upper {1.0, 1.0, 1.0};
  some_bounds.SetBounds(lower, upper);

  std::vector<double> point {0.5, 0.5, 0.5};
  CHECK(some_bounds.IsWithinBounds(point));
  point[1] = -1.0;
  CHECK(!some_bounds.IsWithinBounds(point));
  point[1] = 2.0;
  CHECK(!some_bounds.IsWithinBounds(point));
}

TEST(Bounds_PointWithinBoundsEdgeCases)
{
  Unfit::Bounds some_bounds(3);
  std::vector<double> point {max_positive_num, max_negative_num,
      max_positive_num};
  CHECK(some_bounds.IsWithinBounds(point));
  point[0] *= 10.0;
  CHECK(!some_bounds.IsWithinBounds(point));
  point[0] = max_positive_num;
  point[1] *= 10.0;
  CHECK(!some_bounds.IsWithinBounds(point));
}

TEST(Bounds_ClampPointWithinBounds)
{
  Unfit::Bounds some_bounds;
  std::vector<double> lower {0.0, 0.0, 0.0};
  std::vector<double> upper {1.0, 1.0, 1.0};
  some_bounds.SetBounds(lower, upper);

  // Already within bounds, should do nothing
  std::vector<double> point {0.5, 0.5, 0.5};
  some_bounds.ClampWithinBounds(point);
  CHECK_CLOSE(0.5, point[0], bound_test_tol);
  CHECK_CLOSE(0.5, point[1], bound_test_tol);
  CHECK_CLOSE(0.5, point[2], bound_test_tol);

  // Clamp to lower bound
  point[1] = -1.0;
  some_bounds.ClampWithinBounds(point);
  CHECK_CLOSE(0.5, point[0], bound_test_tol);
  CHECK_CLOSE(0.0, point[1], bound_test_tol);
  CHECK_CLOSE(0.5, point[2], bound_test_tol);

  // Clamp to upper bound
  point[1] = 2.0;
  some_bounds.ClampWithinBounds(point);
  CHECK_CLOSE(0.5, point[0], bound_test_tol);
  CHECK_CLOSE(1.0, point[1], bound_test_tol);
  CHECK_CLOSE(0.5, point[2], bound_test_tol);
}

TEST(Bounds_SetBoundsErrorCases)
{
  // Array sizes don't match
  Unfit::Bounds some_bounds(3);
  std::vector<double> lower {0.0, 0.0};
  std::vector<double> upper {1.0, 1.0, 1.0};
  CHECK(!some_bounds.SetBounds(lower, upper));

  // Negative infinite bound
  lower.push_back(-std::numeric_limits<double>::infinity());
  CHECK(!some_bounds.SetBounds(lower, upper));

  // Positive infinite bound
  lower[2] = 0.0;
  upper[1] = std::numeric_limits<double>::infinity();
  CHECK(!some_bounds.SetBounds(lower, upper));

  // Lower bound > upper bound
  upper[1] = -1.0;
  CHECK(!some_bounds.SetBounds(lower, upper));

  // Now check that the original bounds are intact
  some_bounds.GetBounds(lower, upper);
  for (auto bnd : lower) CHECK_CLOSE(max_negative_num, bnd, bound_test_tol);
  for (auto bnd : upper) CHECK_CLOSE(max_positive_num, bnd, bound_test_tol);
}

TEST(Bounds_SetOneBoundErrorCases)
{
  // Array sizes don't match
  Unfit::Bounds some_bounds(3);

  // Negative infinite bound
  CHECK(!some_bounds.SetBounds(2, -std::numeric_limits<double>::infinity(), 1.0));

  // Positive infinite bound
  CHECK(!some_bounds.SetBounds(2, 0.0, std::numeric_limits<double>::infinity()));

  // Lower bound > upper bound
  CHECK(!some_bounds.SetBounds(2, 1.0, 0.0));

  // Check adding an index greater than the vector size (should scale)
  CHECK(some_bounds.SetBounds(5, max_negative_num, max_positive_num));
  CHECK_EQUAL(6u, some_bounds.GetNumberOfBounds());

  // Now check that the original and extended bounds are intact
  std::vector<double> lower;
  std::vector<double> upper;
  some_bounds.GetBounds(lower, upper);
  for (auto bnd : lower) CHECK_CLOSE(max_negative_num, bnd, bound_test_tol);
  for (auto bnd : upper) CHECK_CLOSE(max_positive_num, bnd, bound_test_tol);
}

TEST(Bounds_GetSetLowerBound)
{
  Unfit::Bounds some_bounds(2);
  // Default bounds
  CHECK_CLOSE(max_negative_num, some_bounds.GetLowerBound(0), bound_test_tol);
  CHECK_CLOSE(max_negative_num, some_bounds.GetLowerBound(1), bound_test_tol);
  // Invalid indices
  CHECK_CLOSE(max_negative_num, some_bounds.GetLowerBound(2), bound_test_tol);
  CHECK_CLOSE(max_negative_num, some_bounds.GetLowerBound(3), bound_test_tol);

  some_bounds.SetLowerBound(0, 1.0);
  some_bounds.SetLowerBound(1, 0.0);
  // Default bounds
  CHECK_CLOSE(1.0, some_bounds.GetLowerBound(0), bound_test_tol);
  CHECK_CLOSE(0.0, some_bounds.GetLowerBound(1), bound_test_tol);
  CHECK_CLOSE(max_positive_num, some_bounds.GetUpperBound(0), bound_test_tol);
  CHECK_CLOSE(max_positive_num, some_bounds.GetUpperBound(1), bound_test_tol);
  // Invalid indices
  CHECK_CLOSE(max_negative_num, some_bounds.GetLowerBound(2), bound_test_tol);
  CHECK_CLOSE(max_negative_num, some_bounds.GetLowerBound(3), bound_test_tol);
}

TEST(Bounds_GetSetUpperBound)
{
  Unfit::Bounds some_bounds(2);
  // Default bounds
  CHECK_CLOSE(max_positive_num, some_bounds.GetUpperBound(0), bound_test_tol);
  CHECK_CLOSE(max_positive_num, some_bounds.GetUpperBound(1), bound_test_tol);
  // Invalid indices
  CHECK_CLOSE(max_positive_num, some_bounds.GetUpperBound(2), bound_test_tol);
  CHECK_CLOSE(max_positive_num, some_bounds.GetUpperBound(3), bound_test_tol);

  some_bounds.SetUpperBound(0, 1.0);
  some_bounds.SetUpperBound(1, 0.0);
  // Default bounds
  CHECK_CLOSE(1.0, some_bounds.GetUpperBound(0), bound_test_tol);
  CHECK_CLOSE(0.0, some_bounds.GetUpperBound(1), bound_test_tol);
  CHECK_CLOSE(max_negative_num, some_bounds.GetLowerBound(0), bound_test_tol);
  CHECK_CLOSE(max_negative_num, some_bounds.GetLowerBound(1), bound_test_tol);
  // Invalid indices
  CHECK_CLOSE(max_positive_num, some_bounds.GetUpperBound(2), bound_test_tol);
  CHECK_CLOSE(max_positive_num, some_bounds.GetUpperBound(3), bound_test_tol);
}
}  // suite UnitTestBounds
}  // namespace UnitTests
}  // namespace Unfit
