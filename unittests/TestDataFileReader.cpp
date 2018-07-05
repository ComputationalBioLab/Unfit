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
#include <string>
#include <vector>
#include "DataFileReader.hpp"
#include "UnitTest++.h"

static const double zero_tol(1.0e-16);
static const std::string path("./unittests/data/");

namespace Unfit
{
namespace UnitTests
{
SUITE(UnitTestDataFileReader)
{
// *** Tests for the SplitLine method ***
TEST(DataFileReader_SplitStringOnSpaceDelimiter)
{
  std::string test_string = "this is a test 1 2 3";
  std::vector<std::string> words;
  Unfit::DataFileReader<int> dfr;
  dfr.SplitLine(test_string, " ", words);
  CHECK_EQUAL(7u, words.size());
  CHECK_EQUAL("this", words[0]);
  CHECK_EQUAL("is", words[1]);
  CHECK_EQUAL("a", words[2]);
  CHECK_EQUAL("test", words[3]);
  CHECK_EQUAL("1", words[4]);
  CHECK_EQUAL("2", words[5]);
  CHECK_EQUAL("3", words[6]);
}

TEST(DataFileReader_SplitStringOnSpaceDelimiterWithLeadingSpaces)
{
  std::string test_string = "   this is a test 1 2 3";
  std::vector<std::string> words;
  Unfit::DataFileReader<int> dfr;
  dfr.SplitLine(test_string, " ", words);
  CHECK_EQUAL(7u, words.size());
  CHECK_EQUAL("this", words[0]);
  CHECK_EQUAL("is", words[1]);
  CHECK_EQUAL("a", words[2]);
  CHECK_EQUAL("test", words[3]);
  CHECK_EQUAL("1", words[4]);
  CHECK_EQUAL("2", words[5]);
  CHECK_EQUAL("3", words[6]);
}

TEST(DataFileReader_SplitStringOnSpaceDelimiterWithTrailingSpaces)
{
  std::string test_string = "this is a test 1 2 3   ";
  std::vector<std::string> words;
  Unfit::DataFileReader<int> dfr;
  dfr.SplitLine(test_string, " ", words);
  CHECK_EQUAL(7u, words.size());
  CHECK_EQUAL("this", words[0]);
  CHECK_EQUAL("is", words[1]);
  CHECK_EQUAL("a", words[2]);
  CHECK_EQUAL("test", words[3]);
  CHECK_EQUAL("1", words[4]);
  CHECK_EQUAL("2", words[5]);
  CHECK_EQUAL("3", words[6]);
}

TEST(DataFileReader_SplitStringOnCommaDelimiter)
{
  std::string test_string = "this,is,a,test,1,2,3";
  std::vector<std::string> words;
  Unfit::DataFileReader<int> dfr;
  dfr.SplitLine(test_string, ",", words);
  CHECK_EQUAL(7u, words.size());
  CHECK_EQUAL("this", words[0]);
  CHECK_EQUAL("is", words[1]);
  CHECK_EQUAL("a", words[2]);
  CHECK_EQUAL("test", words[3]);
  CHECK_EQUAL("1", words[4]);
  CHECK_EQUAL("2", words[5]);
  CHECK_EQUAL("3", words[6]);
}

TEST(DataFileReader_SplitStringOnCommaDelimiterWithLeadingCommas)
{
  std::string test_string = ",,,this,is,a,test,1,2,3";
  std::vector<std::string> words;
  Unfit::DataFileReader<int> dfr;
  dfr.SplitLine(test_string, ",", words);
  CHECK_EQUAL(7u, words.size());
  CHECK_EQUAL("this", words[0]);
  CHECK_EQUAL("is", words[1]);
  CHECK_EQUAL("a", words[2]);
  CHECK_EQUAL("test", words[3]);
  CHECK_EQUAL("1", words[4]);
  CHECK_EQUAL("2", words[5]);
  CHECK_EQUAL("3", words[6]);
}

TEST(DataFileReader_SplitStringOnCommaDelimiterWithTrailingCommas)
{
  std::string test_string = "this,is,a,test,1,2,3,,,";
  std::vector<std::string> words;
  Unfit::DataFileReader<int> dfr;
  dfr.SplitLine(test_string, ",", words);
  CHECK_EQUAL(7u, words.size());
  CHECK_EQUAL("this", words[0]);
  CHECK_EQUAL("is", words[1]);
  CHECK_EQUAL("a", words[2]);
  CHECK_EQUAL("test", words[3]);
  CHECK_EQUAL("1", words[4]);
  CHECK_EQUAL("2", words[5]);
  CHECK_EQUAL("3", words[6]);
}

TEST(DataFileReader_SplitStringOnTabDelimiter)
{
  std::string test_string = "this\tis\ta\ttest\t1\t2\t3";
  std::vector<std::string> words;
  Unfit::DataFileReader<int> dfr;
  dfr.SplitLine(test_string, "\t", words);
  CHECK_EQUAL(7u, words.size());
  CHECK_EQUAL("this", words[0]);
  CHECK_EQUAL("is", words[1]);
  CHECK_EQUAL("a", words[2]);
  CHECK_EQUAL("test", words[3]);
  CHECK_EQUAL("1", words[4]);
  CHECK_EQUAL("2", words[5]);
  CHECK_EQUAL("3", words[6]);
}

TEST(DataFileReader_SplitStringOnTabDelimiterWithLeadingTabs)
{
  std::string test_string = "\t\t\tthis\tis\ta\ttest\t1\t2\t3";
  std::vector<std::string> words;
  Unfit::DataFileReader<int> dfr;
  dfr.SplitLine(test_string, "\t", words);
  CHECK_EQUAL(7u, words.size());
  CHECK_EQUAL("this", words[0]);
  CHECK_EQUAL("is", words[1]);
  CHECK_EQUAL("a", words[2]);
  CHECK_EQUAL("test", words[3]);
  CHECK_EQUAL("1", words[4]);
  CHECK_EQUAL("2", words[5]);
  CHECK_EQUAL("3", words[6]);
}

TEST(DataFileReader_SplitStringOnTabDelimiterWithTrailingTabs)
{
  std::string test_string = "this\tis\ta\ttest\t1\t2\t3\t\t\t";
  std::vector<std::string> words;
  Unfit::DataFileReader<int> dfr;
  dfr.SplitLine(test_string, "\t", words);
  CHECK_EQUAL(7u, words.size());
  CHECK_EQUAL("this", words[0]);
  CHECK_EQUAL("is", words[1]);
  CHECK_EQUAL("a", words[2]);
  CHECK_EQUAL("test", words[3]);
  CHECK_EQUAL("1", words[4]);
  CHECK_EQUAL("2", words[5]);
  CHECK_EQUAL("3", words[6]);
}

TEST(DataFileReader_SplitStringOnNewlineDelimiter)
{
  std::string test_string = "this\nis\na\ntest\n1\n2\n3";
  std::vector<std::string> words;
  Unfit::DataFileReader<int> dfr;
  dfr.SplitLine(test_string, "\n ,", words);
  CHECK_EQUAL(7u, words.size());
  CHECK_EQUAL("this", words[0]);
  CHECK_EQUAL("is", words[1]);
  CHECK_EQUAL("a", words[2]);
  CHECK_EQUAL("test", words[3]);
  CHECK_EQUAL("1", words[4]);
  CHECK_EQUAL("2", words[5]);
  CHECK_EQUAL("3", words[6]);
}

TEST(DataFileReader_SplitStringOnNewlineDelimiterWithLeadingNewlines)
{
  std::string test_string = "\n\n\nthis\nis\na\ntest\n1\n2\n3";
  std::vector<std::string> words;
  Unfit::DataFileReader<int> dfr;
  dfr.SplitLine(test_string, "\n", words);
  CHECK_EQUAL(7u, words.size());
  CHECK_EQUAL("this", words[0]);
  CHECK_EQUAL("is", words[1]);
  CHECK_EQUAL("a", words[2]);
  CHECK_EQUAL("test", words[3]);
  CHECK_EQUAL("1", words[4]);
  CHECK_EQUAL("2", words[5]);
  CHECK_EQUAL("3", words[6]);
}

TEST(DataFileReader_SplitStringOnNewlineDelimiterWithTrailingNewlines)
{
  std::string test_string = "this\nis\na\ntest\n1\n2\n3\n\n\n";
  std::vector<std::string> words;
  Unfit::DataFileReader<int> dfr;
  dfr.SplitLine(test_string, "\n", words);
  CHECK_EQUAL(7u, words.size());
  CHECK_EQUAL("this", words[0]);
  CHECK_EQUAL("is", words[1]);
  CHECK_EQUAL("a", words[2]);
  CHECK_EQUAL("test", words[3]);
  CHECK_EQUAL("1", words[4]);
  CHECK_EQUAL("2", words[5]);
  CHECK_EQUAL("3", words[6]);
}

TEST(DataFileReader_SplitStringOnMultipleDelimiters)
{
  std::string test_string = ",,,this\t is\t,a test,1\t2 \t3,\t\n\t";
  std::vector<std::string> words;
  Unfit::DataFileReader<int> dfr;
  dfr.SplitLine(test_string, "\t ,\n", words);
  CHECK_EQUAL(7u, words.size());
  CHECK_EQUAL("this", words[0]);
  CHECK_EQUAL("is", words[1]);
  CHECK_EQUAL("a", words[2]);
  CHECK_EQUAL("test", words[3]);
  CHECK_EQUAL("1", words[4]);
  CHECK_EQUAL("2", words[5]);
  CHECK_EQUAL("3", words[6]);
}

TEST(DataFileReader_SplitZeroLengthString)
{
  std::string test_string = "";
  std::vector<std::string> words;
  Unfit::DataFileReader<int> dfr;
  dfr.SplitLine(test_string, "\t ,", words);
  CHECK_EQUAL(0u, words.size());
}

// *** Tests for the ReadFile method ***

TEST(DataFileReader_ReadSingleColumnOfIntegers)
{
  std::string file_name = path + "SingleColumnDataFile.txt";
  Unfit::DataFileReader<int> dfr;
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(9u, dfr.data.size());
  CHECK_EQUAL(1u, dfr.data[0].size());
  CHECK_EQUAL(1u, dfr.data[8].size());
  CHECK_EQUAL(1, dfr.data[0][0]);
  CHECK_EQUAL(2, dfr.data[1][0]);
  CHECK_EQUAL(3, dfr.data[2][0]);
  CHECK_EQUAL(4, dfr.data[3][0]);
  CHECK_EQUAL(5, dfr.data[4][0]);
  CHECK_EQUAL(6, dfr.data[5][0]);
  CHECK_EQUAL(7, dfr.data[6][0]);
  CHECK_EQUAL(8, dfr.data[7][0]);
  CHECK_EQUAL(9, dfr.data[8][0]);
}

TEST(DataFileReader_ReadSingleColumnOfDoubles)
{
  std::string file_name = path + "SingleColumnDataFile.txt";
  Unfit::DataFileReader<double> dfr;
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(9u, dfr.data.size());
  CHECK_EQUAL(1u, dfr.data[0].size());
  CHECK_EQUAL(1u, dfr.data[8].size());
  CHECK_CLOSE(1.0, dfr.data[0][0], zero_tol);
  CHECK_CLOSE(2.0, dfr.data[1][0], zero_tol);
  CHECK_CLOSE(3.0, dfr.data[2][0], zero_tol);
  CHECK_CLOSE(4.0, dfr.data[3][0], zero_tol);
  CHECK_CLOSE(5.0, dfr.data[4][0], zero_tol);
  CHECK_CLOSE(6.0, dfr.data[5][0], zero_tol);
  CHECK_CLOSE(7.0, dfr.data[6][0], zero_tol);
  CHECK_CLOSE(8.0, dfr.data[7][0], zero_tol);
  CHECK_CLOSE(9.0, dfr.data[8][0], zero_tol);
}

TEST(DataFileReader_ReadAFileThatDoesNotContainAnyNumbers)
{
  std::string file_name = path + "NoNumbers.txt";
  Unfit::DataFileReader<int> dfr;
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(4u, rc);
  CHECK(dfr.data.empty());
}

TEST(DataFileReader_ReadAnEmptyFile)
{
  std::string file_name = path + "EmptyFile.txt";
  Unfit::DataFileReader<int> dfr;
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(4u, rc);
  CHECK(dfr.data.empty());
}

TEST(DataFileReader_ReadAFileThatContainsOnlyBlankLinesAndSpaces)
{
  std::string file_name = path + "SpacesAndLines.txt";
  Unfit::DataFileReader<int> dfr;
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(4u, rc);
  CHECK(dfr.data.empty());
}

TEST(DataFileReader_ReadADataFileWithAHeaderRow)
{
  std::string file_name = path + "HeaderRow.txt";
  Unfit::DataFileReader<double> dfr;
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(3u, dfr.data.size());
  CHECK_EQUAL(3u, dfr.data[0].size());
  CHECK_EQUAL(3u, dfr.data[1].size());
  CHECK_EQUAL(3u, dfr.data[2].size());
  CHECK_CLOSE(0.0, dfr.data[0][0], zero_tol);
  CHECK_CLOSE(-80.0, dfr.data[0][1], zero_tol);
  CHECK_CLOSE(0.0, dfr.data[0][2], zero_tol);
  CHECK_CLOSE(1.0, dfr.data[1][0], zero_tol);
  CHECK_CLOSE(-75.0, dfr.data[1][1], zero_tol);
  CHECK_CLOSE(1.0, dfr.data[1][2], zero_tol);
  CHECK_CLOSE(2.0, dfr.data[2][0], zero_tol);
  CHECK_CLOSE(-70.0, dfr.data[2][1], zero_tol);
  CHECK_CLOSE(5.0, dfr.data[2][2], zero_tol);
}

TEST(DataFileReader_ReadAFileThatContainsNumbersAndCharacters)
{
  std::string file_name = path + "NumbersAndCharacters.txt";
  Unfit::DataFileReader<int> dfr;
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(4u, dfr.data.size());
  CHECK_EQUAL(2u, dfr.data[0].size());
  CHECK_EQUAL(2u, dfr.data[1].size());
  CHECK_EQUAL(2u, dfr.data[2].size());
  CHECK_EQUAL(2u, dfr.data[3].size());
  CHECK_EQUAL(123, dfr.data[0][0]);
  CHECK_EQUAL(4, dfr.data[0][1]);
  CHECK_EQUAL(5, dfr.data[1][0]);
  CHECK_EQUAL(6, dfr.data[1][1]);
  CHECK_EQUAL(7, dfr.data[2][0]);
  CHECK_EQUAL(8, dfr.data[2][1]);
  CHECK_EQUAL(9, dfr.data[3][0]);
  CHECK_EQUAL(0, dfr.data[3][1]);
}

TEST(DataFileReader_ReadAFileThatDoesNotExist)  // Fails with return code 1
{
  std::string file_name = path + "ThisFileDoesNotExist.txt";
  Unfit::DataFileReader<int> dfr;
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(1u, rc);
  CHECK(dfr.data.empty());
}

TEST(DataFileReader_ReadAFileWithIrregularRows)
{
  std::string file_name = path + "IrregularRows.txt";
  Unfit::DataFileReader<int> dfr;
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(4u, dfr.data.size());
  CHECK_EQUAL(3u, dfr.data[0].size());
  CHECK_EQUAL(1u, dfr.data[1].size());
  CHECK_EQUAL(4u, dfr.data[2].size());
  CHECK_EQUAL(2u, dfr.data[3].size());
  CHECK_EQUAL(1, dfr.data[0][0]);
  CHECK_EQUAL(2, dfr.data[0][1]);
  CHECK_EQUAL(3, dfr.data[0][2]);
  CHECK_EQUAL(4, dfr.data[1][0]);
  CHECK_EQUAL(5, dfr.data[2][0]);
  CHECK_EQUAL(6, dfr.data[2][1]);
  CHECK_EQUAL(7, dfr.data[2][2]);
  CHECK_EQUAL(8, dfr.data[2][3]);
  CHECK_EQUAL(9, dfr.data[3][0]);
  CHECK_EQUAL(0, dfr.data[3][1]);
}

TEST(DataFileReader_ReadAFileWithAlternateBlankRows)
{
  std::string file_name = path + "AlternateBlankRows.txt";
  Unfit::DataFileReader<int> dfr;
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(4u, dfr.data.size());
  CHECK_EQUAL(3u, dfr.data[0].size());
  CHECK_EQUAL(1u, dfr.data[1].size());
  CHECK_EQUAL(4u, dfr.data[2].size());
  CHECK_EQUAL(2u, dfr.data[3].size());
  CHECK_EQUAL(1, dfr.data[0][0]);
  CHECK_EQUAL(2, dfr.data[0][1]);
  CHECK_EQUAL(3, dfr.data[0][2]);
  CHECK_EQUAL(4, dfr.data[1][0]);
  CHECK_EQUAL(5, dfr.data[2][0]);
  CHECK_EQUAL(6, dfr.data[2][1]);
  CHECK_EQUAL(7, dfr.data[2][2]);
  CHECK_EQUAL(8, dfr.data[2][3]);
  CHECK_EQUAL(9, dfr.data[3][0]);
  CHECK_EQUAL(0, dfr.data[3][1]);
}

TEST(DataFileReader_ReadDoubleDataIntoInt)  // Truncates
{
  std::string file_name = path + "DoubleData.txt";
  Unfit::DataFileReader<int> dfr;
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(3u, dfr.data.size());
  CHECK_EQUAL(3u, dfr.data[0].size());
  CHECK_EQUAL(3u, dfr.data[1].size());
  CHECK_EQUAL(3u, dfr.data[2].size());
  CHECK_EQUAL(1, dfr.data[0][0]);
  CHECK_EQUAL(2, dfr.data[0][1]);
  CHECK_EQUAL(3, dfr.data[0][2]);
  CHECK_EQUAL(4, dfr.data[1][0]);
  CHECK_EQUAL(5, dfr.data[1][1]);
  CHECK_EQUAL(6, dfr.data[1][2]);
  CHECK_EQUAL(6, dfr.data[2][0]);
  CHECK_EQUAL(8, dfr.data[2][1]);
  CHECK_EQUAL(9, dfr.data[2][2]);
}

TEST(DataFileReader_ReadSignedIntoUnsigned)  // Fails with return code 3
{
  std::string file_name = path + "SignedData.txt";
  Unfit::DataFileReader<unsigned> dfr;
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(3u, rc);
  CHECK(dfr.data.empty());
}

TEST(DataFileReader_ReadLargeIntegersIntoShort)  // Fails with return code 2
{
  std::string file_name = path + "LargeIntegers.txt";
  Unfit::DataFileReader<short> dfr;
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(2u, rc);
  CHECK(dfr.data.empty());
}

TEST(DataFileReader_ReadLargeNegativeIntegersIntoShort)  // Fails with code 3
{
  std::string file_name = path + "LargeNegativeIntegers.txt";
  Unfit::DataFileReader<short> dfr;
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(3u, rc);
  CHECK(dfr.data.empty());
}

TEST(DataFileReader_ReadLargeDoublesIntoFloat)  // Fails with return code 2
{
  std::string file_name = path + "LargeFloats.txt";
  Unfit::DataFileReader<float> dfr;
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(2u, rc);
  CHECK(dfr.data.empty());
}

TEST(DataFileReader_ReadLargeNegativeDoublesIntoFloat)  // Fails with code 3
{
  std::string file_name = path + "LargeNegativeFloats.txt";
  Unfit::DataFileReader<float> dfr;
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(3u, rc);
  CHECK(dfr.data.empty());
}

// *** Tests for the AddDelimiters method ***

TEST(DataFileReader_ReadNumbersWithAlternateDelimiters)
{
  std::string file_name = path + "AlternateDelimiters.txt";
  Unfit::DataFileReader<int> dfr;
  dfr.AddDelimiters("abcde");
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(2u, dfr.data.size());
  CHECK_EQUAL(4u, dfr.data[0].size());
  CHECK_EQUAL(4u, dfr.data[1].size());
  CHECK_EQUAL(1, dfr.data[0][0]);
  CHECK_EQUAL(2, dfr.data[0][1]);
  CHECK_EQUAL(3, dfr.data[0][2]);
  CHECK_EQUAL(4, dfr.data[0][3]);
  CHECK_EQUAL(5, dfr.data[1][0]);
  CHECK_EQUAL(6, dfr.data[1][1]);
  CHECK_EQUAL(7, dfr.data[1][2]);
  CHECK_EQUAL(8, dfr.data[1][3]);
}

// *** Tests for the RetrieveColumn method ***

TEST(DataFileReader_RetrieveColumnsFromRegularRowsOfIntegers)
{
  std::string file_name = path + "RegularRows.txt";
  Unfit::DataFileReader<int> dfr;
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(3u, dfr.data.size());
  CHECK_EQUAL(3u, dfr.data[0].size());
  CHECK_EQUAL(3u, dfr.data[1].size());
  CHECK_EQUAL(3u, dfr.data[2].size());

  std::vector<int> column;
  rc = dfr.RetrieveColumn(0, column);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(1, column[0]);
  CHECK_EQUAL(4, column[1]);
  CHECK_EQUAL(7, column[2]);
  rc = dfr.RetrieveColumn(1, column);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(2, column[0]);
  CHECK_EQUAL(5, column[1]);
  CHECK_EQUAL(8, column[2]);
  rc = dfr.RetrieveColumn(2, column);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(3, column[0]);
  CHECK_EQUAL(6, column[1]);
  CHECK_EQUAL(9, column[2]);
}

TEST(DataFileReader_RetrieveColumnsFromRegularRowsOfDoubles)
{
  std::string file_name = path + "RegularRows.txt";
  Unfit::DataFileReader<double> dfr;
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(3u, dfr.data.size());
  CHECK_EQUAL(3u, dfr.data[0].size());
  CHECK_EQUAL(3u, dfr.data[1].size());
  CHECK_EQUAL(3u, dfr.data[2].size());

  std::vector<double> column;
  rc = dfr.RetrieveColumn(0, column);
  CHECK_EQUAL(0u, rc);
  CHECK_CLOSE(1.0, column[0], zero_tol);
  CHECK_CLOSE(4.0, column[1], zero_tol);
  CHECK_CLOSE(7.0, column[2], zero_tol);
  rc = dfr.RetrieveColumn(1, column);
  CHECK_EQUAL(0u, rc);
  CHECK_CLOSE(2.0, column[0], zero_tol);
  CHECK_CLOSE(5.0, column[1], zero_tol);
  CHECK_CLOSE(8.0, column[2], zero_tol);
  rc = dfr.RetrieveColumn(2, column);
  CHECK_EQUAL(0u, rc);
  CHECK_CLOSE(3.0, column[0], zero_tol);
  CHECK_CLOSE(6.0, column[1], zero_tol);
  CHECK_CLOSE(9.0, column[2], zero_tol);
}

TEST(DataFileReader_RetrieveColumnsFromIrregularRows)
{
  std::string file_name = path + "MoreIrregularRows.txt";
  Unfit::DataFileReader<int> dfr;
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(4u, dfr.data.size());
  CHECK_EQUAL(4u, dfr.data[0].size());
  CHECK_EQUAL(5u, dfr.data[1].size());
  CHECK_EQUAL(4u, dfr.data[2].size());
  CHECK_EQUAL(3u, dfr.data[3].size());

  std::vector<int> column;
  rc = dfr.RetrieveColumn(0, column);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(1, column[0]);
  CHECK_EQUAL(4, column[1]);
  CHECK_EQUAL(5, column[2]);
  CHECK_EQUAL(9, column[3]);
  rc = dfr.RetrieveColumn(1, column);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(2, column[0]);
  CHECK_EQUAL(9, column[1]);
  CHECK_EQUAL(6, column[2]);
  CHECK_EQUAL(0, column[3]);
  rc = dfr.RetrieveColumn(2, column);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(3, column[0]);
  CHECK_EQUAL(1, column[1]);
  CHECK_EQUAL(7, column[2]);
  CHECK_EQUAL(3, column[3]);
  rc = dfr.RetrieveColumn(3, column);  // Incomplete, fails with rc=3
  CHECK_EQUAL(3u, rc);
  CHECK(column.empty());
  rc = dfr.RetrieveColumn(4, column);  // Incomplete, fails with rc=3
  CHECK_EQUAL(3u, rc);
  CHECK(column.empty());
  rc = dfr.RetrieveColumn(5, column);  // Does not exist, fails with rc=2
  CHECK_EQUAL(2u, rc);
  CHECK(column.empty());
}

TEST(DataFileReader_RetrieveIncompleteColumns)
{
  std::string file_name = path + "MoreIrregularRows.txt";
  Unfit::DataFileReader<int> dfr;
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(4u, dfr.data.size());
  CHECK_EQUAL(4u, dfr.data[0].size());
  CHECK_EQUAL(5u, dfr.data[1].size());
  CHECK_EQUAL(4u, dfr.data[2].size());
  CHECK_EQUAL(3u, dfr.data[3].size());

  std::vector<int> column;
  rc = dfr.RetrieveColumn(0, column);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(1, column[0]);
  CHECK_EQUAL(4, column[1]);
  CHECK_EQUAL(5, column[2]);
  CHECK_EQUAL(9, column[3]);
  rc = dfr.RetrieveColumn(1, column);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(2, column[0]);
  CHECK_EQUAL(9, column[1]);
  CHECK_EQUAL(6, column[2]);
  CHECK_EQUAL(0, column[3]);
  rc = dfr.RetrieveColumn(2, column);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(3, column[0]);
  CHECK_EQUAL(1, column[1]);
  CHECK_EQUAL(7, column[2]);
  CHECK_EQUAL(3, column[3]);
  rc = dfr.RetrieveColumn(3, column, true);  // Incomplete, rc=3
  CHECK_EQUAL(3u, rc);
  CHECK_EQUAL(3u, column.size());
  CHECK_EQUAL(4, column[0]);
  CHECK_EQUAL(8, column[1]);
  CHECK_EQUAL(8, column[2]);
  rc = dfr.RetrieveColumn(4, column, true);  // Incomplete, rc=3
  CHECK_EQUAL(3u, rc);
  CHECK_EQUAL(1u, column.size());
  CHECK_EQUAL(5, column[0]);
  rc = dfr.RetrieveColumn(5, column, true);  // Does not exist, fails with rc=2
  CHECK_EQUAL(2u, rc);
  CHECK(column.empty());
}

TEST(DataFileReader_RetrieveColumnsIfDataIsEmpty)  // No data, fails with rc=1
{
  Unfit::DataFileReader<int> dfr;
  std::vector<int> column;
  auto rc = dfr.RetrieveColumn(0, column);
  CHECK_EQUAL(1u, rc);
  CHECK(column.empty());
}

// *** Tests for optional skip line in ReadFile ***

TEST(DataFileReader_RetrieveColumnsFromRegularRowsOfIntegersOptional)
{
  std::string file_name = path + "RegularRows.txt";
  Unfit::DataFileReader<int> dfr;
  auto rc = dfr.ReadFile(file_name, 1);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(2u, dfr.data.size());
  CHECK_EQUAL(3u, dfr.data[0].size());
  CHECK_EQUAL(3u, dfr.data[1].size());

  std::vector<int> column;
  rc = dfr.RetrieveColumn(0, column);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(4, column[0]);
  CHECK_EQUAL(7, column[1]);
  rc = dfr.RetrieveColumn(1, column);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(5, column[0]);
  CHECK_EQUAL(8, column[1]);
}

TEST(DataFileReader_RetrieveColumnsFromRegularRowsOfDoublesOptional)
{
  std::string file_name = path + "RegularRows.txt";
  Unfit::DataFileReader<double> dfr;
  auto rc = dfr.ReadFile(file_name, 2);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(1u, dfr.data.size());
  CHECK_EQUAL(3u, dfr.data[0].size());

  std::vector<double> column;
  rc = dfr.RetrieveColumn(0, column);
  CHECK_EQUAL(0u, rc);
  CHECK_CLOSE(7.0, column[0], zero_tol);
}

TEST(DataFileReader_RetrieveColumnsFromIrregularRowsOptional)
{
  std::string file_name = path + "MoreIrregularRows.txt";
  Unfit::DataFileReader<int> dfr;
  auto rc = dfr.ReadFile(file_name, 1);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(3u, dfr.data.size());
  CHECK_EQUAL(5u, dfr.data[0].size());
  CHECK_EQUAL(4u, dfr.data[1].size());
  CHECK_EQUAL(3u, dfr.data[2].size());

  std::vector<int> column;
  rc = dfr.RetrieveColumn(1, column);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(9, column[0]);
  CHECK_EQUAL(6, column[1]);
  CHECK_EQUAL(0, column[2]);
  rc = dfr.RetrieveColumn(2, column);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(1, column[0]);
  CHECK_EQUAL(7, column[1]);
  CHECK_EQUAL(3, column[2]);
  rc = dfr.RetrieveColumn(3, column);  // Incomplete, fails with rc=3
  CHECK_EQUAL(3u, rc);
  CHECK(column.empty());
  rc = dfr.RetrieveColumn(4, column);  // Incomplete, fails with rc=3
  CHECK_EQUAL(3u, rc);
  CHECK(column.empty());
  rc = dfr.RetrieveColumn(5, column);  // Does not exist, fails with rc=2
  CHECK_EQUAL(2u, rc);
  CHECK(column.empty());
}

TEST(DataFileReader_RetrieveColumnsFromRegularRowsOfIntegersUnreasonableSkip)
{
  std::string file_name = path + "RegularRows.txt";
  Unfit::DataFileReader<int> dfr;
  auto rc = dfr.ReadFile(file_name, 10);
  CHECK_EQUAL(4u, rc);
  CHECK_EQUAL(0u, dfr.data.size());
  rc = dfr.ReadFile(file_name, 4);
  CHECK_EQUAL(4u, rc);
  CHECK_EQUAL(0u, dfr.data.size());
}

// *** Tests for retrieving the data as a vector ***

TEST(DataFileReader_RetrieveDataRowWiseAsVectorAndReuseDataFileReader)
{
  Unfit::DataFileReader<int> dfr;

  std::string file_name = path + "RegularRows.txt";
  auto rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(3u, dfr.data.size());
  CHECK_EQUAL(3u, dfr.data[0].size());
  std::vector<int> vec = dfr.RetrieveDataRowWiseAsVector();
  CHECK_EQUAL(9u, vec.size());
  CHECK_EQUAL(1, vec[0]);
  CHECK_EQUAL(2, vec[1]);
  CHECK_EQUAL(3, vec[2]);
  CHECK_EQUAL(4, vec[3]);
  CHECK_EQUAL(5, vec[4]);
  CHECK_EQUAL(6, vec[5]);
  CHECK_EQUAL(7, vec[6]);
  CHECK_EQUAL(8, vec[7]);
  CHECK_EQUAL(9, vec[8]);

  file_name = path + "NumbersAndCharacters.txt";
  rc = dfr.ReadFile(file_name);
  CHECK_EQUAL(0u, rc);
  CHECK_EQUAL(4u, dfr.data.size());
  CHECK_EQUAL(2u, dfr.data[0].size());
  vec = dfr.RetrieveDataRowWiseAsVector();
  CHECK_EQUAL(8u, vec.size());
  CHECK_EQUAL(123, vec[0]);
  CHECK_EQUAL(4, vec[1]);
  CHECK_EQUAL(5, vec[2]);
  CHECK_EQUAL(6, vec[3]);
  CHECK_EQUAL(7, vec[4]);
  CHECK_EQUAL(8, vec[5]);
  CHECK_EQUAL(9, vec[6]);
  CHECK_EQUAL(0, vec[7]);
}
}  // suite UnitTestDataFileReader
}  // namespace UnitTests
}  // namespace Unfit
