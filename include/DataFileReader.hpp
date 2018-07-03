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
#ifndef UNFIT_INCLUDE_DATAFILEREADER_HPP_
#define UNFIT_INCLUDE_DATAFILEREADER_HPP_

#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace Unfit
{
/**
 * \brief Reads in numeric data from a file
 *
 * The main goal of this object is to read numeric data from a file and to
 * store it in an array (a vector of vectors). See the documentation for
 * the ReadFile method for details. In addition, the string splitting method
 * that is used internally can be called directly if that functionality is
 * needed. It is also possible to access the data that has been read in and
 * extract a column via the RetrieveColumn method. The key design consideration
 * was that the implementation should be self contained, i.e., it does not use
 * e.g. boost's split function. It is also a header-only implementation,
 * meaning that you only have to include this header file in your code to
 * access all of the functionality. You do not have to link against any
 * extra libraries. The struct has been templated so you can choose to read the
 * data into a container of any numeric data type. Normal truncation rules
 * apply.
 *
 * Extensive unit testing has been performed on this implementation using the
 * UnitTest++ testing framework. The tests and example files can be found in
 * the "unittests" directory. gcov reports 100% line and function coverage,
 * and valgrind reports no memory errors or leaks for the test suite.
 */
template <typename T>
  struct DataFileReader
{
  /**
   * Create a DataFileReader with no data and initialise the default
   * element delimiters to tab, space and comma.
   *
   * Intended use:
   *   DataFileReader data_file_reader;
   */
  DataFileReader();

  /**
   * This method attempts to open the file name that is passed in and will read
   * the (numeric) data inside and put it in a 2D array (a std::vector of
   * std::vectors, called "data". The method makes no attempt to search for the
   * file so the file name will need to include the appropriate path where
   * needed.
   *
   * The lines in the file are read in one at a time and each line is then
   * stored in a vector of numbers, hence the the storage is row-wise.
   * Non-numeric data in the file is ignored, as are blank lines. This means,
   * for example, if your data file has one or more header rows, there is no
   * need to remove them. The method will also read in lines of different
   * lengths with no problem.
   *
   * By default, the data in the file will be split on tabs, spaces and commas.
   * You can add to this list via the AddDelimiter method. Checks are performed
   * to make sure the data being read fits into the data type that has been
   * specified, and exits with a non-zero return code if this is not the case.
   *
   * If you read a second file in with the same data reader, the original data
   * stored by the reader will be lost and replaced with the new data.
   *
   * Intended use:
   *   DataFileReader data_file_reader;
   *   data_file_reader.ReadFile( "file_name" );
   *
   * \param file_name a string containing the file name and path information
   * \param skip an optional parameter to skip the first "skip" lines
   * \return a return code, 0 is success
   *
   * Return codes:
   *   0 = data was read successfully
   *   1 = file could not be opened
   *   2 = number larger than the maximum for the data type
   *   3 = number larger than the maximum negative for the data type,
   *       or writing a negative into an unsigned type
   *   4 = file contains no data (numbers)
   */
  unsigned ReadFile(std::string file_name, unsigned skip = 0);

  /**
   * The default delimiters used to split each line in the input file are
   * tab, space, and comma. If you want to add to this list, call this
   * method and pass in a quoted string. Each character is treated as a
   * separate delimiter.
   *
   * Intended use:
   *   data_file_reader.AddDelimiters( ";:" );
   *
   * \param new_delimiters the additional delimiters to be used
   */
  void AddDelimiters(const std::string &new_delimiters);

  /**
   * This method allows you to get a column from the data that has been read in
   * and copies it to a std::vector, which has to be the same type as the
   * reader (int reader = int vector, double reader = double vector, etc). The
   * first argument is the index of the desired column, starting from zero. If
   * the rows that were read in had irregular lengths, some columns will only
   * be partially populated. In this case this function will flag this with a
   * return code, and will not return the column.
   *
   * Intended use:
   *   int rc = data_file_reader.RetrieveColumn( index, column );
   *
   * \param column_number the index (from zero) of the desired column
   * \param column a vector in which the requested column is placed
   * \param return_incomplete_columns When false (default) returns an empty
   *        vector if the column is incomplete. When true returns the
   *        incomplete column
   * \return a return code, 0 is success
   *
   * Return codes:
   *   0 = success
   *   1 = reader contains no data
   *   2 = requested column does not exist
   *   3 = requested column exists, but is not fully populated
   */
  unsigned RetrieveColumn(std::size_t column_number, std::vector<T> &column,
      bool return_incomplete_columns = false);

 /**
   * This method allows you to return all of the data as a single 1D vector. The
   * result will contain all of the rows appended one after the other. For
   * example, if you have:
   *
   *   [1 2 3]
   *   [7 8 9]
   *   [4 5 6]
   *
   * then the resulting vector will be [1 2 3 7 8 9 4 5 6].
   *
   * \return a vector of the data
   */
  std::vector<T> RetrieveDataRowWiseAsVector();

  /**
   * This method is used internally by the ReadFile method, so you do not need
   * to call this at all if all you want to do is read data from a file.
   * However, it can be very useful if you have a string that you want to split
   * into separate elements based on certain delimiters. Just pass the string
   * of interest and the delimiter(s) you want to split on, and this method
   * will perform the split and give you back a vector of strings. It will
   * ignore multiple sequential delimiters (e.g. three spaces in a row does not
   * create empty elements) and removes any preceding or trailing delimiters.
   * When using multiple delimiters, each character is treated as a separate
   * delimiter.
   *
   * Intended use:
   *   data_file_reader.SplitLine( in_string, "\t ,", vector_of_words );
   *
   * \param line the string to be split
   * \param delimiters a list of the delimiters on which to split the line
   * \param words a vector of strings, one for each word after the split
   */
  void SplitLine(const std::string &line, const std::string &delimiters,
      std::vector<std::string> &words);

  /** An array to store the data that is read in */
  std::vector<std::vector<T>> data;
  /** Each line will be split on these delimiters */
  std::string default_delimiters;
};


template <class T>
  DataFileReader<T>::DataFileReader() :
    data(),
    default_delimiters("\t ,")
{}


template <class T>
  unsigned DataFileReader<T>::ReadFile(std::string file_name, unsigned skip)
{
  // Just in case this is not the first time, make sure we start from scratch
  data.clear();
  // Attempt to open the file, then check it is open
  std::ifstream data_file(file_name.c_str());
  if (!data_file.is_open()) return 1;
  // Read in the data, one line at a time
  while (data_file.good()) {
    // Get a line
    std::string line;
    getline(data_file, line);
    // Split the line based on the delimiters
    std::vector<std::string> words;
    SplitLine(line, default_delimiters, words);
    // Skip the first "skip" numeric lines if requested
    if (skip > 0) {
      --skip;
      continue;
    }
    // Convert to type T then store in a vector
    std::vector<T> numbers;
    for (auto entry = 0u; entry < words.size(); ++entry) {
      std::stringstream ss;
      ss << words[entry].c_str();
      long double value;
      // Check we have a number before including it
      if (ss >> value) {
        // Check the maximum size is okay
        if (value > std::numeric_limits<T>::max()) {
          data.clear();
          return 2;
        }
        // Check the maximum negative size is okay
        // Use std::numeric_limits::lowest when available
        if (std::numeric_limits<T>::is_integer) {
          if (value < std::numeric_limits<T>::min()) {
            data.clear();
            return 3;
          }
        }
        else if (value <
            -static_cast<long double>(std::numeric_limits<T>::max())) {
          data.clear();
          return 3;
        }
        numbers.push_back(static_cast<T>(value));
      }
    }
    // If the line contained numbers, put these into a row of the data vector
    if (!numbers.empty()) data.push_back( numbers );
  }
  if (data.empty()) return 4;
  return 0;
}


template <class T>
  void DataFileReader<T>::AddDelimiters(const std::string &new_delimiters)
{
  default_delimiters += new_delimiters;
}


template <class T>
  unsigned DataFileReader<T>::RetrieveColumn(std::size_t column_number,
    std::vector<T> &column, bool return_incomplete_columns)
{
  column.clear();
  if (data.empty()) return 1;
  std::size_t min_row_length = std::numeric_limits<std::size_t>::max();
  std::size_t max_row_length = 0u;
  for (auto d : data) {
    if (d.size() < min_row_length) min_row_length = d.size();
    if (d.size() > max_row_length) max_row_length = d.size();
  }
  if (column_number >= max_row_length) return 2;  // Column does not exist
  if (column_number >= min_row_length) {  // Column is incomplete
    if (return_incomplete_columns) {
      for (auto i = 0u; i < data.size(); ++i) {
        if (data[i].size() > column_number) {
          column.push_back(data[i][column_number]);
        }
      }
    }
    return 3;
  }
  else {
    column.resize(data.size());
    for (auto i = 0u; i < data.size(); ++i) {
      column[i] = data[i][column_number];
    }
    return 0;
  }
}


template <class T>
  std::vector<T> DataFileReader<T>::RetrieveDataRowWiseAsVector()
{
  // This does a range insert for every row of data in order
  std::vector<T> result;
  for (auto d : data) {
    result.insert(end(result), begin(d), end(d));
  }
  return result;
}


template <class T>
  void DataFileReader<T>::SplitLine(const std::string &line,
    const std::string &delimiters, std::vector<std::string> &words)
{
  // Get the first word from the line
  auto word_start = line.find_first_not_of(delimiters, 0);
  auto word_end = line.find_first_of(delimiters, word_start);

  // While we have not hit the end of the string
  while (word_start != std::string::npos || word_end != std::string::npos) {
    // Add the current word to the list (vector)
    words.push_back(line.substr(word_start, (word_end-word_start)));
    // Get the next word in the line
    word_start = line.find_first_not_of(delimiters, word_end);
    word_end = line.find_first_of(delimiters, word_start);
  }
}

}  // namespace Unfit

#endif
