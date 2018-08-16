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
#ifndef UNFIT_INCLUDE_GENERICMODELCOSTFUNCTION_HPP_
#define UNFIT_INCLUDE_GENERICMODELCOSTFUNCTION_HPP_

#include <vector>
#include "GenericModel.hpp"
#include "GenericCostFunction.hpp"

namespace Unfit
{
/**
 * This class acts as a cost function for any model, i.e., a model that
 * derives from the GenericModel class. For example, you can pass your model
 * and data to the constructor and a cost function will be created. You can
 * then, for example, pass the result straight to any of the Unfit optimisers to
 * find your model parameters. Of course you could also write your own cost
 * function, or embed your model in a cost function, but why not let this class
 * do the work for you. In terms of maths, if we write something like
 * y = f(x0, x1, x2, ... ,xN), your model will implement f(x0, x1, x2, ... ,xN).
 * Here we assume you have some data, y, x0, x1, ..., and you want to see if
 * your model is a good fit to your y data. To do this we introduce a residual.
 * For every data point, i, we will calculate the residual r_i = y_i - f(x_i).
 * The cost of the model is the vector of these residuals.
 */
class GenericModelCostFunction : public Unfit::GenericCostFunction
{
 public:
  /**
   * We require the model to be passed in when we create one of these cost
   * functions, so we must disallow use of the default constructor.
   */
  GenericModelCostFunction() = delete;

  /**
   * Create a cost function for a GenericModel. This requires that we pass in
   * the model itself (model). If you use this constructor, make sure you pass
   * the data to the cost function using the SetData function. Without data, it
   * is not possible to calculate the cost of the model.
   *
   * \param model Your model, made by deriving from GenericModel
   */
  GenericModelCostFunction(GenericModel &model);

  /**
   * Create a cost function for a GenericModel. Here we pass in the model
   * itself (model), and the data that we want to use to evaluate the model. For
   * any model, the data takes the form of one or more independent variables
   * (x0, x1, ... , xN) and a dependent variable (y). These could represent any
   * type of data you like, e.g. strain and stress, or time and distance. Note
   * that the data is not checked for consistency during construction, but it is
   * checked before the data is used. You need to make sure that your data
   * vectors are the same length (i.e., y, x0, x1, ... are the same length).
   * If they are not, you will get an error when you try to use the data with
   * your model.
   *
   * \param model Your model, made by deriving from GenericModel
   * \param x A vector of vectors containing the independent variable values
   *          for your data
   * \param y A vector containing the dependent variable values for your data
   */
  GenericModelCostFunction(GenericModel &model,
      const std::vector<std::vector<double>> &x, const std::vector<double> &y);

  /**
   * In case anyone wants to derive from this class in the future the destructor
   * the should be virtual. The default destructor is fine.
   */
  virtual ~GenericModelCostFunction() = default;

  /**
   * For any model that follows the interface prescribed by GenericModel, this
   * method takes in a vector of model parameters (c), evaluates the model using
   * these parameters, and calculates the vector of residuals between the model
   * output and the stored data. For each data point, i, the residual is
   * calculated as r_i = y_i - f(x_i), where f(x_i) is the model evaluated at
   * the point x_i with the provided parameters.
   *
   * \param c A vector of model parameters/constants
   * \return A vector containing the residuals (costs) between the model and the
   *           data
   */
  std::vector<double> operator()(const std::vector<double> &c);

  /**
   * This method takes in a set of data (independent = x, dependent = y) and
   * checks that the data is valid before we try to use it. It checks that none
   * of the vectors are empty and that they are all the same length (i.e., y,
   * x0, x1, ... are the same length). If everything is okay then it will
   * return true. If any issues are found it will return false.
   *
   * \param x A vector of vectors containing the independent variable values
   * \param y A vector containing the dependent variable values
   * \return true if the data is valid, otherwise false
   */
  bool CheckData(const std::vector<std::vector<double>> &x,
      const std::vector<double> &y);

  /**
   * A simple method to extract a copy of the data that is currently stored
   * in the cost function.
   *
   * \param x A vector of vectors containing the independent variable values
   * \param y A vector containing the dependent variable values
   */
  void GetData(std::vector<std::vector<double>> &x, std::vector<double> &y);

  /**
   * This method takes in a set of data (independent = x, dependent = y) that
   * then overwrites any data that is stored within the cost function. However,
   * before anything is overwritten, CheckData is called to make sure the new
   * data is valid. If the data vectors are empty, or are different lengths then
   * the method will return false and nothing will be overwritten.
   *
   * \param x A vector of vectors containing the independent variable values
   * \param y A vector containing the dependent variable values
   * \return true if the data was successfully stored, otherwise false (in which
   *           case the method did nothing)
   */
  bool SetData(const std::vector<std::vector<double>> &x,
      const std::vector<double> &y);

 private:
  /** Store a reference to the model. */
  GenericModel &model_;
  /** Store the independent variable data so we don't pass it in each time. */
  std::vector<std::vector<double>> x_;
  /** Store the dependent variable data so we don't pass it in each time. */
  std::vector<double> y_;
  /** Store if valid data has been successfully added to the cost function */
  bool has_data_;
};

}  // namespace Unfit

#endif  // UNFIT_INCLUDE_GENERICMODELCOSTFUNCTION_HPP_
