#ifndef MYFIRSTMODEL_HPP_
#define MYFIRSTMODEL_HPP_

#include <cmath>
#include <vector>
#include "GenericModel.hpp"

class MyFirstModel : public Unfit::GenericModel
{
 public:
  std::vector<double> operator()(const std::vector<double> &c,
      const std::vector<std::vector<double>> &t)
  {
    std::vector<double> model(t[0].size());
    for (auto i = 0u; i < t[0].size(); ++i) {
      model[i] = c[0]*exp(-c[1]*t[0][i]);
    }
    return model;
  }
};

#endif  // MYFIRSTMODEL_HPP_
