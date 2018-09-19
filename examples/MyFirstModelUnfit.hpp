#ifndef MYFIRSTMODELUNFIT_HPP_
#define MYFIRSTMODELUNFIT_HPP_

#include <cmath>
#include <vector>
#include "GenericModel.hpp"

class MyFirstModelUnfit : public Unfit::GenericModel
{
 public:
  std::vector<double> operator()(const std::vector<double> &c,
      const std::vector<std::vector<double>> &t)
  {
    auto model = t[0];
    for (auto &m : model){
      m = c[0]*exp(-c[1]*m);
    }

//    // Here is an explicit way of writing the same thing (inefficient)
//    std::vector<double> times = t[0];
//    std::vector<double> model(times.size(), 0.0);
//    for (auto i = 0u; i < times.size(); ++i) {
//      model[i] = c[0]*exp(-c[1]*times[i]);
//    }

    return model;
  }
};

#endif  // MYFIRSTMODELUNFIT_HPP_
