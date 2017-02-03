/*
 * Plot.hpp
 *
 *  Created on: 2 sept. 2016
 *      Author: tajchman
 */

#ifndef PLOT_HPP_
#define PLOT_HPP_

#include "values.hpp"
#include "Poisson_Parameters.hpp"

struct PlotInternal;

class Plot {
public:
  Plot();
  virtual ~Plot();

  void init(size_t n, size_t m, size_t p, Poisson_Parameters & P);

  void process(const Values &v);
  void save(const Values &v, double time, int order);
  
  void terminate();

protected:

  void save();
  PlotInternal * m_plot;
  std::string m_resultPath;
  int m_order;
};

#endif /* PLOT_HPP_ */
