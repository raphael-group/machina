/*
 *  tensor.cpp
 *
 *   Created on: 30-jun-2015
 *       Author: M. El-Kebir
 */

#include "tensor.h"
#include "utils.h"
#include <sstream>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

Tensor::Tensor()
  : _k(0)
  , _m(0)
  , _n(0)
  , _rowLabel()
  , _colLabel()
  , _labelsSet(false)
{
}

Tensor::Tensor(int k, int m, int n)
  : _k(k)
  , _m(m)
  , _n(n)
  , _rowLabel(m, "")
  , _colLabel(n, "")
  , _labelsSet(false)
{
}

void Tensor::setRowLabels(const std::string& labels)
{
  _labelsSet = true;
  std::stringstream ss(labels);
  for (int i = 0; i < _m; ++i)
  {
    ss >> _rowLabel[i];
  }
}
  
void Tensor::setColLabels(const std::string& labels)
{
  typedef std::vector<std::string> StringVector;
  
  _labelsSet = true;
  
  StringVector s;
  boost::split(s, labels, boost::is_any_of(" "));
  
  assert(s.size() >= _n);
  
  std::stringstream ss(labels);
  for (int j = 0; j < _n; ++j)
  {
    _colLabel[j] = s[j];
  }
}
  
void Tensor::setLabels(std::istream& in)
{
  _labelsSet = true;
  std::string line;
  getline(in, line);
  
  getline(in, line);
  setRowLabels(line);
  
  getline(in, line);
  setColLabels(line);
}
