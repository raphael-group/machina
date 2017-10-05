/*
 *  matrix.cpp
 *
 *   Created on: 22-dec-2014
 *       Author: M. El-Kebir
 */

#include "matrix.h"
#include "utils.h"
#include <sstream>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

Matrix::Matrix()
  : _m(0)
  , _n(0)
  , _rowLabel()
  , _colLabel()
  , _labelsSet(false)
{
}

Matrix::Matrix(int m, int n)
  : _m(m)
  , _n(n)
  , _rowLabel(m, "")
  , _colLabel(n, "")
  , _labelsSet(false)
{
}

void Matrix::setRowLabels(const std::string& labels)
{
  _labelsSet = true;
  std::stringstream ss(labels);
  for (int i = 0; i < _m; ++i)
  {
    ss >> _rowLabel[i];
  }
}
  
void Matrix::setColLabels(const std::string& labels)
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
  
void Matrix::setLabels(std::istream& in)
{
  _labelsSet = true;
  std::string line;
  getline(in, line);
  
  getline(in, line);
  setRowLabels(line);
  
  getline(in, line);
  setColLabels(line);
}
