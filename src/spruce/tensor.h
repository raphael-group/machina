/*
 *  tensor.h
 *
 *   Created on: 30-jun-2015
 *       Author: M. El-Kebir
 */

#ifndef TENSOR_H
#define TENSOR_H

#include "utils.h"

class Tensor
{
public:
  typedef std::vector<std::string> StringVector;
  
  Tensor();
  
  Tensor(int k, int m, int n);
  
  int k() const
  {
    return _k;
  }
  
  int m() const
  {
    return _m;
  }
  
  int n() const
  {
    return _n;
  }

  void setLabels(std::istream& in);
  
  const std::string& getRowLabel(int i) const
  {
    assert(0 <= i && i < _m);
    return _rowLabel[i];
  }
  
  const StringVector& getRowLabels() const
  {
    return _rowLabel;
  }
  
  const StringVector& getColLabels() const
  {
    return _colLabel;
  }
  
  const std::string& getColLabel(int j) const
  {
    assert(0 <= j && j < _n);
    return _colLabel[j];
  }
  
  void setRowLabel(int i, const std::string& label)
  {
    assert(0 <= i && i < _m);
    _rowLabel[i] = label;
    _labelsSet = true;
  }
  
  void setColLabel(int j, const std::string& label)
  {
    assert(0 <= j && j < _n);
    _colLabel[j] = label;
    _labelsSet = true;
  }
  
protected:
  int _k;
  int _m;
  int _n;
  StringVector _rowLabel;
  StringVector _colLabel;
  bool _labelsSet;
  
private:
  void setRowLabels(const std::string& labels);
  void setColLabels(const std::string& labels);
};

#endif // TENSOR_H
