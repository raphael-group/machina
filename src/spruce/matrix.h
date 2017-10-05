/*
 *  matrix.h
 *
 *   Created on: 22-dec-2014
 *       Author: M. El-Kebir
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <string>
#include <cassert>

class Matrix
{
public:
  typedef std::vector<std::string> StringVector;
  
  Matrix();
  
  Matrix(int m, int n);
  
  int getNrRows() const
  {
    return _m;
  }
  
  int getNrCols() const
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
  int _m;
  int _n;
  StringVector _rowLabel;
  StringVector _colLabel;
  bool _labelsSet;
  
private:
  void setRowLabels(const std::string& labels);
  void setColLabels(const std::string& labels);
};

#endif // MATRIX_H
