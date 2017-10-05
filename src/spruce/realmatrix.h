/*
 *  realmatrix.h
 *
 *   Created on: 22-dec-2014
 *       Author: M. El-Kebir
 */

#ifndef REALMATRIX_H
#define REALMATRIX_H

#include "matrix.h"
#include "utils.h"

class RealMatrix : public Matrix
{
public:
  RealMatrix();
  
  RealMatrix(int m, int n);
  
  const DoubleMatrix& getMatrix() const
  {
    return _C;
  }
  
  double operator()(int row, int col) const
  {
    assert(0 <= row && row < _m);
    assert(0 <= col && col < _n);
    
    return _C[row][col];
  }
  
  void set(int row, int col, double val)
  {
    assert(0 <= row && row < _m);
    assert(0 <= col && col < _n);
    
    _C[row][col] = val;
  }
  
  bool operator==(const RealMatrix& other) const
  {
    return _m == other._m && _n == other._n && _C == other._C;
  }
  
  bool operator!=(const RealMatrix& other) const
  {
    return !this->operator==(other);
  }
  
  RealMatrix subMatrix(const IntVector& columns) const;
  
  friend std::ostream& operator<<(std::ostream& out,
                                  const RealMatrix& matrix);
  
  friend std::istream& operator>>(std::istream& in,
                                  RealMatrix& matrix);
  
protected:
  DoubleMatrix _C;
};

#endif // REALMATRIX_H
