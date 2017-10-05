/*
 *  realmatrix.cpp
 *
 *   Created on: 22-dec-2014
 *       Author: M. El-Kebir
 */

#include "realmatrix.h"

RealMatrix::RealMatrix()
  : Matrix()
  , _C()
{
}
  
RealMatrix::RealMatrix(int m, int n)
  : Matrix(m, n)
  , _C(m, DoubleVector(n, 0))
{
}
  
std::ostream& operator<<(std::ostream& out,
                         const RealMatrix& matrix)
{
  out << matrix._C;
  if (matrix._labelsSet)
  {
    out << std::endl;
    
    for (int i = 0; i < matrix._m; ++i)
    {
      out << matrix._rowLabel[i] << " ";
    }
    out << std::endl;
    
    for (int j = 0; j < matrix._n; ++j)
    {
      out << matrix._colLabel[j] << " ";
    }
    out << std::endl;
  }
  return out;
}
  
std::istream& operator>>(std::istream& in,
                         RealMatrix& matrix)
{
  in >> matrix._C;
  
  matrix._m = matrix._C.size();
  matrix._n = matrix._C.empty() ? 0 : matrix._C.front().size();
  matrix._rowLabel.resize(matrix._m);
  matrix._colLabel.resize(matrix._n);
  
  return in;
}
  
RealMatrix RealMatrix::subMatrix(const IntVector& columns) const
{
  int n = columns.size();
  
  RealMatrix res(_m, n);
  res._rowLabel = _rowLabel;
  
  for (int j = 0; j < n; ++j)
  {
    int jj = columns[j];
    res._colLabel[j] = _colLabel[jj];
    for (int i = 0; i < _m; ++i)
    {
      res._C[i][j] = _C[i][jj];
    }
  }
  
  res._labelsSet = _labelsSet;
  
  return res;
}
