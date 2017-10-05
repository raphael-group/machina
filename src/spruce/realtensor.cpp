/*
 *  realtensor.cpp
 *
 *   Created on: 1-jul-2015
 *       Author: M. El-Kebir
 */

#include "realtensor.h"

RealTensor::RealTensor()
  : Tensor()
  , _C()
  , _trackChanges(false)
  , _deltas()
{
}
  
RealTensor::RealTensor(int k, int m, int n)
  : Tensor(k, m, n)
  , _C(k, DoubleMatrix(m, DoubleVector(n, 0)))
  , _trackChanges(false)
  , _deltas()
{
}
  
void RealTensor::rollBack()
{
  for (DeltaListRevIt it = _deltas.rbegin(); it != _deltas.rend(); ++it)
  {
    const Delta& d = *it;
    _C[d._i][d._p][d._c] = d._value;
  }
  _deltas.clear();
}
  
std::ostream& operator<<(std::ostream& out,
                         const RealTensor& tensor)
{
  out << tensor._C;
  if (tensor._labelsSet)
  {
    out << std::endl;
    
    for (int i = 0; i < tensor._m; ++i)
    {
      out << tensor._rowLabel[i] << " ";
    }
    out << std::endl;
    
    for (int j = 0; j < tensor._n; ++j)
    {
      out << tensor._colLabel[j] << " ";
    }
    out << std::endl;
  }
  return out;
}
  
std::istream& operator>>(std::istream& in,
                         RealTensor& tensor)
{
  in >> tensor._C;
  
  tensor._k = tensor._C.size();
  tensor._m = tensor._k == 0 ? 0 : tensor._C.front().size();
  tensor._n = tensor._m == 0 ? 0 : tensor._C.front().front().size();
  tensor._rowLabel.resize(tensor._m);
  tensor._colLabel.resize(tensor._n);
  
  return in;
}
  
RealTensor RealTensor::subTensor(const IntVector& columns) const
{
  int n = columns.size();
  
  RealTensor res(_k, _m, n);
  res._rowLabel = _rowLabel;
  
  for (int c = 0; c < n; ++c)
  {
    int cc = columns[c];
    res._colLabel[c] = _colLabel[cc];
    
    for (int i = 0; i < _k; ++i)
    {
      for (int p = 0; p < _m; ++p)
      {
        res._C[i][p][c] = _C[i][p][cc];
      }
    }
  }
  
  res._labelsSet = _labelsSet;
  
  return res;
}
