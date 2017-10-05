/*
 *  realtensor.h
 *
 *   Created on: 1-jul-2015
 *       Author: M. El-Kebir
 */

#ifndef REALTENSOR_H
#define REALTENSOR_H

#include "tensor.h"
#include "utils.h"
#include <list>

class RealTensor : public Tensor
{
public:
  RealTensor();
  
  RealTensor(int k, int m, int n);
  
  const DoubleTensor& getTensor() const
  {
    return _C;
  }
  
  double operator()(int slice, int row, int col) const
  {
    assert(0 <= slice && slice < _k);
    assert(0 <= row && row < _m);
    assert(0 <= col && col < _n);
    
    return _C[slice][row][col];
  }
  
  void set(int slice, int row, int col, double val)
  {
    assert(0 <= slice && slice < _k);
    assert(0 <= row && row < _m);
    assert(0 <= col && col < _n);
    
    if (_C[slice][row][col] != val)
    {
      if (_trackChanges)
      {
        Delta d;
        d._i = slice;
        d._p = row;
        d._c = col;
        d._value = _C[slice][row][col];
        
        _deltas.push_back(d);
      }
      _C[slice][row][col] = val;
    }
  }
  
  double getCumFreq(int p, int c, const IntSet& D) const
  {
    double res = 0;
    for (IntSetIt it = D.begin(); it != D.end(); ++it)
    {
      res += (*this)(*it, p, c);
    }
    return res;
  }
  
  bool operator==(const RealTensor& other) const
  {
    return _k == other._k && _m == other._m
      && _n == other._n && _C == other._C;
  }
  
  bool operator!=(const RealTensor& other) const
  {
    return !this->operator==(other);
  }
  
  RealTensor subTensor(const IntVector& columns) const;
  
  friend std::ostream& operator<<(std::ostream& out,
                                  const RealTensor& tensor);
  
  friend std::istream& operator>>(std::istream& in,
                                  RealTensor& tensor);
  
  bool getTrackChanges() const
  {
    return _trackChanges;
  }
  
  void setTrackChanges(bool trackChanges)
  {
    _trackChanges = trackChanges;
  }
  
  void rollBack();
  
protected:
  DoubleTensor _C;
  bool _trackChanges;
  
  struct Delta
  {
    int _p;
    int _c;
    int _i;
    double _value;
  };
  
  typedef std::list<Delta> DeltaList;
  typedef DeltaList::const_iterator DeltaListIt;
  typedef DeltaList::const_reverse_iterator DeltaListRevIt;
  
  DeltaList _deltas;
};

#endif // REALTENSOR_H
