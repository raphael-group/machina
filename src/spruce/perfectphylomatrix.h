/*
 *  perfectphylomatrix.h
 *
 *   Created on: 25-sep-2015
 *       Author: M. El-Kebir
 */

#ifndef PERFECTPHYLOMATRIX_H
#define PERFECTPHYLOMATRIX_H

#include "utils.h"

namespace gm {
  
class PerfectPhyloMatrix
{
public:
  PerfectPhyloMatrix();
  
  PerfectPhyloMatrix(int n, int k);
  
  int n() const
  {
    return _n;
  }
  
  int k() const
  {
    return _k;
  }
  
  int operator()(const IntPair& ci, int d) const
  {
    return this->operator()(ci.first, ci.second, d);
  }
  
  int operator()(int c, int i, int d) const
  {
    assert(0 <= c && c < _n);
    assert(0 <= d && d < _n);
    assert(0 <= i && i < _k);
    
    if (_A[c][i].empty())
      return -1;
    else
      return _A[c][i][d];
  }
  
  void set(int c, int i, int d, int j)
  {
    assert(0 <= c && c < _n);
    assert(0 <= d && d < _n);
    assert(0 <= i && i < _k);
    assert(0 <= j && j < _k);
    
    // c == d => i == j
    assert( !(c == d) || (i == j));
    
    if (_A[c][i].empty())
      _A[c][i] = IntVector(_n, 0);
    
    _A[c][i][d] = j;
  }
  
  void set(const IntPair& ci, const IntPair& dj)
  {
    set(ci.first, ci.second, dj.first, dj.second);
  }
  
  void set(const IntPair& ci, int d, int j)
  {
    set(ci.first, ci.second, d, j);
  }
  
  bool defined(int c, int i) const
  {
    return !_A[c][i].empty();
  }
  
  bool operator==(const PerfectPhyloMatrix& other) const
  {
    return _n == other._n && _k == other._k && _A == other._A;
  }
  
  bool operator!=(const PerfectPhyloMatrix& other) const
  {
    return !this->operator==(other);
  }
  
  int hammingDist(const IntPair& ci, const IntPair& dj) const
  {
    return hammingDist(ci.first, ci.second, dj.first, dj.second);
  }
  
  int hammingDist(int c, int i, int d, int j) const
  {
    int dist = 0;
    for (int e = 0; e < _n; ++e)
    {
      if (this->operator()(c, i, e) != this->operator()(d, j, e))
        ++dist;
    }
    return dist;
  }
  
  friend std::ostream& operator<<(std::ostream& out, const PerfectPhyloMatrix& A);
  friend std::istream& operator>>(std::istream& in, PerfectPhyloMatrix& A);
  
private:
  typedef std::vector<IntMatrix> StlInt3Matrix;
  
  int _n;
  int _k;
  StlInt3Matrix _A;
  
  int toRowIdx(int c, int i) const
  {
    return i == 0 ? 0 : _k * (i - 1) + c + 1;
  }
};
  
std::ostream& operator<<(std::ostream& out, const PerfectPhyloMatrix& A);
std::istream& operator>>(std::istream& in, PerfectPhyloMatrix& A);
  
} // namespace gm

#endif // PERFECTPHYLOMATRIX_H
