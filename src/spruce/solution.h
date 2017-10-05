/*
 *  solution.h
 *
 *   Created on: 1-oct-2015
 *       Author: M. El-Kebir
 */

#ifndef SOLUTION_H
#define SOLUTION_H

#include "utils.h"
#include "perfectphylomatrix.h"
#include "realmatrix.h"
#include "realtensor.h"
#include "statetree.h"

namespace gm {
  
// forward class declaration
class RootedCladisticEnumeration;

class Solution
{
public:
  typedef std::vector<StateTree> StateTreeVector;
  
  Solution(const RealTensor& F,
           const StateTreeVector& S,
           int m, int n, int k);
  
  Solution();
  
  const PerfectPhyloMatrix& A() const
  {
    return _A;
  }
  
  PerfectPhyloMatrix& A()
  {
    return _A;
  }
  
  const RealMatrix& U() const
  {
    return _U;
  }
  
  RealMatrix& U()
  {
    return _U;
  }
  
  const RealTensor& observedF() const
  {
    return _observedF;
  }
  
  const StateTreeVector& S() const
  {
    return _S;
  }
  
  RealTensor& observedF()
  {
    return _observedF;
  }
  
  StateTreeVector& S()
  {
    return _S;
  }
  
  const RealTensor& inferredF() const
  {
    return _inferredF;
  }
  
  RealTensor& inferredF()
  {
    return _inferredF;
  }
  
  double distance() const
  {
    return _distance;
  }
  
  double& distance()
  {
    return _distance;
  }
  
  friend std::ostream& operator<<(std::ostream& out, const Solution& sol);
  friend std::istream& operator>>(std::istream& in, Solution& sol);
  
  friend class RootedCladisticEnumeration;
  friend class RootedCladisticNoisyEnumeration;
  
private:
  RealTensor _observedF;
  StateTreeVector _S;
  PerfectPhyloMatrix _A;
  RealMatrix _U;
  RealTensor _inferredF;
  double _distance;
};
  
std::ostream& operator<<(std::ostream& out, const Solution& sol);
std::istream& operator>>(std::istream& in, Solution& sol);
  
} // namespace gm

#endif // SOLUTION_H