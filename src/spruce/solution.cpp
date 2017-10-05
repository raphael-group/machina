/*
 *  solution.cpp
 *
 *   Created on: 1-oct-2015
 *       Author: M. El-Kebir
 */

#include "solution.h"

namespace gm {
  
Solution::Solution(const RealTensor& F,
                   const StateTreeVector& S,
                   int m, int n, int k)
  : _observedF(F)
  , _S(S)
  , _A(n, k)
  , _U(m, n*(k-1) + 1)
  , _inferredF(k, m, n)
  , _distance(0)
{
}
  
Solution::Solution()
  : _observedF()
  , _S()
  , _A()
  , _U()
  , _inferredF()
  , _distance(0)
{
}
  
std::ostream& operator<<(std::ostream& out, const Solution& sol)
{
  out << sol._observedF << std::endl;
  
  const int n = sol._S.size();
  for (int c = 0; c < n; ++c)
  {
    out << sol._S[c];
  }
  out << std::endl;
  
  out << sol._A;
  out << sol._U;
  out << sol._inferredF;
  out << "#distance = " << sol.distance() << std::endl;
  
  return out;
}
  
std::istream& operator>>(std::istream& in, Solution& sol)
{
  std::string line;
  
  in >> sol._observedF;
  sol._observedF.setLabels(in);
  
  getline(in, line);
  
  const int n = sol._observedF.n();
  const int k = sol._observedF.k();
  sol._S = Solution::StateTreeVector(n, StateTree(k));
  for (int c = 0; c < n; ++c)
  {
    in >> sol._S[c];
  }
  
  getline(in, line);

  in >> sol._A;
  in >> sol._U;
  in >> sol._inferredF;
  sol._inferredF.setLabels(in);
  
  return in;
}
  
} // namespace gm
