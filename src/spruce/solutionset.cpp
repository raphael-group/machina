/*
 *  solutionset.cpp
 *
 *   Created on: 1-oct-2015
 *       Author: M. El-Kebir
 */

#include "solutionset.h"
#include "rootedcladisticancestrygraph.h"
#include "perfectphylotree.h"
#include <boost/unordered_set.hpp>

namespace gm {
  
SolutionSet::SolutionSet()
  : _sol()
{
}
  
void SolutionSet::sort()
{
  std::sort(_sol.begin(), _sol.end(), Compare());
}
  
void SolutionSet::add(const SolutionSet& sols)
{
  int n = sols.solutionCount();
  for (int i = 0; i < n; ++i)
  {
    add(sols.solution(i));
  }
}

void SolutionSet::add(const Solution& sol)
{
  _sol.push_back(sol);
}
  
int SolutionSet::unique()
{
  typedef std::pair<IntPair, IntPair> IntPairPair;
  typedef std::set<IntPairPair> ArcSet;
  typedef std::set<ArcSet> ArcSetSet;
  
  int res = 0;
  ArcSetSet AS;
  for (SolutionVectorNonConstIt it = _sol.begin(); it != _sol.end();)
  {
    const Solution& sol = *it;
    PerfectPhyloTree T(sol.A(), sol.S());
    const Digraph& TT = T.T();
    
    ArcSet arcs;
    for (Digraph::ArcIt a_cidj(TT); a_cidj != lemon::INVALID; ++a_cidj)
    {
      Digraph::Node v_ci = TT.source(a_cidj);
      Digraph::Node v_dj = TT.target(a_cidj);
      
      const IntPair& ci = T.nodeToCharState(v_ci);
      const IntPair& dj = T.nodeToCharState(v_dj);
      
      arcs.insert(std::make_pair(ci, dj));
    }
    
    if (AS.find(arcs) == AS.end())
    {
      AS.insert(arcs);
      ++it;
    }
    else
    {
      it = _sol.erase(it);
      ++res;
    }
  }
  return res;
}
  
std::ostream& operator<<(std::ostream& out, const SolutionSet& sols)
{
  int solCount = sols.solutionCount();
  out << solCount << " # solutions" << std::endl << std::endl;
  
  for (int idx = 0; idx < solCount; ++idx)
  {
    out << sols.solution(idx);
  }
  
  return out;
}
  
std::istream& operator>>(std::istream& in, SolutionSet& sols)
{
  int solCount = -1;
  
  std::string line;
  getline(in, line);
  std::stringstream ss(line);
  ss >> solCount;
  
  getline(in, line);
  
  if (solCount <= 0)
  {
    throw std::runtime_error(getLineNumber() + "Error: solCount should be nonnegative");
  }
  
  sols._sol.clear();
  for (int idx = 0; idx < solCount; ++idx)
  {
    Solution sol;
    in >> sol;
    sols.add(sol);
    
    getline(in, line);
  }
  
  return in;
}
  
} // namespace gm
