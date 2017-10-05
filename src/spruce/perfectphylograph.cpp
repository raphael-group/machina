/*
 *  perfectphylograph.cpp
 *
 *   Created on: 25-sep-2015
 *       Author: M. El-Kebir
 */

#include "perfectphylograph.h"

namespace gm {
  
PerfectPhyloGraph::PerfectPhyloGraph(const PerfectPhyloMatrix& A)
  : _A(A)
  , _G()
  , _root(lemon::INVALID)
  , _charStateToNode(_A.n(), NodeVector(_A.k(), lemon::INVALID))
  , _nodeToCharState(_G)
{
  const int n = _A.n();
  const int k = _A.k();
  
  // add root node
  _root = _G.addNode();
  _nodeToCharState[_root] = std::make_pair(0, 0);
  for (int c = 0; c < n; ++c)
  {
    _charStateToNode[c][0] = _root;
  }
  
  // add the other n(k-1) nodes
  for (int c = 0; c < n; ++c)
  {
    for (int i = 1; i < k; ++i)
    {
      if (_A.defined(c, i))
      {
        Node v_ci = _G.addNode();
        _charStateToNode[c][i] = v_ci;
        _nodeToCharState[v_ci] = std::make_pair(c, i);
      }
    }
  }
  
  // add edges
  for (NodeIt v_ci(_G); v_ci != lemon::INVALID; ++v_ci)
  {
    const IntPair& ci = _nodeToCharState[v_ci];
    for (NodeIt v_dj = v_ci; v_dj != lemon::INVALID; ++v_dj)
    {
      if (v_ci == v_dj)
        continue;
      
      const IntPair& dj = _nodeToCharState[v_dj];
      if (_A.hammingDist(ci, dj) == 1)
      {
        _G.addEdge(v_ci, v_dj);
      }
    }
  }
}
  
void PerfectPhyloGraph::writeDOT(std::ostream& out) const
{
  out << "graph G {" << std::endl;
  
  for (NodeIt v_ci(_G); v_ci != lemon::INVALID; ++v_ci)
  {
    const IntPair& ci = _nodeToCharState[v_ci];
    
    out << "\t" << _G.id(v_ci) << " [label=\"";
    if (v_ci == _root)
      out << "(*,0)";
    else
      out << "(" << ci.first << "," << ci.second << ")";
    out << "\"]" << std::endl;
  }
  
  for (EdgeIt e(_G); e != lemon::INVALID; ++e)
  {
    out << "\t" << _G.id(_G.u(e)) << " -- " << _G.id(_G.v(e)) << std::endl;
  }
  
  out << "}" << std::endl;
}
  
} // namespace gm