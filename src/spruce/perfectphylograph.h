/*
 *  perfectphylograph.h
 *
 *   Created on: 25-sep-2015
 *       Author: M. El-Kebir
 */

#ifndef PERFECTPHYLOGRAPH_H
#define PERFECTPHYLOGRAPH_H

#include "utils.h"
#include "perfectphylomatrix.h"

namespace gm {
  
class PerfectPhyloGraph
{
public:
  GRAPH_TYPEDEFS(Graph);
  typedef Graph::NodeMap<IntPair> IntPairNodeMap;
  typedef std::vector<Node> NodeVector;
  typedef std::vector<NodeVector> NodeMatrix;
  
  PerfectPhyloGraph(const PerfectPhyloMatrix& A);
  
  void writeDOT(std::ostream& out) const;
  
  const Graph& G() const
  {
    return _G;
  }
  
  Node root() const
  {
    return _root;
  }
  
  Node charStateToNode(int c, int i) const
  {
    assert(0 <= c && c < _A.n());
    assert(0 <= i && i < _A.k());
    
    return _charStateToNode[c][i];
  }
  
  const IntPair& nodeToCharState(Node v_ci) const
  {
    return _nodeToCharState[v_ci];
  }
  
private:
  const PerfectPhyloMatrix& _A;
  Graph _G;
  Node _root;
  
  NodeMatrix _charStateToNode;
  IntPairNodeMap _nodeToCharState;
};

}

#endif // PERFECTPHYLOGRAPH_H