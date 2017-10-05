/*
 *  perfectphylotree.h
 *
 *   Created on: 25-sep-2015
 *       Author: M. El-Kebir
 */

#ifndef PERFECTPHYLOTREE_H
#define PERFECTPHYLOTREE_H

#include "utils.h"
#include "perfectphylomatrix.h"
#include "perfectphylograph.h"
#include "statetree.h"
#include "realtensor.h"
#include "realmatrix.h"

namespace gm {

class PerfectPhyloTree
{
public:
  DIGRAPH_TYPEDEFS(Digraph);
  typedef Digraph::NodeMap<IntPair> IntPairNodeMap;
  typedef std::vector<Node> NodeVector;
  typedef std::vector<NodeVector> NodeMatrix;
  
  typedef std::map<std::string, std::string> StringToStringMap;
  
  typedef std::vector<StateTree> StateTreeVector;
  
  PerfectPhyloTree(const PerfectPhyloMatrix& A,
                   const StateTreeVector& S);
  
  PerfectPhyloTree(const PerfectPhyloTree& other);
  
  const PerfectPhyloMatrix& A() const
  {
    return _A;
  }
  
  const StateTree& S(int c) const
  {
    assert(0 <= c && c < _A.n());
    return _S[c];
  }
  
  int numOfVertices() const
  {
    return lemon::countNodes(_T);
  }
  
  int numOfEdges() const
  {
    return lemon::countArcs(_T);
  }
  
  int edgeRecall(const PerfectPhyloTree& otherTree) const;

  const Digraph& T() const
  {
    return _T;
  }
  
  Node root() const
  {
    return _root;
  }
  
  Node parent(Node v_dj) const
  {
    if (v_dj == lemon::INVALID)
    {
      return lemon::INVALID;
    }
    else
    {
      InArcIt a_cidj(_T, v_dj);
      if (a_cidj == lemon::INVALID)
      {
        return lemon::INVALID;
      }
      else
      {
        return _T.source(a_cidj);
      }
    }
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
  
  bool isAncestral(Node v_ci, Node v_dj) const
  {
    if (v_ci == v_dj)
    {
      return true;
      
    }
    if (v_ci == _root)
    {
      return true;
    }
    
    while (v_dj != _root)
    {
      v_dj = _T.source(InArcIt(_T, v_dj));
      
      if (v_ci == v_dj)
      {
        return true;
      }
    }
    
    return false;
  }
  
  bool isMutationVertex(Node v_ci) const
  {
    // is this a mutation vertex?
    const IntPair& ci = nodeToCharState(v_ci);
    int pi_i = S(ci.first).parent(ci.second);
    
    if (pi_i < 0)
    {
      const std::string& label_i = S(ci.first).label(ci.second);

      int x_i = -1, y_i = -1, z_i = -1;
      sscanf(label_i.c_str(), "(%d,%d,%d)", &x_i, &y_i, &z_i);
      return z_i == 1;
    }
    else
    {
      const std::string& label_i = S(ci.first).label(ci.second);
      const std::string& label_pi_i = S(ci.first).label(pi_i);
      
      int x_i = -1, y_i = -1, z_i = -1;
      sscanf(label_i.c_str(), "(%d,%d,%d)", &x_i, &y_i, &z_i);
      
      int x_pi_i = -1, y_pi_i = -1, z_pi_i = -1;
      sscanf(label_pi_i.c_str(), "(%d,%d,%d)", &x_pi_i, &y_pi_i, &z_pi_i);
      
      return x_i == x_pi_i && y_i == y_pi_i;
    }
  }
  
  bool isIncomparable(Node v_ci, Node v_dj) const
  {
    return !isAncestral(v_ci, v_dj) && !isAncestral(v_dj, v_ci);
  }
  
  bool isStateComplete(int c) const
  {
    const int n = _A.n();
    const int k = _A.k();
    
    assert(0 <= c && c < n);
    for (int i = 0; i < k; ++i)
    {
      if (_S[c].isPresent(i) &&
          _charStateToNode[c][i] == lemon::INVALID)
      {
        return false;
      }
    }
    
    return true;
  }
  
  bool isStateComplete() const
  {
    const int n = _A.n();
    
    // todo: adapt to check only characters that are present!
    // are all states of S_c present for all characters c?
    for (int c = 0; c < n; ++c)
    {
      if (!isStateComplete(c))
      {
        return false;
      }
    }
    
    return true;
  }
  
  void writeDOT(std::ostream& out) const;
  void writeDOT(const RealTensor& F,
                const RealMatrix& U,
                const StringToStringMap& label2color,
                std::ostream& out) const;
  
  friend std::ostream& operator<<(std::ostream& out, const PerfectPhyloTree& T);
  friend std::istream& operator>>(std::istream& in, PerfectPhyloTree& T);
  
  const std::string& label(Node v_ci) const
  {
    return _label[v_ci];
  }
  
  void setLabels(const RealTensor& F);
  
private:
  PerfectPhyloMatrix _A;
  StateTreeVector _S;
  
  Digraph _T;
  Node _root;
  
  NodeMatrix _charStateToNode;
  IntPairNodeMap _nodeToCharState;
  
  StringNodeMap _label;
  
  void initArcs(const PerfectPhyloGraph& G,
                PerfectPhyloGraph::Node v_ci,
                PerfectPhyloGraph::BoolNodeMap& visited);
};
  
std::ostream& operator<<(std::ostream& out, const PerfectPhyloTree& T);
std::istream& operator>>(std::istream& in, PerfectPhyloTree& T);
  
} // namespace gm

#endif // PERFECTPHYLOTREE_H