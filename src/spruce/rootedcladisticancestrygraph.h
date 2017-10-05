/*
 *  rootedcladisticancestrygraph.h
 *
 *   Created on: 25-sep-2015
 *       Author: M. El-Kebir
 */

#ifndef ROOTEDCLADISTICANCESTRYGRAPH_H
#define ROOTEDCLADISTICANCESTRYGRAPH_H

#include "utils.h"
#include "realtensor.h"
#include "statetree.h"
#include <lemon/list_graph.h>
#include <lemon/adaptors.h>
#include <set>
#include <list>

namespace gm {
  
class RootedCladisticAncestryGraph
{
public:
  DIGRAPH_TYPEDEFS(Digraph);
  typedef std::vector<Node> NodeVector;
  typedef std::vector<NodeVector> NodeMatrix;
  
  typedef std::vector<StateTree> StateTreeVector;

  typedef std::set<IntPair> IntPairSet;
  typedef std::list<IntPair> IntPairList;
  
  typedef Digraph::NodeMap<IntPairSet> IntPairNodeMap;
  typedef Digraph::NodeMap<IntPairList> IntPairListNodeMap;
  
  typedef std::set<int> IntSet;
  typedef IntSet::const_iterator IntSetIt;
  typedef std::set<IntSet> IntSetFamily;
  typedef IntSetFamily::const_iterator IntSetFamilyIt;
  typedef IntSetFamily::iterator IntSetFamilyNonConstIt;
  typedef Digraph::ArcMap<IntSetFamily> IntSetFamilyArcMap;
  typedef std::pair<int, IntSetFamily> CharIntSetFamilyPair;
  typedef std::vector<IntSetFamily> IntSetFamilyVector;
  typedef Digraph::NodeMap<IntSetFamilyVector> IntSetFamilyVectorMap;
  typedef Digraph::NodeMap<CharIntSetFamilyPair> CharIntSetFamilyPairNodeMap;
  typedef std::set<Arc> ArcSet;
  typedef ArcSet::const_iterator ArcSetIt;
  typedef std::set<Node> NodeSet;
  typedef NodeSet::const_iterator NodeSetIt;
  typedef Digraph::NodeMap<IntSet> IntSetNodeMap;
  typedef lemon::SubDigraph<const Digraph> SubDigraph;
  typedef SubDigraph::ArcIt SubArcIt;
  typedef SubDigraph::NodeIt SubNodeIt;
  typedef SubDigraph::OutArcIt SubOutArcIt;
  typedef SubDigraph::InArcIt SubInArcIt;
  
  typedef std::pair<IntSet, IntSet> IntSetPair;
  typedef std::list<IntSetPair> IntSetPairList;
  typedef IntSetPairList::const_iterator IntSetPairListIt;
  typedef std::vector<IntSetPairList> IntSetPairListVector;
  typedef IntSetPairListVector::const_iterator IntSetPairListVectorIt;
  typedef std::vector<IntSetPairListVector> IntSetPairListMatrix;
  typedef IntSetPairListMatrix::const_iterator IntSetPairListMatrixIt;
  typedef Digraph::ArcMap<IntSetPairList> IntSetPairListArcMap;
    
  RootedCladisticAncestryGraph(const RealTensor& F,
                               const StateTreeVector& S);
  
  virtual void writeDOT(std::ostream& out) const;
  
  int numOfVertices() const
  {
    return lemon::countNodes(_G);
  }
  
  const RealTensor& F() const
  {
    return _F;
  }
  
  const StateTree& S(int c) const
  {
    assert(0 <= c && c < _F.n());
    return _S[c];
  }
  
  const StateTreeVector& S() const
  {
    return _S;
  }
  
  const Digraph& G() const
  {
    return _G;
  }
   
  Node root() const
  {
    return _root;
  }
  
  Node charStateToNode(int c, int i) const
  {
    assert(0 <= c && c < _F.n());
    assert(0 <= i && i < _F.k());
    
    return _charStateToNode[c][i];
  }
  
  const IntPairSet& nodeToCharState(Node v_ci) const
  {
    return _nodeToCharState[v_ci];
  }
  
  const IntPairList& nodeToCharStateList(Node v_ci) const
  {
    return _nodeToCharStateList[v_ci];
  }
  
  const std::string& label(Node v_ci) const
  {
    return _label[v_ci];
  }
  
  void setLabels(const RealTensor& F)
  {
    for (NodeIt v_ci(_G); v_ci != lemon::INVALID; ++v_ci)
    {
      if (v_ci == _root)
      {
        _label[v_ci] = "(*,0)";
      }
      else
      {
        const IntPairSet& X_ci = _nodeToCharState[v_ci];
        std::string lbl;
        
        bool first = true;
        for (auto it = X_ci.begin(); it != X_ci.end(); ++it)
        {
          const IntPair& ci = *it;
          if (first)
            first = false;
          else
            lbl += ";";
          
          lbl += "(" + F.getColLabel(ci.first) + "," + _S[ci.first].label(ci.second) + ")";
          _label[v_ci] = lbl;
        }
      }
    }
  }
  
  virtual void init();
  
protected:
  const RealTensor& _F;
  const StateTreeVector& _S;
  
  Digraph _G;
  Node _root;
  NodeMatrix _charStateToNode;
  IntPairNodeMap _nodeToCharState;
  IntPairListNodeMap _nodeToCharStateList;
  StringNodeMap _label;
  
  struct Compare
  {
  public:
    Compare(const StateTreeVector& S)
      : _S(S)
    {
    }
    
    bool operator()(const IntPair& ci, const IntPair& dj)
    {
      if (ci.first < dj.first)
      {
        return true;
      }
      else if (ci.first > dj.first)
      {
        return false;
      }
      else
      {
        assert(ci != dj);
        // is ci an ancestor of dj?
        if (_S[ci.first].D(ci.second).count(dj.second))
        {
          return true;
        }
        else
        {
          return false;
        }
      }
    }
    
  private:
    const StateTreeVector& _S;
  };
};
  
} // namespace gm

#endif // ROOTEDCLADISTICANCESTRYGRAPH_H
