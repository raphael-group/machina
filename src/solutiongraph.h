/*
 * solutiongraph.h
 *
 *  Created on: 21-aug-2017
 *      Author: M. El-Kebir
 */

#ifndef SOLUTIONGRAPH_H
#define SOLUTIONGRAPH_H

#include "utils.h"
#include "clonetree.h"

class SolutionGraph
{
public:
  SolutionGraph(const CloneTree& T,
                const IntPairNodeMap& leafLabelingT,
                const IntPairToNodeMap& lcaT,
                const IntPairSetNodeMap& vertexToStateSet,
                const IntPairToNodeSetMap& X,
                int primaryIndex,
                const StringVector& indexToAnatomicalSite);
  
  ~SolutionGraph();
  
  void writeDOT(std::ostream& out) const;
  
  void writeDOT(std::ostream& out, const SubDigraph& S) const;
  
  int getNrStateTrees() const
  {
    return _stateTrees.size();
  }
  
  const SubDigraph& getStateTree(int index) const
  {
    assert(0 <= index && index < _stateTrees.size());
    return _stateTrees[index]->subT();
  }
  
  const CloneTree& getRefinedCloneTree(int index) const
  {
    assert(0 <= index && index < _refinements.size());
    return _refinements[index]->_Tprime;
  }
  
  const StringNodeMap& getVertexLabeling(int index) const
  {
    assert(0 <= index && index < _refinements.size());
    return _refinements[index]->_lPlus;
  }

private:
  void initG();
  
  void enumerate();
  
  void refine(const SubDigraph& S);
  
  void refine(const SubDigraph& S,
              Node v_inT,
              Digraph& Tprime,
              Node v_inTprime,
              StringNodeMap& label,
              StringNodeMap& lPlus);
  
  struct Refinement
  {
  public:
    Refinement(const Digraph& T,
               Node root,
               const StringNodeMap& id,
               const StringNodeMap& lPlus)
      : _Tprime(T, root, id, lPlus)
      , _lPlus(_Tprime.tree())
    {
      for (NodeIt v(T); v != lemon::INVALID; ++v)
      {
        const std::string& label_v = id[v];
        Node vv = _Tprime.getNodeByLabel(label_v);
        _lPlus[vv] = lPlus[v];
      }
    }
    
    CloneTree _Tprime;
    StringNodeMap _lPlus;
  };
  
  typedef std::vector<Refinement*> RefinementVector;
  
  struct SubTree
  {
  public:
    SubTree(const Digraph& G)
      : _G(G)
      , _filterNodes(_G, false)
      , _filterArcs(_G, false)
      , _subT(_G, _filterNodes, _filterArcs)
    {
    }
    
    SubDigraph& subT()
    {
      return _subT;
    }

  private:
    const Digraph& _G;
    BoolNodeMap _filterNodes;
    BoolArcMap _filterArcs;
    SubDigraph _subT;
  };
  
  typedef std::vector<SubTree*> SubTreeVector;
  
  bool isConnected(const SubDigraph& S,
                   const IntPairSet& Sigma_u,
                   IntPair& root_sc) const;
  
  bool areSiblings(const SubDigraph& S,
                   const IntPairSet& Sigma_u,
                   IntPair& parent_sc) const;
  
private:
  const CloneTree& _T;
  const IntPairNodeMap& _leafLabelingT;
  IntPairToNodeMap _lcaT;
  const IntPairSetNodeMap& _vertexToStateSet;
  const IntPairToNodeSetMap& _X;
  const int _primaryIndex;
  const StringVector& _indexToAnatomicalSite;
  Digraph _G;
  Node _rootG;
  IntPairNodeMap _vertexLabelingG;
  IntPairToNodeMap _scToG;
  SubTreeVector _stateTrees;
  RefinementVector _refinements;
  
};

#endif // SOLUTIONGRAPH_H
