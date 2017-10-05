/*
 *  statetree.h
 *
 *   Created on: 28-sep-2015
 *       Author: M. El-Kebir
 */

#ifndef STATETREE_H
#define STATETREE_H

#include "utils.h"
#include "migrationtree.h"

class StateTree
{
private:
  StateTree(); // hide default constructor
  
public:
  DIGRAPH_TYPEDEFS(Digraph);
  
  StateTree(int k);
  
  StateTree(const StateTree& other);
  
  StateTree& operator=(const StateTree& other);
  
  StateTree(const IntVector& pi);
  
  StateTree(const MigrationTree& T,
            const StringVector& indexToAnatomicalSite);
  
  int k() const
  {
    return _k;
  }
  
  int parent(int i) const
  {
    assert(0 <= i && i < _k);
    Node v_i = _stateToNode[i];
    Arc a = InArcIt(_S, v_i);
    if (a == lemon::INVALID)
    {
      if (v_i == _root)
      {
        return -1;
      }
      else
      {
        return -2;
      }
    }
    else
    {
      return _nodeToState[_S.source(a)];
    }
  }
  
  int numVertices() const
  {
    return lemon::countNodes(_S);
  }
  
  bool isPresent(int i) const
  {
    assert(0 <= i && i < _k);
    return parent(i) != -2;
  }
  
  bool isParent(int i, int j) const
  {
    assert(0 <= i && i < _k);
    assert(1 <= j && j < _k);

    Node v_i = _stateToNode[i];
    Node v_j = _stateToNode[j];
    
    Arc a = InArcIt(_S, v_j);
    if (a == lemon::INVALID)
    {
      return false;
    }
    else
    {
      return v_i == _S.source(a);
    }
  }
  
  bool areSibblings(int i, int j) const
  {
    assert(1 <= i && i < _k);
    assert(1 <= j && j < _k);
    
    Node v_i = _stateToNode[i];
    Node v_j = _stateToNode[j];
    
    Arc a1 = InArcIt(_S, v_i);
    if (a1 == lemon::INVALID)
      return false;
    
    Node v_pi_i = _S.source(a1);
    
    Arc a2 = InArcIt(_S, v_j);
    if (a2 == lemon::INVALID)
      return false;
    Node v_pi_j = _S.source(a2);
    
    return v_pi_i == v_pi_j;
  }
  
  bool isAncestor(int i, int j) const
  {
    assert(0 <= i && i < _k);
    assert(0 <= j && j < _k);
    
    Node v_i = _stateToNode[i];
    
    Node v = _stateToNode[j];
    while (v != v_i && v != _root)
    {
      v = _S.source(InArcIt(_S, v));
    }
    
    return v == v_i;
  }
  
  bool isDescendant(int i, int j) const
  {
    return isAncestor(j, i);
  }
  
  bool isIncomparable(int i, int j) const
  {
    return !isAncestor(i, j) && !isAncestor(j, i);
  }
  
  const IntSet& D(int i) const
  {
    assert(0 <= i && i < _k);
    return _D[_stateToNode[i]];
  }
  
  const std::string& label(int i) const
  {
    assert(0 <= i && i < _k);
    return _label[_stateToNode[i]];
  }
  
  void setLabel(int i, const std::string& l)
  {
    _label[_stateToNode[i]] = l;
  }
  
  const Digraph& S() const
  {
    return _S;
  }
  
  int state(Node v_i) const
  {
    assert(v_i != lemon::INVALID);
    return _nodeToState[v_i];
  }
  
  void writeEdgeList(std::ostream& out) const;
  
  void writeDOT(std::ostream& out) const;
  
  friend std::ostream& operator<<(std::ostream& out, const StateTree& S);
  friend std::istream& operator>>(std::istream& in, StateTree& S);
  
private:
  typedef std::vector<Node> NodeVector;
  typedef Digraph::NodeMap<IntSet> IntSetNodeMap;
  typedef Digraph::NodeMap<std::string> StringNodeMap;
  
private:
  int _k;
  Digraph _S;
  Node _root;
  StringNodeMap _label;
  
  NodeVector _stateToNode;
  IntNodeMap _nodeToState;
  
  IntSetNodeMap _D;
  
  void initD(Node v_i);
  void init(const IntVector& pi);
};

typedef std::vector<StateTree> StateTreeVector;
  
std::ostream& operator<<(std::ostream& out, const StateTree& S);
std::istream& operator>>(std::istream& in, StateTree& S);

#endif // STATETREE_H
