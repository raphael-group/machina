/*
 *  statetree.cpp
 *
 *   Created on: 28-sep-2015
 *       Author: M. El-Kebir
 */

#include "statetree.h"

#include <lemon/connectivity.h>
#include <stdio.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

StateTree::StateTree(const MigrationTree& T,
                     const StringVector& indexToAnatomicalSite)
  : _k(T.getNrSamples())
  , _S()
  , _root(lemon::INVALID)
  , _label(_S)
  , _stateToNode(_k, lemon::INVALID)
  , _nodeToState(_S)
  , _D(_S)
{
  assert(indexToAnatomicalSite[0] == T.label(T.root()));
  NodeNodeMap T2S(T.tree(), lemon::INVALID);
  
  StringToIntMap invNodeMap;
  _root = _S.addNode();
  _stateToNode[0] = _root;
  _nodeToState[_root] = 0;
  _label[_root] = T.label(T.root());
  T2S[T.root()] = _root;
  invNodeMap[_label[_root]] = 0;

  for (int i = 0; i < indexToAnatomicalSite.size(); ++i)
  {
    const std::string& sStr = indexToAnatomicalSite[i];
    if (sStr == _label[_root])
    {
      continue;
    }
    
    Node vv = _S.addNode();
    _label[vv] = sStr;
    _nodeToState[vv] = i;
    _stateToNode[i] = vv;
    invNodeMap[_label[vv]] = i;
  }
  
  for (ArcIt uv(T.tree()); uv != lemon::INVALID; ++uv)
  {
    Node u = T.tree().source(uv);
    Node v = T.tree().target(uv);
    
    assert(invNodeMap.count(T.label(u)) == 1);
    assert(invNodeMap.count(T.label(v)) == 1);
    
    int i = invNodeMap[T.label(u)];
    int j = invNodeMap[T.label(v)];
    
    Node uu = _stateToNode[i];
    Node vv = _stateToNode[j];
    
    _S.addArc(uu, vv);
  }
  
  initD(_root);
  
  assert(lemon::dag(_S));
}

StateTree::StateTree(int k)
  : _k(k)
  , _S()
  , _root(lemon::INVALID)
  , _label(_S)
  , _stateToNode(_k, lemon::INVALID)
  , _nodeToState(_S)
  , _D(_S)
{
  IntVector pi(_k, -1);
  for (int i = 1; i < _k; ++i)
  {
    pi[i] = i - 1;
  }
  
  init(pi);
}

void StateTree::init(const IntVector& pi)
{
  // state 0 should be the root
  assert(pi[0] == -1);
  char buf[1024];
  
  _root = _S.addNode();
  _stateToNode[0] = _root;
  _nodeToState[_root] = 0;
  snprintf(buf, 1024, "%d", 0);
  _label[_root] = buf;
  
  // init nodes of S_c
  for (int i = 1; i < _k; ++i)
  {
//    if (0 <= pi[i] && pi[i] < _k)
//    if (pi[i] != -2)
    {
      Node v_i = _S.addNode();
      _stateToNode[i] = v_i;
      _nodeToState[v_i] = i;
      snprintf(buf, 1024, "%d", i);
      _label[v_i] = buf;
    }
  }
  
  // init edges of S_c
  for (int i = 1; i < _k; ++i)
  {
    int pi_i = pi[i];
    
    //assert(0 <= pi_i < _k);
    if (0 <= pi_i && pi_i < _k)
    {
      _S.addArc(_stateToNode[pi_i], _stateToNode[i]);
    }
  }
  
  initD(_root);
  
  assert(lemon::dag(_S));
}
  
StateTree::StateTree(const IntVector& pi)
  : _k(pi.size())
  , _S()
  , _root(lemon::INVALID)
  , _label(_S, "")
  , _stateToNode(_k, lemon::INVALID)
  , _nodeToState(_S)
  , _D(_S)
{
  init(pi);
}
  
StateTree::StateTree(const StateTree& other)
  : _k(other._k)
  , _S()
  , _root(lemon::INVALID)
  , _label(_S)
  , _stateToNode(_k, lemon::INVALID)
  , _nodeToState(_S)
  , _D(_S)
{
  lemon::digraphCopy(other._S, _S)
    .node(other._root, _root)
    .nodeMap(other._nodeToState, _nodeToState)
    .nodeMap(other._D, _D)
    .nodeMap(other._label, _label)
    .run();
  
  for (NodeIt v_i(_S); v_i != lemon::INVALID; ++v_i)
  {
    _stateToNode[_nodeToState[v_i]] = v_i;
  }
}
  
StateTree& StateTree::operator=(const StateTree& other)
{
  if (this != &other)
  {
    _k = other._k;
    _stateToNode = NodeVector(_k, lemon::INVALID);
    
    lemon::digraphCopy(other._S, _S)
      .node(other._root, _root)
      .nodeMap(other._nodeToState, _nodeToState)
      .nodeMap(other._label, _label)
      .nodeMap(other._D, _D)
      .run();
    
    for (NodeIt v_i(_S); v_i != lemon::INVALID; ++v_i)
    {
      _stateToNode[_nodeToState[v_i]] = v_i;
    }
  }
  return *this;
}
  
void StateTree::writeEdgeList(std::ostream& out) const
{
  bool first = true;
  for (ArcIt a_ij(_S); a_ij != lemon::INVALID; ++a_ij)
  {
    if (first)
    {
      first = false;
    }
    else
    {
      out << " ; ";
    }
    
    Node v_i = _S.source(a_ij);
    Node v_j = _S.target(a_ij);
    
    out << _label[v_i] << " -> " << _label[v_j];
  }
}
  
void StateTree::writeDOT(std::ostream& out) const
{
  out << "digraph S {" << std::endl;
  
  for (NodeIt v_i(_S); v_i != lemon::INVALID; ++v_i)
  {
    int i = _nodeToState[v_i];
    
    out << "\t" << i << " [label=\"" << _label[v_i] << "\"]" << std::endl;
  }
  
  for (ArcIt a_ij(_S); a_ij != lemon::INVALID; ++a_ij)
  {
    int i = _nodeToState[_S.source(a_ij)];
    int j = _nodeToState[_S.target(a_ij)];
    
    out << "\t" << i << " -> " << j << std::endl;
  }
  
  out << "}" << std::endl;
}
  
void StateTree::initD(Node v_i)
{
  IntSet& D_i = _D[v_i];
  D_i.clear();
  
  D_i.insert(_nodeToState[v_i]);
  
  for (OutArcIt a(_S, v_i); a != lemon::INVALID; ++a)
  {
    Node v_j = _S.target(a);
    
    initD(v_j);
    D_i.insert(_D[v_j].begin(), _D[v_j].end());
  }
}
  
std::ostream& operator<<(std::ostream& out, const StateTree& S)
{
  // output pi
  out << -1;
  for (int i = 1; i < S.k(); ++i)
  {
    out << " " << S.parent(i);
  }
  out << std::endl;
  
  // output labels
  out << S.label(0);
  for (int i = 1; i < S.k(); ++i)
  {
    out << " " << S.label(i);
  }
  out << std::endl;

  return out;
}
  
std::istream& operator>>(std::istream& in, StateTree& S)
{
  typedef std::vector<std::string> StringVector;
  
  std::string line;
  getline(in, line);
  std::stringstream ss(line);
  
  StringVector s;
  boost::split(s, line, boost::is_any_of(" \t"));
  
  IntVector pi(s.size(), -1);
  for (int i = 0; i < s.size(); ++i)
  {
    ss >> pi[i];
  }
  
  S = StateTree(pi);
  
  getline(in, line);
  ss.clear();
  ss.str(line);
  for (int i = 0; i < S.k(); ++i)
  {
    ss >> S._label[S._stateToNode[i]];
  }
  
  return in;
}
