/*
 * charactertree.cpp
 *
 *  Created on: 13-oct-2016
 *      Author: M. El-Kebir
 */

#include "charactertree.h"

CharacterTree::CharacterTree()
  : BaseTree()
  , _stateVector(_tree)
  , _nrStates()
  , _indexToCharacter()
  , _characterToIndex()
{
}

CharacterTree::CharacterTree(const CloneTree& T,
                             const std::string& primary,
                             StringVector& stateToAnatomicalSite,
                             StringToIntMap& anatomicalSiteToState)
  : BaseTree(T)
  , _stateVector(_tree, IntVector(1, -1))
  , _nrStates(1, T.getNrAnatomicalSites())
  , _indexToCharacter()
  , _characterToIndex()
{
  _characterToIndex["SITE"] = 0;
  _indexToCharacter.push_back("SITE");
  
  // set the states, starting with the primary
  assert(T.getAnatomicalSites().count(primary) == 1);
  stateToAnatomicalSite.push_back(primary);
  anatomicalSiteToState[primary] = 0;
  
  for (const std::string& s : T.getAnatomicalSites())
  {
    if (s == primary)
      continue;
    
    anatomicalSiteToState[s] = stateToAnatomicalSite.size();
    stateToAnatomicalSite.push_back(s);
  }
  
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (_isLeaf[u])
    {
      Node uu = T.getNodeByLabel(label(u));
      assert(uu != lemon::INVALID);
      
      const std::string& s = T.l(uu);
      assert(anatomicalSiteToState.count(s) == 1);
      _stateVector[u] = IntVector(1, anatomicalSiteToState[s]);
    }
    else
    {
      _stateVector[u] = IntVector(1, -1);
    }
  }
}

CharacterTree::CharacterTree(const CharacterTree& other)
  : BaseTree(other)
  , _stateVector(_tree)
  , _nrStates(other._nrStates)
  , _indexToCharacter(other._indexToCharacter)
  , _characterToIndex(other._characterToIndex)
{
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    const std::string& lbl = label(u);
    Node other_u = other.getNodeByLabel(lbl);
    
    _stateVector[u] = other._stateVector[other_u];
  }
}

void CharacterTree::eraseLeaf(Node v)
{
  assert(isLeaf(v));
  
  Node u = parent(v);
  
  // 1. remove v
  const std::string& lbl_v = label(v);
  _idToNode.erase(lbl_v);
  _tree.erase(v);
  
  if (lemon::countOutArcs(_tree, u) == 1)
  {
    Node w = _tree.target(OutArcIt(_tree, u));
   
    // 2. collapse u and w
    const std::string lbl_u = label(u);
    const std::string& lbl_w = label(w);
    
    _idToNode[lbl_w] = u;
    _nodeToId[u] = lbl_w;
    _idToNode.erase(lbl_u);
    
    for (OutArcIt a(_tree, w); a != lemon::INVALID; ++a)
    {
      Node w_child = _tree.target(a);
      _tree.addArc(u, w_child);
    }
    _stateVector[u] = _stateVector[w];
    _tree.erase(w);
  }
  
  init();
}

bool CharacterTree::readLeafLabeling(std::istream& in)
{
  typedef Digraph::NodeMap<StringSet> StringSetNodeMap;
  
  StringSetNodeMap characterSet(_tree);
  StringSet characters;

  // first line is the anatomical site set
  std::string line;
  getline(in, line);

  boost::split(_indexToCharacter, line, boost::is_any_of("\t "));
  const int n = _indexToCharacter.size();
  
  _characterToIndex.clear();
  for (int i = 0; i < n; ++i)
  {
    const std::string& c = _indexToCharacter[i];
    _characterToIndex[c] = i;
  }
  
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    _stateVector[u] = IntVector(n, 0);
  }

  while (in.good())
  {
    getline(in, line);
    if (line.empty())
      continue;
    
    StringVector s;
    boost::split(s, line, boost::is_any_of("\t "));
    
    std::string label_u = s[0];
    
    if (_idToNode.count(label_u) == 0)
    {
      std::cerr << "Error: vertex with label '" << label_u << "' does not exist" << std::endl;
      return false;
    }
    
    Node u = _idToNode[label_u];
    if (!_isLeaf[u])
    {
      std::cerr << "Error: vertex with label '" << label_u << "' is not a leaf" << std::endl;
      return false;
    }
    
    for (int i = 0; i < n; ++i)
    {
      int state = boost::lexical_cast<int>(s[i+1]);
      _stateVector[u][i] = state;
    }
  }
  
  _nrStates = IntVector(n, 0);
  for (int i = 0; i < n; ++i)
  {
    for (NodeIt u(_tree); u != lemon::INVALID; ++u)
    {
      _nrStates[i] = std::max(_nrStates[i], _stateVector[u][i]);
    }
  }
  for (int i = 0; i < n; ++i)
  {
    ++_nrStates[i];
  }
  
  return true;
}
