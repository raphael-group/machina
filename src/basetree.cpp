/*
 * basetree.cpp
 *
 *  Created on: 24-oct-2016
 *      Author: M. El-Kebir
 */

#include "basetree.h"
#include <lemon/bfs.h>

BaseTree::BaseTree()
  : _tree()
  , _nodeToId(_tree)
  , _idToNode()
  , _isLeaf(_tree, false)
  , _leafSet()
  , _leafSubset(_tree)
  , _level(_tree)
{
}

BaseTree::BaseTree(const BaseTree& other)
  : _tree()
  , _nodeToId(_tree)
  , _idToNode()
  , _isLeaf(_tree, false)
  , _leafSet()
  , _leafSubset(_tree)
  , _level(_tree)
{
  lemon::digraphCopy(other._tree, _tree)
    .node(other._root, _root)
    .nodeMap(other._nodeToId, _nodeToId)
    .run();
  
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    const std::string& str = _nodeToId[u];
    _idToNode[str] = u;
  }
  
  init();
}

BaseTree::BaseTree(const Digraph& T,
                   Node root,
                   const StringNodeMap& label)
  : _tree()
  , _nodeToId(_tree)
  , _idToNode()
  , _isLeaf(_tree, false)
  , _leafSet()
  , _leafSubset(_tree)
  , _level(_tree)
{
  lemon::digraphCopy(T, _tree)
    .node(root, _root)
    .nodeMap(label, _nodeToId)
    .run();
  
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    const std::string& str = _nodeToId[u];
    assert(_idToNode.count(str) == 0);
    _idToNode[str] = u;
  }
  
  init();
}

void BaseTree::write(std::ostream& out) const
{
  for (ArcIt a(_tree); a != lemon::INVALID; ++a)
  {
    Node u = _tree.source(a);
    Node v = _tree.target(a);
    
    out << _nodeToId[u] << " " << _nodeToId[v] << std::endl;
  }
}

bool BaseTree::read(std::istream& in)
{
  _idToNode.clear();
  
  while (in.good())
  {
    std::string line;
    getline(in, line);
    
    if (line.empty())
      break;
    
    StringVector s;
    boost::split(s, line, boost::is_any_of("\t "));
    
    if (s.size() != 2)
    {
      std::cerr << "Error: line '" << line << "' incorrect number of vertices" << std::endl;
      return false;
    }
    
    std::string label_u = s[0];
    std::string label_v = s[1];
    
    if (_idToNode.count(label_u) == 0)
    {
      Node u = _tree.addNode();
      _idToNode[label_u] = u;
      _nodeToId[u] = label_u;
    }
    if (_idToNode.count(label_v) == 0)
    {
      Node v = _tree.addNode();
      _idToNode[label_v] = v;
      _nodeToId[v] = label_v;
    }
    
    Node u = _idToNode[label_u];
    Node v = _idToNode[label_v];
    
    _tree.addArc(u, v);
  }
  
  _root = lemon::INVALID;
  for (NodeIt node(_tree); node != lemon::INVALID; ++node)
  {
    if (InArcIt(_tree, node) == lemon::INVALID)
    {
      if (_root != lemon::INVALID)
      {
        std::cerr << "Error: multiple root node '" << _nodeToId[node]
                  << "' and '" << _nodeToId[_root] << "'" << std::endl;
        return false;
      }
      _root = node;
    }
  }
  
  if (_idToNode.empty())
  {
    std::cerr << "Error: empty tree" << std::endl;
    return false;
  }
  
  init();
  
  if (!isValid())
  {
    std::cerr << "Error: tree is not a valid tree" << std::endl;
    return false;
  }
  
  return true;
}

NodeList BaseTree::pathFromRoot(const Digraph& T, Node u)
{
  NodeList result;
  while ((InArcIt(T, u) != lemon::INVALID))
  {
    result.push_front(u);
    u = T.source(InArcIt(T, u));
  }
  result.push_front(u);
  return result;
}

NodeList BaseTree::path(Node u, Node v) const
{
  NodeList result;
  
  if (!isAncestor(u, v))
  {
    return result;
  }
  
  while (u != v)
  {
    result.push_front(v);
    v = _tree.source(InArcIt(_tree, v));
  }
  result.push_front(u);
  
  return result;
}

void BaseTree::init()
{
  _leafSet.clear();
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    _leafSubset[u].clear();
    
    OutArcIt a(_tree, u);
    if (a == lemon::INVALID)
    {
      _leafSet.insert(u);
      _isLeaf[u] = true;
    }
  }
  
  lemon::bfs(_tree).distMap(_level).run(_root);
  
  initLeafSubset(_root);
}

void BaseTree::initLeafSubset(Node u)
{
  NodeSet merged;
  for (OutArcIt a(_tree, u); a != lemon::INVALID; ++a)
  {
    Node v = _tree.target(a);
    initLeafSubset(v);
    merged.insert(_leafSubset[v].begin(), _leafSubset[v].end());
  }
  
  if (merged.empty())
  {
    // node is a leaf
    merged.insert(u);
  }
  _leafSubset[u] = merged;
}

bool BaseTree::readVertexLabeling(std::istream& in,
                                  const BaseTree& T,
                                  StringNodeMap& label)
{
  const Digraph& tree = T.tree();
  for (NodeIt u(tree); u != lemon::INVALID; ++u)
  {
    label[u] = "";
  }
  
  while (in.good())
  {
    std::string line;
    getline(in, line);
    
    if (line.empty())
      continue;
    
    StringVector s;
    boost::split(s, line, boost::is_any_of("\t "));
    
    std::string label_u = s[0];
    std::string label_s = s[1];
    
    Node u = T.getNodeByLabel(label_u);
    if (u == lemon::INVALID)
    {
      std::cerr << "Error: clone-tree vertex with label '" << label_u << "' does not exist" << std::endl;
      return false;
    }
    
    label[u] = label_s;
  }
  
  for (NodeIt u(tree); u != lemon::INVALID; ++u)
  {
    if (label[u].empty())
    {
      std::cerr << "Error: vertex '" << T.label(u) << "' left unlabeled" << std::endl;
      return false;
    }
  }
  
  return true;
}

bool BaseTree::readColorMap(std::istream& in,
                            StringToIntMap& colorMap)
{
  colorMap.clear();
  
  while (in.good())
  {
    std::string line;
    getline(in, line);
    
    if (line.empty())
      continue;
    
    StringVector s;
    boost::split(s, line, boost::is_any_of("\t "));
    
    std::string anatomicalSite = s[0];
    int color = boost::lexical_cast<int>(s[1]);
    
    if (colorMap.count(anatomicalSite) != 0)
    {
      std::cerr << "Error: anatomical site '" << anatomicalSite << "' is already assigned a color" << std::endl;
      return false;
    }

    colorMap[anatomicalSite] = color;
  }
  
  return true;
}

Node BaseTree::getLCA(const Digraph& T, const NodeSet& nodes)
{
  if (nodes.size() == 1)
  {
    return *nodes.begin();
  }
  
  NodeListVector allPaths;
  NodeListItVector iteratorList;
  for (Node node : nodes)
  {
    allPaths.push_back(pathFromRoot(T, node));
    iteratorList.push_back(allPaths.back().begin());
  }
  
  // all first elements in iteratorList are set to the root and hence are the same
  Node lca;
  bool same = true;
  while (same)
  {
    lca = *iteratorList.front();
    
    for (int i = 0; i < iteratorList.size(); ++i)
    {
      if (++iteratorList[i] == allPaths[i].end())
        same = false;
    }
    
    if (same)
    {
      Node first = *iteratorList.front();
      for (const NodeListIt& it : iteratorList)
      {
        if (*it != first)
        {
          same = false;
        }
      }
    }
  }
  
  return lca;
}

bool BaseTree::isValid() const
{
  if (_root == lemon::INVALID)
    return false;
  
  // check connectivity
  lemon::Bfs<Digraph> bfs(_tree);
  bfs.run(_root);
  
  for (NodeIt v(_tree); v != lemon::INVALID; ++v)
  {
    if (!bfs.reached(v))
      return false;
  }
  
  return true;
}

bool BaseTree::isConnected(const NodeSet& nodes) const
{
  Node lca = getLCA(nodes);
  
  if (nodes.count(lca) == 0)
  {
    return false;
  }
  
  for (Node u : nodes)
  {
    NodeList P = path(lca, u);
    for (Node v : P)
    {
      if (nodes.count(v) == 0)
      {
        return false;
      }
    }
  }
  
  return true;
}

void BaseTree::writeDOT(std::ostream& out) const
{
  out << "digraph T {" << std::endl;
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (_isLeaf[u])
    {
      out << "\t\t" << _tree.id(u) << " [label=\"" << _nodeToId[u] << "\"]" << std::endl;
    }
  }
  out << "\t}" << std::endl;
  
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (!_isLeaf[u])
    {
      out << "\t" << _tree.id(u) << " [label=\"" << _nodeToId[u] << "\"]" << std::endl;
    }
  }
  
  for (ArcIt a(_tree); a != lemon::INVALID; ++a)
  {
    out << "\t" << _tree.id(_tree.source(a)) << " -> " << _tree.id(_tree.target(a)) << std::endl;
  }
  
  out << "}" << std::endl;
}

