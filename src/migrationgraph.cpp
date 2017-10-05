/*
 * migrationgraph.cpp
 *
 *  Created on: 24-oct-2016
 *      Author: M. El-Kebir
 */

#include "migrationgraph.h"
#include <lemon/bfs.h>

const int MigrationGraph::_nrPatterns = 4;

MigrationGraph::MigrationGraph()
  : _G()
  , _root(lemon::INVALID)
  , _nodeToId(_G)
  , _idToNode()
{
}

MigrationGraph::MigrationGraph(const Digraph& G,
                               Node root,
                               const StringNodeMap& id)
: _G()
, _root(lemon::INVALID)
, _nodeToId(_G)
, _idToNode()
{
  lemon::digraphCopy(G, _G)
    .node(root, _root)
    .nodeMap(id, _nodeToId)
    .run();
  
  for (NodeIt u(_G); u != lemon::INVALID; ++u)
  {
    const std::string& str = _nodeToId[u];
    assert(_idToNode.count(str) == 0);
    _idToNode[str] = u;
  }
}

MigrationGraph::MigrationGraph(const MigrationGraph& other)
  : _G()
  , _root(lemon::INVALID)
  , _nodeToId(_G)
  , _idToNode()
{
  lemon::digraphCopy(other._G, _G)
    .node(other._root, _root)
    .nodeMap(other._nodeToId, _nodeToId)
    .run();
  
  for (NodeIt u(_G); u != lemon::INVALID; ++u)
  {
    const std::string& str = _nodeToId[u];
    _idToNode[str] = u;
  }
}

MigrationGraph::MigrationGraph(const CloneTree& T,
                               const StringNodeMap& lPlus)
  : _G()
  , _root(lemon::INVALID)
  , _nodeToId(_G)
  , _idToNode()
{
  StringSet Sigma = T.getAnatomicalSites();
  
  for (const std::string& s : Sigma)
  {
    Node x = _G.addNode();
    _nodeToId[x] = s;
    _idToNode[s] = x;
  }
  
  assert(_idToNode.count(lPlus[T.root()]) == 1);
  _root = _idToNode[lPlus[T.root()]];
  
  for (ArcIt a(T.tree()); a != lemon::INVALID; ++a)
  {
    Node u = T.tree().source(a);
    Node v = T.tree().target(a);
    
    assert(_idToNode.count(lPlus[u]) == 1);
    assert(_idToNode.count(lPlus[v]) == 1);

    Node x = _idToNode[lPlus[u]];
    Node y = _idToNode[lPlus[v]];
    
    if (x != y)
    {
      _G.addArc(x, y);
    }
  }
}

bool MigrationGraph::read(std::istream& in)
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
      Node u = _G.addNode();
      _idToNode[label_u] = u;
      _nodeToId[u] = label_u;
    }
    if (_idToNode.count(label_v) == 0)
    {
      Node v = _G.addNode();
      _idToNode[label_v] = v;
      _nodeToId[v] = label_v;
    }
    
    Node u = _idToNode[label_u];
    Node v = _idToNode[label_v];
    
    _G.addArc(u, v);
  }
  
  _root = lemon::INVALID;
  for (NodeIt node(_G); node != lemon::INVALID; ++node)
  {
    if (InArcIt(_G, node) == lemon::INVALID)
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

  if (!isConnected())
  {
    std::cerr << "Error: migration graph is not connected" << std::endl;
    return false;
  }

  return true;
}

bool MigrationGraph::isConnected() const
{
  if (_root == lemon::INVALID)
    return false;
  
  // check connectivity
  lemon::Bfs<Digraph> bfs(_G);
  bfs.run(_root);
  
  for (NodeIt v(_G); v != lemon::INVALID; ++v)
  {
    if (!bfs.reached(v))
      return false;
  }
  
  return true;
}

void MigrationGraph::write(std::ostream& out) const
{
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    Node u = _G.source(a);
    Node v = _G.target(a);
    
    out << _nodeToId[u] << " " << _nodeToId[v] << std::endl;
  }
}

void MigrationGraph::writeDOT(std::ostream& out) const
{
  StringSet Sigma = getAnatomicalSites();
  
  StringToIntMap colorMap;
  int idx = 0;
  for (const std::string& s : Sigma)
  {
    colorMap[s] = ++idx;
  }

  writeDOT(out, colorMap);
}

void MigrationGraph::writeDOT(std::ostream& out, const StringToIntMap& colorMap) const
{
  out << "digraph barS {" << std::endl;
  
  /// first root
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  out << "\t" << _G.id(_root)
      << " [shape=box,penwidth=3,colorscheme=set19,color="
      << colorMap.find(l(_root))->second
      << ",label=\"" << l(_root) << "\"]" << std::endl;
  out << "\t}" << std::endl;
  
  // then the leaves
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt x(_G); x != lemon::INVALID; ++x)
  {
    if (lemon::countOutArcs(_G, x) == 0)
    {
      out << "\t\t" << _G.id(x)
          << " [shape=box,penwidth=3,colorscheme=set19,color="
          << colorMap.find(l(x))->second
          << ",label=\"" << l(x) << "\"]" << std::endl;
    }
  }
  out << "\t}" << std::endl;
  
  // finally the rest of the nodes
  for (NodeIt x(_G); x != lemon::INVALID; ++x)
  {
    if (lemon::countOutArcs(_G, x) != 0 && x != _root)
    {
      out << "\t" << _G.id(x)
      << " [shape=box,penwidth=3,colorscheme=set19,color="
      << colorMap.find(l(x))->second
      << ",label=\"" << l(x) << "\"]" << std::endl;
    }
  }
  
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    Node x = _G.source(a);
    Node y = _G.target(a);
    
    const std::string& s_x = l(x);
    const std::string& s_y = l(y);
    
    out << "\t" << _G.id(x) << " -> " << _G.id(y);
    if (s_x == s_y)
    {
      out << " [penwidth=3,colorscheme=set19,color=" << colorMap.find(s_x)->second << "]";
    }
    else
    {
      out << " [penwidth=3,colorscheme=set19,color=\"" << colorMap.find(s_x)->second << ";0.5:" << colorMap.find(s_y)->second << "\"]";
    }
    
    out << std::endl;
  }
  
  out << "}" << std::endl;
}
