/*
 * clonetree.cpp
 *
 *  Created on: 19-oct-2016
 *      Author: M. El-Kebir
 */

#include "clonetree.h"

CloneTree::CloneTree()
  : BinaryTree()
  , _l(_tree)
{
}

CloneTree::CloneTree(const Digraph& T,
                     Node root,
                     const StringNodeMap& id,
                     const StringNodeMap& l)
  : BinaryTree(T, root, id)
  , _l(_tree)
{
  for (NodeIt v(T); v != lemon::INVALID; ++v)
  {
    if (OutArcIt(T, v) == lemon::INVALID)
    {
      const std::string& id_v = id[v];
      _l[getNodeByLabel(id_v)] = l[v];
    }
  }
}

CloneTree::CloneTree(const CloneTree& other)
  : BinaryTree(other)
  , _l(_tree)
{
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    const std::string& lbl = label(u);
    Node other_u = other.getNodeByLabel(lbl);
    
    _l[u] = other._l[other_u];
  }
}

CloneTree::CloneTree(const BoolVector& dyckWord, int nrSamples)
  : BinaryTree(dyckWord)
  , _l(_tree)
{
  char buf[1024];
  StringVector samples;
  samples.push_back("P");
  
  const int n = _leafSet.size();
  assert(n >= nrSamples);
  
  for (int i = 1; i < nrSamples; ++i)
  {
    snprintf(buf, 1024, "M_%d", i);
    samples.push_back(buf);
  }
  
  std::uniform_int_distribution<> dis(0, nrSamples - 1);
  for (int i = nrSamples; i < n; ++i)
  {
    int j = dis(g_rng);
    samples.push_back(samples[j]);
  }

  std::shuffle(samples.begin(), samples.end(), g_rng);
  for (Node x : _leafSet)
  {
    assert(!samples.empty());
    _l[x] = samples.back();
    samples.pop_back();
  }
  
  assert(samples.empty());
}

bool CloneTree::readLeafLabeling(std::istream& in)
{
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    _l[u] = "";
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
    
    if (_idToNode.count(label_u) == 0)
    {
      std::cerr << "Error: clone-tree vertex with label '" << label_u << "' does not exist" << std::endl;
      return false;
    }
    
    Node u = _idToNode[label_u];
    if (!_isLeaf[u])
    {
      std::cerr << "Error: clone-tree vertex with label '" << label_u << "' is not a leaf" << std::endl;
      return false;
    }
    
    _l[u] = label_s;
  }
  
  for (Node u : _leafSet)
  {
    if (_l[u].empty())
    {
      std::cerr << "Error: leaf '" << label(u) << "' left unlabeled" << std::endl;
      return false;
    }
  }
  
  return true;
}

void CloneTree::writeLeafLabeling(std::ostream& out) const
{
  for (Node u : _leafSet)
  {
    out << _nodeToId[u] << " " << _l[u] << std::endl;
  }
}

void CloneTree::writeVertexLabeling(std::ostream& out,
                                    const StringNodeMap& lPlus) const
{
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    out << _nodeToId[u] << " " << lPlus[u] << std::endl;
  }
}

void CloneTree::writeDOT(std::ostream& out) const
{
  StringToIntMap colorMap = generateColorMap();
  writeDOT(out, colorMap);
}

void CloneTree::writeDOT(std::ostream& out,
                         const StringNodeMap& lPlus) const
{
  StringToIntMap colorMap = generateColorMap();
  writeDOT(out, lPlus, colorMap);
}

void CloneTree::writeDOT(std::ostream& out,
                         const StringNodeMap& lPlus,
                         const StringToIntMap& colorMap) const
{
  out << "digraph T {" << std::endl;
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (_isLeaf[u])
    {
      out << "\t\t" << _tree.id(u) << " [penwidth=3,colorscheme=set19,color="
          << colorMap.find(lPlus[u])->second
          << ",label=\"" << _nodeToId[u] << "\\n" << lPlus[u] << "\"]" << std::endl;
    }
  }
  out << "\t}" << std::endl;
  
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (!_isLeaf[u])
    {
      out << "\t" << _tree.id(u) << " [penwidth=3,colorscheme=set19,color="
          << colorMap.find(lPlus[u])->second << ",label=\"" << _nodeToId[u] << "\"]" << std::endl;
    }
  }
  
  for (ArcIt a(_tree); a != lemon::INVALID; ++a)
  {
    Node u = _tree.source(a);
    Node v = _tree.target(a);
    
    const std::string& s_u = lPlus[u];
    const std::string& s_v = lPlus[v];
    
    out << "\t" << _tree.id(u) << " -> " << _tree.id(v);
    if (s_u == s_v)
    {
      out << " [penwidth=3,colorscheme=set19,color=" << colorMap.find(s_u)->second << "]";
    }
    else
    {
      out << " [penwidth=3,colorscheme=set19,color=\"" << colorMap.find(s_u)->second << ";0.5:" << colorMap.find(s_v)->second << "\"]";
    }
    
    out << std::endl;
  }
  
  out << "}" << std::endl;
}

void CloneTree::writeDOT(std::ostream& out,
                         const StringToIntMap& colorMap) const
{
  out << "digraph T {" << std::endl;
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (_isLeaf[u])
    {
      out << "\t\t" << _tree.id(u) << " [penwidth=3,colorscheme=set19,color="
          << colorMap.find(l(u))->second << ",label=\"" << _nodeToId[u] << "\\n" << _l[u] << "\"]" << std::endl;
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

