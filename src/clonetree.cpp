/*
 * clonetree.cpp
 *
 *  Created on: 19-oct-2016
 *      Author: M. El-Kebir
 */

#include "clonetree.h"
#include <iomanip>

CloneTree::CloneTree()
  : BaseTree()
  , _l(_tree)
{
}

CloneTree::CloneTree(const CloneTree& other)
  : BaseTree(other)
  , _l(_tree)
{
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    const std::string& lbl = label(u);
    Node other_u = other.getNodeByLabel(lbl);
    
    _l[u] = other._l[other_u];
  }
}

CloneTree& CloneTree::operator =(const CloneTree& other)
{
  if (this != &other)
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
    
    for (NodeIt u(_tree); u != lemon::INVALID; ++u)
    {
      const std::string& lbl = label(u);
      Node other_u = other.getNodeByLabel(lbl);
      
      _l[u] = other._l[other_u];
    }
  }
  
  return *this;
}

CloneTree::CloneTree(const Digraph& T,
                     Node root,
                     const StringNodeMap& label,
                     const StringNodeMap& l)
  : BaseTree(T, root, label)
  , _l(_tree)
{
  for (NodeIt v(T); v != lemon::INVALID; ++v)
  {
    if (OutArcIt(T, v) == lemon::INVALID)
    {
      const std::string& label_v = label[v];
      _l[getNodeByLabel(label_v)] = l[v];
    }
  }
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
//    if (!_isLeaf[u])
//    {
//      std::cerr << "Error: clone-tree vertex with label '" << label_u << "' is not a leaf" << std::endl;
//      return false;
//    }
    
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
                         const StringNodeMap& lPlus,
                         const StringToIntMap& colorMap,
                         const DoubleNodeMap& U,
                         const IntNodeMap& characterLabel) const
{
  out << "digraph T {" << std::endl;
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (_isLeaf[u])
    {
      out << "\t\t" << _tree.id(u) << " [width=0.2,height=0.2,penwidth=3,colorscheme=set19,color="
          << colorMap.find(lPlus[u])->second
          << ",label=\"";
      out << std::setprecision(3) << 100 * U[u] << "%";
      out << "\"]" << std::endl;
    }
  }
  out << "\t}" << std::endl;
  
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (!_isLeaf[u])
    {
      out << "\t" << _tree.id(u) << " [width=0.2,height=0.2,penwidth=3,colorscheme=set19,color="
      << colorMap.find(lPlus[u])->second << ",label=\"\"]" << std::endl;
    }
  }
  out << "\tinv [style=\"invis\"]" << std::endl;
  out << "\tinv -> "
  << _tree.id(_root) <<"[penwidth=3,colorscheme=set19,color="
  << colorMap.find(lPlus[_root])->second << ",label=\"" << characterLabel[_root]
  << "\",fontsize=18]" << std::endl;
  
  for (ArcIt a(_tree); a != lemon::INVALID; ++a)
  {
    Node u = _tree.source(a);
    Node v = _tree.target(a);
    
    const std::string& s_u = lPlus[u];
    const std::string& s_v = lPlus[v];
    
    out << "\t" << _tree.id(u) << " -> " << _tree.id(v);
    if (s_u == s_v)
    {
      out << " [fontsize=18,penwidth=3,colorscheme=set19,color=" << colorMap.find(s_u)->second;
    }
    else
    {
      out << " [fontsize=18,penwidth=3,colorscheme=set19,color=\"" << colorMap.find(s_u)->second << ";0.5:" << colorMap.find(s_v)->second << "\"";
    }
    if (characterLabel[_tree.target(a)] != -1)
    {
      out << ",label=\"" << characterLabel[_tree.target(a)] << "\"";
    }
    out << "]";
    out << std::endl;
  }
  
  out << "}" << std::endl;
}


void CloneTree::writeDOT(std::ostream& out,
                         const StringNodeMap& lPlus,
                         const StringToIntMap& colorMap,
                         const DoubleVectorNodeMap& U,
                         const StringNodeMap& characterLabel) const
{
  out << "digraph T {" << std::endl;
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (_isLeaf[u])
    {
      out << "\t\t" << _tree.id(u) << " [width=0.2,height=0.2,penwidth=3,colorscheme=set19,color="
          << colorMap.find(lPlus[u])->second
          << ",label=\"";
      
      const DoubleVector& vecU = U[u];
      bool first = true;
      for (int p = 0; p < vecU.size(); ++p)
      {
//        if (g_tol.nonZero(vecU[p]))
        {
          if (first)
            first = false;
          else
            out << "\\n";
          
          out << std::setprecision(3) << 100 * vecU[p] << "%";
        }
      }
      out << "\"]" << std::endl;
    }
  }
  out << "\t}" << std::endl;
  
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (!_isLeaf[u])
    {
      out << "\t" << _tree.id(u) << " [width=0.2,height=0.2,penwidth=3,colorscheme=set19,color="
          << colorMap.find(lPlus[u])->second << ",label=\"\"]" << std::endl;
    }
  }
  out << "\tinv [style=\"invis\"]" << std::endl;
  out << "\tinv -> "
      << _tree.id(_root) <<"[penwidth=3,colorscheme=set19,color="
      << colorMap.find(lPlus[_root])->second << ",label=\"" << characterLabel[_root]
      << "\",fontsize=18]" << std::endl;
  
  for (ArcIt a(_tree); a != lemon::INVALID; ++a)
  {
    Node u = _tree.source(a);
    Node v = _tree.target(a);
    
    const std::string& s_u = lPlus[u];
    const std::string& s_v = lPlus[v];
    
    out << "\t" << _tree.id(u) << " -> " << _tree.id(v);
    if (s_u == s_v)
    {
      out << " [fontsize=18,penwidth=3,colorscheme=set19,color=" << colorMap.find(s_u)->second;
    }
    else
    {
      out << " [fontsize=18,penwidth=3,colorscheme=set19,color=\"" << colorMap.find(s_u)->second << ";0.5:" << colorMap.find(s_v)->second << "\"";
    }
    if (!characterLabel[_tree.target(a)].empty())
    {
      out << ",label=\"" << characterLabel[_tree.target(a)] << "\"";
    }
    out << "]";
    out << std::endl;
  }
  
  out << "}" << std::endl;
}

void CloneTree::writeDOT(std::ostream& out,
                         const StringNodeMap& lPlus,
                         const StringToIntMap& colorMap,
                         const DoubleVectorNodeMap& F,
                         const DoubleNodeMap& U) const
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
          << ",label=\"" << _nodeToId[u];
      out << "\\n" << U[u] << "\\n" << lPlus[u] << "\"]" << std::endl;
    }
  }
  out << "\t}" << std::endl;
  
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (!_isLeaf[u])
    {
      out << "\t" << _tree.id(u) << " [penwidth=3,colorscheme=set19,color="
          << colorMap.find(lPlus[u])->second << ",label=\"" << _nodeToId[u];
      for (int s = 0; s < F[u].size(); ++s)
      {
        out << "\\n" << F[u][s];
      }
      out << "\"]" << std::endl;
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

void CloneTree::writeDOT(std::ostream& out) const
{
  StringToIntMap colorMap = generateColorMap();
  writeDOT(out, colorMap);
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

void CloneTree::mergeSameSiblingLeaves()
{
  typedef std::map<std::string, NodeList> AnatomicalSiteToNodeList;
  char buf[1024];
  
  for (NodeIt v(_tree); v != lemon::INVALID; ++v)
  {
    if (_isLeaf[v])
    {
      continue;
    }
    
    AnatomicalSiteToNodeList map;
    for (OutArcIt a(_tree, v); a != lemon::INVALID; ++a)
    {
      Node w = _tree.target(a);
      if (_isLeaf[w])
      {
        map[l(w)].push_back(w);
      }
    }
    
    int j = 0;
    for (const auto& sV : map)
    {
      const NodeList& V = sV.second;
      Node grandParent = _tree.source(InArcIt(_tree, V.front()));
      
      if (V.size() > 1)
      {
        // make a backbone
        NodeVector backbone;
        for (int i = 0; i < V.size() - 1; ++i)
        {
          backbone.push_back(_tree.addNode());
          
          snprintf(buf, 1024, "%s^%d^%d", _nodeToId[grandParent].c_str(), j, i);
          std::string lbl_parent =  buf;
          _nodeToId[backbone.back()] = lbl_parent;
          _idToNode[lbl_parent] = backbone.back();
        }
        
        _tree.addArc(grandParent, backbone.front());
        
        for (int i = 1; i < V.size() - 1; ++i)
        {
          _tree.addArc(backbone[i-1], backbone[i]);
        }
        
        int i = 0;
        for (Node w : V)
        {
          _tree.erase(InArcIt(_tree, w));
          _tree.addArc(backbone[i], w);
          if (i < V.size() - 2)
            ++i;
        }
        ++j;
      }
    }
  }
  
  init();
}

void CloneTree::writeDOT(std::ostream& out,
                         const StringNodeMap& lPlus,
                         const StringToIntMap& colorMap,
                         const DoubleNodeMap& U,
                         const StringNodeMap& characterLabel) const
{
  DoubleVectorNodeMap newU(_tree);
  for (NodeIt v(_tree); v != lemon::INVALID; ++v)
  {
    newU[v] = DoubleVector(1, U[v]);
  }
  
  writeDOT(out, lPlus, colorMap, newU, characterLabel);
}

void CloneTree::writeDOT(std::ostream& out,
                         const StringToIntMap& colorMap,
                         const DoubleNodeMap& U) const
{
  out << "digraph T {" << std::endl;
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (_isLeaf[u])
    {
      out << "\t\t" << _tree.id(u) << " [penwidth=3,colorscheme=set19,color="
          << colorMap.find(l(u))->second << ",label=\""
          << _nodeToId[u] << "\\n"
          << _l[u] << "\\n"
          << std::setprecision(2) << 100 * U[u]
          << "%"<< "\"]" << std::endl;
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

void CloneTree::writeDOT(std::ostream& out,
                         const StringNodeMap& lPlus,
                         const StringToIntMap& colorMap,
                         const DoubleVectorNodeMap& U) const
{
  out << "digraph T {" << std::endl;
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (_isLeaf[u])
    {
      out << "\t\t" << _tree.id(u) << " [penwidth=3,colorscheme=set19,color="
          << colorMap.find(l(u))->second << ",label=\""
          << _nodeToId[u] << "\\n"
          << _l[u] << "\\n";
      
      bool first = true;
      for (double usage : U[u])
      {
        if (first)
          first = false;
        else
          out << " ";
        out << std::setprecision(2) << 100 * usage << "%";
      }
      out << "\"]" << std::endl;
    }
  }
  out << "\t}" << std::endl;
  
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (!_isLeaf[u])
    {
      out << "\t" << _tree.id(u) << " [penwidth=3,colorscheme=set19,color="
          << colorMap.find(lPlus[u])->second << ",label=\""
          << _nodeToId[u] << "\"]" << std::endl;
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
                         const StringToIntMap& colorMap,
                         const DoubleVectorNodeMap& U) const
{
  out << "digraph T {" << std::endl;
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (_isLeaf[u])
    {
      out << "\t\t" << _tree.id(u) << " [penwidth=3,colorscheme=set19,color="
      << colorMap.find(l(u))->second << ",label=\""
      << _nodeToId[u] << "\\n"
      << _l[u] << "\\n";
      
      bool first = true;
      for (double usage : U[u])
      {
        if (first)
          first = false;
        else
          out << " ";
        out << std::setprecision(2) << 100 * usage << "%";
      }
      out << "\"]" << std::endl;
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

SplitSet CloneTree::getSplits() const
{
  SplitSet splits;
  const NodeSet allLeaves = leafSet();
  for (NodeIt v(_tree); v != lemon::INVALID; ++v)
  {
    const NodeSet& L_v = leafSubset(v);
    NodeSet L_complement_v;
    std::set_difference(allLeaves.begin(), allLeaves.end(),
                        L_v.begin(), L_v.end(),
                        std::inserter(L_complement_v,
                                      L_complement_v.begin()));
    
    StringSet split_v;
    StringSet complement_v;
    for (Node u : L_v)
    {
      const std::string& sStr = l(u);
      split_v.insert(sStr);
    }
    for (Node u : L_complement_v)
    {
      const std::string& sStr = l(u);
      complement_v.insert(sStr);
    }
    
    Split split;
    split.insert(split_v);
    split.insert(complement_v);
    splits.insert(split);
  }
  
  return splits;
}
