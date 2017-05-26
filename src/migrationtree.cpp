/*
 * migrationtree.cpp
 *
 *  Created on: 18-oct-2016
 *      Author: M. El-Kebir
 */

#include "migrationtree.h"
#include "gabowmyers.h"

MigrationTree::MigrationTree()
  : BaseTree()
{
}

MigrationTree::MigrationTree(const MigrationTree& other)
  : BaseTree()
{
  NodeNodeMap other2new(other._tree);
  
  lemon::digraphCopy(other._tree, _tree)
    .node(other._root, _root)
    .nodeMap(other._nodeToId, _nodeToId)
    .nodeRef(other2new)
    .run();
  
  init();
   
  for (const auto& kv : other._idToNode)
  {
    _idToNode[kv.first] = other2new[kv.second];
  }
}

MigrationTree::MigrationTree(const SubDigraph& T,
                             Node root,
                             const StringNodeMap& label)
  : BaseTree()
{
//  NodeNodeMap other2new(T);
 
  lemon::digraphCopy(T, _tree)
    .node(root, _root)
    .nodeMap(label, _nodeToId)
//    .nodeRef(other2new)
    .run();
  
  init();

  for (NodeIt x(_tree); x != lemon::INVALID; ++x)
  {
    _idToNode[_nodeToId[x]] = x;
  }
}

MigrationTree::~MigrationTree()
{
}

bool MigrationTree::readVertexLabeling(std::istream& in)
{
  return true;
}

void MigrationTree::writeVertexLabeling(std::ostream& out) const
{
}

void MigrationTree::writeDOT(std::ostream& out) const
{
  // construct sample to color scheme
  std::set<std::string> samples;
  for (NodeIt x(_tree); x != lemon::INVALID; ++x)
  {
    const std::string& s = l(x);
    samples.insert(s);
  }
  std::map<std::string, int> colorMap;
  int idx = 0;
  for (const std::string& s : samples)
  {
    colorMap[s] = ++idx;
  }
  
  out << "digraph barS {" << std::endl;
  
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt x(_tree); x != lemon::INVALID; ++x)
  {
    if (_isLeaf[x])
    {
      out << "\t\t" << _tree.id(x)
      << " [penwidth=3,colorscheme=set19,color="
      << colorMap[l(x)]
      << ",label=\"" << l(x) << "\"]" << std::endl;
    }
  }
  out << "\t}" << std::endl;
  
  for (NodeIt x(_tree); x != lemon::INVALID; ++x)
  {
    if (!_isLeaf[x])
    {
      out << "\t" << _tree.id(x)
      << " [penwidth=3,colorscheme=set19,color="
      << colorMap[l(x)]
      << ",label=\"" << l(x) << "\"]" << std::endl;
    }
  }
  
  for (ArcIt a(_tree); a != lemon::INVALID; ++a)
  {
    Node x = _tree.source(a);
    Node y = _tree.target(a);
    
    const std::string& s_x = l(x);
    const std::string& s_y = l(y);
    
    out << "\t" << _tree.id(x) << " -> " << _tree.id(y);
    if (s_x == s_y)
    {
      out << " [penwidth=3,colorscheme=set19,color=" << colorMap[s_x] << "]";
    }
    else
    {
      out << " [penwidth=3,colorscheme=set19,color=\"" << colorMap[s_x] << ";0.5:" << colorMap[s_y] << "\"]";
    }
    
    out << std::endl;
  }
  
  out << "}" << std::endl;
}

void MigrationTree::enumerate(const std::string& primary,
                              const StringSet& metastases,
                              MigrationTreeList& listBarS)
{
  // construct G
  Digraph G;
  StringNodeMap label(G);
  Node root = G.addNode();
  label[root] = primary;
  
  NodeVector met_nodes;
  for (const std::string& M : metastases)
  {
    Node met_node = G.addNode();
    label[met_node] = M;
    
    met_nodes.push_back(met_node);
    G.addArc(root, met_node);
  }
  
  int i = 0;
  for (const std::string& M1 : metastases)
  {
    int j = 0;
    for (const std::string& M2 : metastases)
    {
      if (M1 != M2)
      {
        G.addArc(met_nodes[i], met_nodes[j]);
      }
      ++j;
    }
    ++i;
  }
  
  GabowMyers gm(G, root);
//  gm.writeDOT(label, std::cout);
  gm.run();
  
  BoolNodeMap filterNodesT(G, false);
  BoolArcMap filterArcsT(G, false);
  SubDigraph T(G, filterNodesT, filterArcsT);
  
  const int n = gm.getNrTrees();
  
  for (int i = 0; i < n; ++i)
  {
    gm.result(T, i);
    MigrationTree barS(T, root, label);
    listBarS.push_back(barS);
  }
}

std::ostream& operator<<(std::ostream& out, const MigrationTreeList& listBarS)
{
  if (listBarS.empty())
  {
    return out;
  }
  
  out << listBarS.front().getNrSamples() << " #samples" << std::endl;
  out << listBarS.size() << " #condensed sample trees" << std::endl;
  
  for (const MigrationTree& barS : listBarS)
  {
    barS.write(out);
    out << std::endl;
  }
  
  return out;
}

std::istream& operator>>(std::istream& in, MigrationTreeList& listBarS)
{
  std::string line;
  getline(in, line);
  
  int nrSamples = -1;
  std::stringstream ss(line);
  ss >> nrSamples;
  
  if (nrSamples < 1)
  {
    throw std::runtime_error("Error: number of samples should be greater than 0");
  }
  
  int nrMigrationTrees = -1;
  getline(in, line);
  ss.str(line);
  ss >> nrMigrationTrees;
  
  if (nrMigrationTrees < 1)
  {
    throw std::runtime_error("Error: number of condensed sample trees should be greater than 0");
  }
  
  for (int i = 0; i < nrMigrationTrees; ++i)
  {
    MigrationTree barS;
    if (!barS.read(in))
    {
      throw std::runtime_error("Error: reading condensed sample tree failed");
    }
//    getline(in, line);
    listBarS.push_back(barS);
  }
  
  return in;
}
