/*
 *  perfectphylotree.cpp
 *
 *   Created on: 28-sep-2015
 *       Author: M. El-Kebir
 */

#include "perfectphylotree.h"
#include "perfectphylograph.h"
#include <lemon/connectivity.h>
#include <boost/range/adaptor/reversed.hpp>

namespace gm {

PerfectPhyloTree::PerfectPhyloTree(const PerfectPhyloMatrix& A,
                                   const StateTreeVector& S)
  : _A(A)
  , _S(S)
  , _T()
  , _root(lemon::INVALID)
  , _charStateToNode(_A.n(), NodeVector(_A.k(), lemon::INVALID))
  , _nodeToCharState(_T, std::make_pair(-1, -1))
  , _label(_T)
{
  const int n = _A.n();
  const int k = _A.k();
  
  // add root node
  _root = _T.addNode();
  _nodeToCharState[_root] = std::make_pair(0, 0);
  for (int c = 0; c < n; ++c)
  {
    _charStateToNode[c][0] = _root;
  }
  
  // add the other n(k-1) nodes
  for (int c = 0; c < n; ++c)
  {
    for (int i = 1; i < k; ++i)
    {
      if (_A.defined(c, i))
      {
        Node v_ci = _T.addNode();
        _charStateToNode[c][i] = v_ci;
        _nodeToCharState[v_ci] = std::make_pair(c, i);
      }
    }
  }
  
  // add edges
  PerfectPhyloGraph G(_A);
//  G.writeDOT(std::cout);
  PerfectPhyloGraph::BoolNodeMap visited(G.G(), false);
  initArcs(G, G.root(), visited);
}
  
PerfectPhyloTree::PerfectPhyloTree(const PerfectPhyloTree& other)
  : _A(other._A)
  , _S(other._S)
  , _T()
  , _root(lemon::INVALID)
  , _charStateToNode(_A.n(), NodeVector(_A.k(), lemon::INVALID))
  , _nodeToCharState(_T, std::make_pair(-1, -1))
  , _label(_T)
{
  lemon::digraphCopy(other._T, _T)
    .nodeMap(other._nodeToCharState, _nodeToCharState)
    .nodeMap(other._label, _label)
    .node(other._root, _root)
    .run();
  
  for (NodeIt v_ci(_T); v_ci != lemon::INVALID; ++v_ci)
  {
    const IntPair& ci = _nodeToCharState[v_ci];
    _charStateToNode[ci.first][ci.second] = v_ci;
  }
}
  
void PerfectPhyloTree::setLabels(const RealTensor& F)
{
  char buf[1024];
  for (NodeIt v_ci(_T); v_ci != lemon::INVALID; ++v_ci)
  {
    const IntPair& ci = _nodeToCharState[v_ci];
    if (ci.first == 0 && ci.second == 0)
    {
      snprintf(buf, 1024, "(*,%s)",
               _S[ci.first].label(ci.second).c_str());
    }
    else
    {
      snprintf(buf, 1024, "(%s,%s)",
               F.getColLabel(ci.first).c_str(),
               _S[ci.first].label(ci.second).c_str());
    }
    _label[v_ci] = buf;
  }
}
  
int PerfectPhyloTree::edgeRecall(const PerfectPhyloTree& otherTree) const
{
  int recall = 0;
  for (NodeIt v_dj(_T); v_dj != lemon::INVALID; ++v_dj)
  {
    const std::string& label_v_dj = _label[v_dj];
    Node v_parent_dj = parent(v_dj);
    if (v_parent_dj != lemon::INVALID)
    {
      const std::string& label_v_parent_dj = _label[v_parent_dj];
      for (NodeIt u_dj(otherTree.T()); u_dj != lemon::INVALID; ++u_dj)
      {
        const std::string label_u_dj = otherTree.label(u_dj);
        if (label_v_dj == label_u_dj)
        {
          Node u_parent_dj = otherTree.parent(u_dj);
          
          if (u_parent_dj != lemon::INVALID &&
              otherTree.label(u_parent_dj) == label_v_parent_dj)
          {
            ++recall;
          }
        }
      }
    }
  }
  return recall;
}
  
void PerfectPhyloTree::initArcs(const PerfectPhyloGraph& G,
                                PerfectPhyloGraph::Node v_ci,
                                PerfectPhyloGraph::BoolNodeMap& visited)
{
  assert(!visited[v_ci]);
  
  const IntPair& ci = G.nodeToCharState(v_ci);
  
  visited[v_ci] = true;
  
  typedef std::vector<IntPair> IntPairVector;
  typedef std::vector<IntPairVector> IntPairMatrix;
  
  for (PerfectPhyloGraph::IncEdgeIt e(G.G(), v_ci); e != lemon::INVALID; ++e)
  {
    PerfectPhyloGraph::Node v_dj = G.G().oppositeNode(v_ci, e);
    if (visited[v_dj])
      continue;
    
    const IntPair& dj = G.nodeToCharState(v_dj);

    // find lowest ancestor v_dl of v_dj
    Node v_dl = _charStateToNode[ci.first][ci.second];
    while (v_dl != _root)
    {
      if (_nodeToCharState[v_dl].first == dj.first)
      {
        break;
      }
      else
      {
        v_dl = _T.source(InArcIt(_T, v_dl));
      }
    }
    
    if (v_dl == _root)
    {
      if (_S[dj.first].isParent(0, dj.second))
      {
        _T.addArc(_charStateToNode[ci.first][ci.second], _charStateToNode[dj.first][dj.second]);
        initArcs(G, v_dj, visited);
      }
    }
    else
    {
      const IntPair& dl = _nodeToCharState[v_dl];
      if (_S[dl.first].isParent(dl.second, dj.second))
      {
        _T.addArc(_charStateToNode[ci.first][ci.second], _charStateToNode[dj.first][dj.second]);
        initArcs(G, v_dj, visited);
      }
    }
  }
}

void PerfectPhyloTree::writeDOT(const RealTensor& F,
                                const RealMatrix& U,
                                const StringToStringMap& label2color,
                                std::ostream& out) const
{
  out << "digraph T {" << std::endl;
  
  // collapse, assign representative
  typedef Digraph::NodeMap<Node> NodeNodeMap;
  NodeNodeMap representative(_T, lemon::INVALID);
  BoolNodeMap used(_T, false);
  
  for (NodeIt v_ci(_T); v_ci != lemon::INVALID; ++v_ci)
  {
    const IntPair& ci =_nodeToCharState[v_ci];
    const int row_idx = ci.first == 0 && ci.second == 0 ? 0 : F.n() * (ci.second - 1) + ci.first + 1;
    
    for (int p = 0; p < F.m(); ++p)
    {
      if (U(p, row_idx) >= 0.05)
      {
        used[v_ci] = true;
        representative[v_ci] = v_ci;
        break;
      }
    }
    
    if (lemon::countOutArcs(_T, v_ci) > 1)
    {
      representative[v_ci] = v_ci;
    }
  }
  
  for (NodeIt v_ci(_T); v_ci != lemon::INVALID; ++v_ci)
  {
    if (v_ci == _root)
      continue;
    
    if (representative[v_ci] != lemon::INVALID)
    {
      // assign v_ci as the representative of all direct ancestors, with out-deg 1
      Node v_dj = _T.source(InArcIt(_T, v_ci));
      while (representative[v_dj] == lemon::INVALID)
      {
        representative[v_dj] = v_ci;
        v_dj = _T.source(InArcIt(_T, v_dj));
      }
    }
  }
  
  for (NodeIt v_ci(_T); v_ci != lemon::INVALID; ++v_ci)
  {
    if (representative[v_ci] != v_ci)
      continue;
    
    StringVector label_v_ci;
    bool mutation = false;
    std::string color = "black";
    Node v_dj = v_ci;
    while (representative[v_dj] == v_ci)
    {
      if (v_dj == _root)
      {
        label_v_ci.push_back("*  (1,1,0)");
        mutation = true;
        break;
      }
      else
      {
        const IntPair& dj =_nodeToCharState[v_dj];
        if (!isMutationVertex(v_dj))
        {
          label_v_ci.push_back("<font color=\"maroon\">" + F.getColLabel(dj.first) + "  " + _S[dj.first].label(dj.second) + "</font>");
        }
        else
        {
          label_v_ci.push_back(F.getColLabel(dj.first) + "  " + _S[dj.first].label(dj.second));
          mutation = true;
          auto label_it = label2color.find(F.getColLabel(dj.first));
          if (label_it != label2color.end())
          {
            color = label_it->second;
          }
        }
      }
      v_dj = _T.source(InArcIt(_T, v_dj));
    }
    
    out << "\t" << _T.id(v_ci) << " [label=<";
    bool first = true;
    for (const std::string& label : boost::adaptors::reverse(label_v_ci))
    {
      if (first)
        first = false;
      else
        out << "<BR/>";
      
      out << label;
    }
    out << ">";
    
    if (mutation)
    {
      out << ",penwidth=5,color=" << color;
    }
    else
    {
      out << ",penwidth=3,style=dashed";
    }

//    const IntPair& ci =_nodeToCharState[v_ci];
//    if (v_ci == _root)
//    {
//      out << "\t" << _T.id(v_ci) << " [label=\"*\\n" << _S[ci.first].label(ci.second) << "\",penwidth=5";
//    }
//    else
//    {
//      out << "\t" << _T.id(v_ci) << " [label=\"" << F.getColLabel(ci.first) << "  " << _S[ci.first].label(ci.second) << "\"";
//      
//      if (!isMutationVertex(v_ci))
//      {
//        out << ",fontcolor=maroon,penwidth=3,style=dashed";
//      }
//      else
//      {
//        auto label_it = label2color.find(F.getColLabel(ci.first));
//        if (label_it != label2color.end())
//        {
//          out << ",penwidth=5,color=" << label_it->second;
//        }
//      }
//    }
//    
//    if (!used[v_ci])
//    {
//      out << ",shape=box";
//    }
    
    out << "]" << std::endl;
  }
  
  for (ArcIt a_cidj(_T); a_cidj != lemon::INVALID; ++a_cidj)
  {
    Node v_ci = _T.source(a_cidj);
    if (representative[v_ci] != v_ci)
      continue;
    
    Node v_dj = _T.target(a_cidj);
    out << "\t" << _T.id(v_ci) << " -> " << _T.id(representative[v_dj]) << std::endl;
//    out << "\t" << _T.id(v_ci) << " -> " << _T.id(v_dj) << std::endl;
  }
  
  out << "}" << std::endl;
}
  
void PerfectPhyloTree::writeDOT(std::ostream& out) const
{
  out << "digraph T {" << std::endl;
  
  for (NodeIt v_ci(_T); v_ci != lemon::INVALID; ++v_ci)
  {
    const IntPair& ci =_nodeToCharState[v_ci];
    if (v_ci == _root)
    {
      out << "\t" << _T.id(v_ci) << " [label=\"*\\n" << _S[ci.first].label(ci.second) << "\",penwidth=5";
    }
    else
    {
      out << "\t" << _T.id(v_ci) << " [label=\"" << ci.first << "\\n" << _S[ci.first].label(ci.second) << "\"";
      
      if (!isMutationVertex(v_ci))
      {
        out << ",color=red";
      }
    }
    
    out << "]" << std::endl;
  }
  
  for (ArcIt a_cidj(_T); a_cidj != lemon::INVALID; ++a_cidj)
  {
    Node v_ci = _T.source(a_cidj);
    Node v_dj = _T.target(a_cidj);
    out << "\t" << _T.id(v_ci) << " -> " << _T.id(v_dj) << std::endl;
  }
  
  out << "}" << std::endl;
}
  
std::ostream& operator<<(std::ostream& out, const PerfectPhyloTree& T)
{
  out << T.A();
  for (int c = 0; c < T.A().n(); ++c)
  {
    out << T.S(c);
  }
  return out;
}
  
std::istream& operator>>(std::istream& in, PerfectPhyloTree& T)
{
  PerfectPhyloMatrix A;
  in >> A;
  
  // who's owning S? T is.
  return in;
}
  
} // namespace gm
