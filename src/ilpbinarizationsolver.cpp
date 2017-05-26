/*
 * ilpbinarizationsolver.cpp
 *
 *  Created on: 3-feb-2016
 *      Author: M. El-Kebir
 */

#include "ilpbinarizationsolver.h"

IlpBinarizationSolver::IlpBinarizationSolver(const NonBinaryCloneTree& T,
                                             const std::string& primary,
                                             MigrationGraph::Pattern pattern,
                                             const std::string& gurobiLogFilename,
                                             const StringPairList& forcedComigrations)
  : IlpSolver(T, primary, pattern, gurobiLogFilename, forcedComigrations)
  , _G()
  , _label(_G)
  , _TtoG(_T.tree(), lemon::INVALID)
  , _GtoT(_G, lemon::INVALID)
  , _w()
  , _pResCloneTree(NULL)
  , _pResLPlus(NULL)
{
}

void IlpBinarizationSolver::initIndices()
{
  // init _G
  constructG(_T.tree(), _T.root());

  IlpSolver::initIndices();
}

void IlpBinarizationSolver::initVariables()
{
  IlpSolver::initVariables();
  
  const Digraph& G = tree();
  
  const int nrArcs = _indexToArc.size();
  
  char buf[1024];
  
  _w = VarArray(nrArcs);
  for (int ij = 0; ij < nrArcs; ++ij)
  {
    Arc a_ij = _indexToArc[ij];
    Node v_i = G.source(a_ij);
    Node v_j = G.target(a_ij);
    
    snprintf(buf, 1024, "w_%s_%s",
             label(v_i).c_str(),
             label(v_j).c_str());
    _w[ij] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
  }
  
  _model.update();
}

void IlpBinarizationSolver::initConstraints()
{
  IlpSolver::initConstraints();
  
  const Digraph& G = tree();
  const int nrArcs = _indexToArc.size();
  
  for (int ij = 0; ij < nrArcs; ++ij)
  {
    _model.addConstr(_y[ij] - _w[ij] <= 0);
  }
  
  GRBLinExpr sum;
  
  // unique parent: in-degree 1
  for (NodeIt v_j(G); v_j != lemon::INVALID; ++v_j)
  {
    if (v_j == root())
    {
      continue;
    }
    
    for (InArcIt a_ij(G, v_j); a_ij != lemon::INVALID; ++a_ij)
    {
      int ij = (*_pArcToIndex)[a_ij];
      sum += _w[ij];
    }
    
    _model.addConstr(sum == 1);
    sum.clear();
  }
  
  // out-degree 2
  for (NodeIt v_i(G); v_i != lemon::INVALID; ++v_i)
  {
    if (isLeaf(v_i))
    {
      continue;
    }
    
    if (lemon::countOutArcs(G, v_i) == 1)
    {
      continue;
    }
    
    for (OutArcIt a_ij(G, v_i); a_ij != lemon::INVALID; ++a_ij)
    {
      int ij = (*_pArcToIndex)[a_ij];
      sum += _w[ij];
    }
    
    _model.addConstr(sum == 2);
    sum.clear();
  }
}

void IlpBinarizationSolver::processSolution()
{
  // get vertex labeling
  IlpSolver::processSolution();
  
  // get binarization
  BoolArcMap binarization(_G, false);
  
  for (ArcIt a_ij(_G); a_ij != lemon::INVALID; ++a_ij)
  {
    int ij = (*_pArcToIndex)[a_ij];
    if (_w[ij].get(GRB_DoubleAttr_X) >= 0.4)
    {
      binarization[a_ij] = true;
    }
  }
  
  lemon::FilterArcs<Digraph> TT(_G, binarization);
  Digraph resT;
  StringNodeMap resLPlus(resT);
  StringNodeMap resLabel(resT);
  Node resRoot;
  lemon::digraphCopy(TT, resT)
    .nodeMap(*_pLPlus, resLPlus)
    .nodeMap(_label, resLabel)
    .node(root(), resRoot)
    .run();
  
  _pResCloneTree = new CloneTree(resT,
                                 resRoot,
                                 resLabel,
                                 resLPlus);
  
  _pResLPlus = new StringNodeMap(_pResCloneTree->tree());
  for (NodeIt v(resT); v != lemon::INVALID; ++v)
  {
    Node new_v = _pResCloneTree->getNodeByLabel(resLabel[v]);
    _pResLPlus->set(new_v, resLPlus[v]);
  }
}

void IlpBinarizationSolver::constructG(const Digraph& T,
                                       Node v)
{
  char buf[1024];
  
  int out_deg = lemon::countOutArcs(T, v);
  if (out_deg > 2)
  {
    NodeVector nodes;
    for (int i = 0; i < out_deg - 1; ++i)
    {
      Node vv = _G.addNode();
      _GtoT[vv] = v;
      
      snprintf(buf, 1024, "%s_%d", _T.label(v).c_str(), i);
      _label[vv] = buf;
      
      for (int j = 0; j < i; ++j)
      {
        _G.addArc(nodes[j], vv);
      }
      if (i == 0)
      {
        _TtoG[v] = vv;
      }
      
      nodes.push_back(vv);
    }
    
    for (OutArcIt a(T, v); a != lemon::INVALID; ++a)
    {
      Node w = T.target(a);
      constructG(T, w);
      
      for (int i = 0; i < out_deg - 1; ++i)
      {
        _G.addArc(nodes[i], _TtoG[w]);
      }
    }
  }
  else
  {
    Node vv = _G.addNode();
    _label[vv] = _T.label(v);
    _TtoG[v] = vv;
    _GtoT[vv] = v;

    for (OutArcIt a(T, v); a != lemon::INVALID; ++a)
    {
      Node w = T.target(a);
      constructG(T, w);
      
      _G.addArc(vv, _TtoG[w]);
    }
  }
}

void IlpBinarizationSolver::writeDOT(std::ostream& out,
                                     const StringToIntMap& colorMap) const
{
  out << "digraph G {" << std::endl;
  
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt uu(_G); uu != lemon::INVALID; ++uu)
  {
    Node u = _GtoT[uu];
    if (OutArcIt(_G, uu) == lemon::INVALID)
    {
      out << "\t\t" << _G.id(uu) << " [penwidth=3,colorscheme=set19,color="
          << colorMap.find(_T.l(u))->second << ",label=\"" << _label[uu]
          << "\\n" << _T.l(u) << "\"]" << std::endl;
    }
  }
  out << "\t}" << std::endl;
  
  for (NodeIt uu(_G); uu != lemon::INVALID; ++uu)
  {
    if (OutArcIt(_G, uu) != lemon::INVALID)
    {
      out << "\t" << _G.id(uu) << " [label=\"" << _label[uu] << "\"]"
          << std::endl;
    }
  }
  
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    out << "\t" << _G.id(_G.source(a)) << " -> " << _G.id(_G.target(a))
        << std::endl;
  }
  
  out << "}" << std::endl;
}
