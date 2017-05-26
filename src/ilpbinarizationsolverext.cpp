/*
 * ilpbinarizationsolverext.cpp
 *
 *  Created on: 15-apr-2017
 *      Author: M. El-Kebir
 */

#include "ilpbinarizationsolverext.h"

IlpBinarizationSolverExt::IlpBinarizationSolverExt(const NonBinaryCloneTree& T,
                                                   const FrequencyMatrix& F,
                                                   const std::string& primary,
                                                   Mode mode,
                                                   const std::string& gurobiLogFilename,
                                                   const StringPairList& forcedComigrations)
  : IlpSolverExt(T, F, primary, mode, gurobiLogFilename, forcedComigrations)
  , _isSampleNode(_G, false)
  , _sampleNodeToArcMap(_G, lemon::INVALID)
  , _w()
{
}

void IlpBinarizationSolverExt::initVariables()
{
  IlpSolverExt::initVariables();

  const int nrArcs = _indexToArc.size();
  
  char buf[1024];
  
  _w = VarArray(nrArcs);
  for (int ij = 0; ij < nrArcs; ++ij)
  {
    Arc a_ij = _indexToArc[ij];
    Node v_i = _G.source(a_ij);
    Node v_j = _G.target(a_ij);
    
    snprintf(buf, 1024, "w_%s_%s",
             label(v_i).c_str(),
             label(v_j).c_str());
    _w[ij] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
  }
  
  _model.update();
}

void IlpBinarizationSolverExt::initConstraints()
{
  IlpSolverExt::initConstraints();
  
  const int nrArcs = _indexToArc.size();
  const int nrSamples = _indexToSample.size();
  
  for (int ij = 0; ij < nrArcs; ++ij)
  {
    _model.addConstr(_y[ij] - _w[ij] <= 0);
  }
  
  GRBLinExpr sum, sum2;
  
  // unique parent: in-degree 1
  for (NodeIt v_j(_G); v_j != lemon::INVALID; ++v_j)
  {
    if (v_j == root())
    {
      continue;
    }
    
    for (InArcIt a_ij(_G, v_j); a_ij != lemon::INVALID; ++a_ij)
    {
      int ij = _arcToIndex[a_ij];
      sum += _w[ij];
    }
    
    if (_isSampleNode[v_j])
    {
      _model.addConstr(sum == 1);
      sum.clear();
    }
  }
  
  // out-degree 2
  for (NodeIt v_i(_G); v_i != lemon::INVALID; ++v_i)
  {
    if (isLeaf(v_i))
    {
      continue;
    }
    
    if (lemon::countOutArcs(_G, v_i) == 1)
    {
      continue;
    }
    
    for (OutArcIt a_ij(_G, v_i); a_ij != lemon::INVALID; ++a_ij)
    {
      int ij = _arcToIndex[a_ij];
      sum += _w[ij];
    }

    _model.addConstr(sum == 2);
    sum.clear();
  }
  
  // sum of out-degrees is #enabled arcs (do this iteratively)
  // move disabled arcs to last node
  // if leaf is unused, disable corresponding arc.
  // assign disabled arcs last node of the expansion
}

void IlpBinarizationSolverExt::processSolution()
{
  const int nrSamples = _indexToSample.size();
  const int nrCharacters = _F.getNrCharacters();
  
  for (NodeIt v_i(_G); v_i != lemon::INVALID; ++v_i)
  {
    int i = _nodeToIndex[v_i];
    bool found = false;
    for (int s = 0; s < nrSamples; ++s)
    {
      if (_x[i][s].get(GRB_DoubleAttr_X) >= 0.4)
      {
        _lPlus[v_i] = _indexToSample[s];
        found = true;
        break;
      }
    }
  }
  
  for (NodeIt v_i(_T.tree()); v_i != lemon::INVALID; ++v_i)
  {
    int mapped_i = _F.characterToIndex(_T.label(v_i));
    for (int s = 0; s < nrSamples; ++s)
    {
      // set U and F
      assert(0 <= mapped_i && mapped_i < nrCharacters);
      _resU[v_i][s] = _u[s][mapped_i].get(GRB_DoubleAttr_X);
      _resF[v_i][s] = _f[s][mapped_i].get(GRB_DoubleAttr_X);
    }
  }
  
  // remove nodes and arcs that are unused
  BoolNodeMap filterNodeMap(_G, true);
  BoolArcMap filterArcMap(_G, true);
  for (ArcIt a_ij(_G); a_ij != lemon::INVALID; ++a_ij)
  {
    Node v_i = _G.source(a_ij);
    Node v_j = _G.target(a_ij);
    int ij = _arcToIndex[a_ij];
    
    if (isLeaf(v_j))
    {
      int s = _sampleToIndex[_l[v_j]];
      if (_resU[_GtoT[v_i]][s] == 0)
      {
        filterNodeMap[v_j] = false;
        filterArcMap[a_ij] = false;
      }
    }
    if (_w[ij].get(GRB_DoubleAttr_X) < 0.4)
    {
        filterArcMap[a_ij] = false;
    }
  }
  
  lemon::SubDigraph<Digraph> TT(_G, filterNodeMap, filterArcMap);
  Digraph resT;
  StringNodeMap resLPlus(resT);
  StringNodeMap resLabel(resT);
  IntNodeMap resNodeToIndex(resT);
  Node resRoot;
  lemon::digraphCopy(TT, resT)
    .nodeMap(_lPlus, resLPlus)
    .nodeMap(_label, resLabel)
    .nodeMap(_nodeToIndex, resNodeToIndex)
    .node(root(), resRoot)
    .run();
  
  IntVector characterCount(_F.getNrCharacters(), 0);
  for (NodeIt v_i(resT); v_i != lemon::INVALID; ++v_i)
  {
    int mapped_i = getCharacterIndex(resT, resNodeToIndex, v_i);
    ++characterCount[mapped_i];
  }
  
  // PRE: unused sample nodes have been removed
  // delete a vertex if it has out-degree 1
  bool done = false;
  do
  {
    bool changed = false;
    for (NodeIt v_i(resT); v_i != lemon::INVALID && !changed; ++v_i)
    {
      int mapped_i = getCharacterIndex(resT, resNodeToIndex, v_i);
      if (characterCount[mapped_i] == 1)
        continue;
      
      int i = resNodeToIndex[v_i];
      Node org_v_i = _indexToNode[i];
      if (OutArcIt(resT, v_i) == lemon::INVALID && _GtoT[org_v_i] != lemon::INVALID)
      {
        // out deg 0 and not a (used) sample node
        resT.erase(v_i);
        --characterCount[mapped_i];
        changed = true;
      }
      else if (lemon::countOutArcs(resT, v_i) == 1)
      {
        // out deg 1
        Node v_j = resT.target(OutArcIt(resT, v_i));
        if (v_i == resRoot && resLPlus[v_i] == resLPlus[v_j])
        {
          // v_i is root vertex
          resT.erase(v_i);
          --characterCount[mapped_i];
          resRoot = v_j;
          changed = true;
        }
        else if (v_i != resRoot)
        {
          // v_i is not the root
          Node v_pi_i = resT.source(InArcIt(resT, v_i));
          if (resLPlus[v_i] == resLPlus[v_j]
              || resLPlus[v_i] == resLPlus[v_pi_i])
          {
            resT.erase(v_i);
            --characterCount[mapped_i];
            resT.addArc(v_pi_i, v_j);
            changed = true;
          }
        }
      }
    }
    if (!changed)
      done = true;
  } while (!done);
  
  // infer edge labeling
  IntNodeMap resCharacterLabel(resT, -1);
  BoolVector resVisited(_F.getNrCharacters(), false);
  labelEdges(resT,
             resNodeToIndex,
             resRoot,
             resCharacterLabel,
             resVisited);
  
  _pResCloneTree = new NonBinaryCloneTree(resT, resRoot, resLabel, resLPlus);
  _pResLPlus = new StringNodeMap(_pResCloneTree->tree());
  _pResF = new DoubleVectorNodeMap(_pResCloneTree->tree());
  _pResU = new DoubleNodeMap(_pResCloneTree->tree());
  _pResCharacterLabel = new IntNodeMap(_pResCloneTree->tree());
  for (NodeIt v(resT); v != lemon::INVALID; ++v)
  {
    Node new_v = _pResCloneTree->getNodeByLabel(resLabel[v]);
    _pResLPlus->set(new_v, resLPlus[v]);
    Node org_v = _indexToNode[resNodeToIndex[v]];
    _pResCharacterLabel->set(new_v, resCharacterLabel[v]);
    
    Node v_in_T = _GtoT[org_v];
    if (v_in_T == lemon::INVALID)
    {
      Node parent_org_v = _G.source(InArcIt(_G, org_v));
      v_in_T = _GtoT[parent_org_v];
    }
    assert(v_in_T != lemon::INVALID);
    _pResF->set(new_v, _resF[v_in_T]);
    
    if (_isSampleNode[org_v])
    {
      Node parent_org_v = _G.source(InArcIt(_G, org_v));
//      int i = _nodeToIndex[parent_org_v];
      int s = _sampleToIndex[_l[org_v]];
      
      _pResU->set(new_v, _resU[_GtoT[parent_org_v]][s]);
    }
  }
  
  assert(_lPlus[root()] == _primary);
}

void IlpBinarizationSolverExt::labelEdges(const Digraph& resT,
                                          const IntNodeMap& resNodeToIndex,
                                          Node v_i,
                                          IntNodeMap& characterLabel,
                                          BoolVector& visited)
{
  Node org_v_i = _indexToNode[resNodeToIndex[v_i]];
  Node tree_v_i = _GtoT[org_v_i];
  if (tree_v_i == lemon::INVALID)
  {
    org_v_i = _G.source(InArcIt(_G, org_v_i));
    tree_v_i = _GtoT[org_v_i];
  }
  int characterIndex = _F.characterToIndex(_T.label(tree_v_i));
  if (!visited[characterIndex])
  {
    visited[characterIndex] = true;
    characterLabel[v_i] = characterIndex;
  }
  
  for (OutArcIt a(resT, v_i); a != lemon::INVALID; ++a)
  {
    Node v_j = resT.target(a);
    labelEdges(resT, resNodeToIndex, v_j, characterLabel, visited);
  }
}

void IlpBinarizationSolverExt::constructG()
{
  constructG(_T.root());
}

void IlpBinarizationSolverExt::constructG(Node v)
{
  const int m = _F.getNrSamples();
  const Digraph& T = _T.tree();
  
  char buf[1024];

  int out_deg = lemon::countOutArcs(T, v);
  NodeVector sampleNodes;
  // attach sample s to vertex i if Umax[s][i] > 0
  for (int s = 0; s < m; ++s)
  {
    const std::string& sStr = _F.indexToSample(s);
    int mapped_i = _F.characterToIndex(_T.label(v));
    if ((_T.isLeaf(v) && _Fmin[s][mapped_i] > 0) ||
        (!_T.isLeaf(v) && _Umax[s][mapped_i] > 0))
    {
      Node v_is = _G.addNode();
      _isSampleNode[v_is] = true;
      _GtoT[v_is] = lemon::INVALID;
      _label[v_is] = _T.label(v) + "_" + sStr;
      _l[v_is] = sStr;
      _nodeToIndex[v_is] = _indexToNode.size();
      _indexToNode.push_back(v_is);
      sampleNodes.push_back(v_is);
    }
  }
  
  out_deg += sampleNodes.size();
  if (out_deg > 2)
  {
    NodeVector nodes;
    for (int i = 0; i < out_deg - 1; ++i)
    {
      Node vv = _G.addNode();
      _isSampleNode[vv] = false;
      _sampleNodeToArcMap[vv] = lemon::INVALID;
      _GtoT[vv] = v;
      _nodeToIndex[vv] = _indexToNode.size();
      _indexToNode.push_back(vv);
      
      snprintf(buf, 1024, "%s_%d", _T.label(v).c_str(), i);
      _label[vv] = buf;
      
      for (int j = 0; j < i; ++j)
      {
        _G.addArc(nodes[j], vv);
      }
      if (i == 0)
      {
        _TtoG[v] = vv;
        if (v == _T.root())
        {
          _root = vv;
        }
      }
      
      nodes.push_back(vv);
    }
    
    for (OutArcIt a(T, v); a != lemon::INVALID; ++a)
    {
      Node w = T.target(a);
      constructG(w);
      
      for (int i = 0; i < nodes.size(); ++i)
      {
        _G.addArc(nodes[i], _TtoG[w]);
      }
    }
    
    // attach samples
    int is = 0;
    for (Node w : sampleNodes)
    {
      for (int i = 0; i < nodes.size(); ++i)
      {
        _G.addArc(nodes[i], w);
      }
      ++is;
    }
  }
  else
  {
    Node vv = _G.addNode();
    _isSampleNode[vv] = false;
    _sampleNodeToArcMap[vv] = lemon::INVALID;
    _label[vv] = _T.label(v);
    _TtoG[v] = vv;
    _GtoT[vv] = v;
    _nodeToIndex[vv] = _indexToNode.size();
    _indexToNode.push_back(vv);
    if (v == _T.root())
    {
      _root = vv;
    }
    
    for (OutArcIt a(T, v); a != lemon::INVALID; ++a)
    {
      Node w = T.target(a);
      constructG(w);
      
      _G.addArc(vv, _TtoG[w]);
    }
    
    // attach samples
    for (Node w : sampleNodes)
    {
      Arc a = _G.addArc(vv, w);
      _sampleNodeToArcMap[w] = a;
    }
  }
}
