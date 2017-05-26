/*
 * ilpsolverext.cpp
 *
 *  Created on: 13-apr-2017
 *      Author: M. El-Kebir
 */

#include "ilpsolverext.h"
#include <lemon/time_measure.h>

IlpSolverExt::IlpSolverExt(const NonBinaryCloneTree& T,
                           const FrequencyMatrix& F,
                           const std::string& primary,
                           MigrationGraph::Pattern pattern,
                           const std::string& gurobiLogFilename,
                           const StringPairList& forcedComigrations)
  : IlpSolver(T, primary, pattern, gurobiLogFilename, forcedComigrations)
  , _F(F)
  , _Fmin(DoubleMatrix(F.getNrSamples(), DoubleVector(F.getNrCharacters(), 0)))
  , _Fmax(DoubleMatrix(F.getNrSamples(), DoubleVector(F.getNrCharacters(), 0)))
  , _Umax(DoubleMatrix(F.getNrSamples(), DoubleVector(F.getNrCharacters(), 0)))
  , _G()
  , _label(_G)
  , _l(_G)
  , _lPlus(_G)
  , _TtoG(_T.tree(), lemon::INVALID)
  , _GtoT(_G, lemon::INVALID)
  , _isSampleNode(_G, false)
  , _f()
  , _u()
  , _w()
  , _resF(_T.tree(), DoubleVector(F.getNrSamples(), 0))
  , _resU(_T.tree(), DoubleVector(F.getNrSamples(), 0))
  , _pResCloneTree(NULL)
  , _pResLPlus(NULL)
  , _pResF(NULL)
  , _pResU(NULL)
  , _pResCharacterLabel(NULL)
{
  _pArcToIndex = new IntArcMap(_G);
  _pNodeToIndex = new IntNodeMap(_G);
}

IlpSolverExt::~IlpSolverExt()
{
  delete _pResLPlus;
  delete _pResF;
  delete _pResU;
  delete _pResCharacterLabel;
  delete _pResCloneTree;
}

void IlpSolverExt::init(double upperBound)
{
  computeFmin(_T.root());
  computeFmax(_T.root());
  computeUmax();
  
  IlpSolver::init(upperBound);
}

void IlpSolverExt::constructG()
{
  constructG(_T.root());
}

void IlpSolverExt::computeFmin(Node v_i)
{
  const int m = _F.getNrSamples();
  const int i = _F.characterToIndex(_T.label(v_i));
  const Digraph& T = _T.tree();
  
  if (OutArcIt(T, v_i) == lemon::INVALID)
  {
    for (int s = 0; s < m; ++s)
    {
      _Fmin[s][i] = _F.min(s,i);
    }
  }
  else
  {
    for (OutArcIt a_ij(T, v_i); a_ij != lemon::INVALID; ++a_ij)
    {
      Node v_j = T.target(a_ij);
      computeFmin(v_j);
    }

    for (int s = 0; s < m; ++s)
    {
      double sum_s = 0;
      for (OutArcIt a_ij(T, v_i); a_ij != lemon::INVALID; ++a_ij)
      {
        Node v_j = T.target(a_ij);
        int j = _F.characterToIndex(_T.label(v_j));
        sum_s += _Fmin[s][j];
      }
      _Fmin[s][i] = std::max(_F.min(s, i), sum_s);
    }
  }
}

void IlpSolverExt::constructG(Node v)
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
      _pNodeToIndex->set(v_is, _indexToNode.size());
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
      _GtoT[vv] = v;
      _pNodeToIndex->set(vv, _indexToNode.size());
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
    _label[vv] = _T.label(v);
    _TtoG[v] = vv;
    _GtoT[vv] = v;
    _pNodeToIndex->set(vv, _indexToNode.size());
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
      _G.addArc(vv, w);
    }
  }
}

void IlpSolverExt::computeFmax(Node v_i)
{
  const int m = _F.getNrSamples();
  const int i = _F.characterToIndex(_T.label(v_i));
  const Digraph& T = _T.tree();

  if (v_i == _T.root())
  {
    for (int s = 0; s < m; ++s)
    {
      _Fmax[s][i] = _F.max(s,i);
    }
  }
  else
  {
    Node v_pi_i = T.source(InArcIt(T, v_i));
    int pi_i = _F.characterToIndex(_T.label(v_pi_i));
    
    for (int s = 0; s < m; ++s)
    {
      double sum_s = _Fmax[s][pi_i];
      for (OutArcIt a_pi_i_j(T, v_pi_i); a_pi_i_j != lemon::INVALID; ++a_pi_i_j)
      {
        Node v_j = T.target(a_pi_i_j);
        if (v_j == v_i) continue;

        int j = _F.characterToIndex(_T.label(v_j));
        sum_s -= _Fmin[s][j];
      }
      
      _Fmax[s][i] = std::min(_F.max(s, i), sum_s);
      assert(_Fmax[s][i] >= 0);
    }
  }
  
  for (OutArcIt a_ij(T, v_i); a_ij != lemon::INVALID; ++a_ij)
  {
    Node v_j = T.target(a_ij);
    computeFmax(v_j);
  }
}

void IlpSolverExt::computeUmax()
{
  const int m = _F.getNrSamples();
  const int n = _F.getNrCharacters();
  const Digraph& T = _T.tree();
  
  // compute Umax
  for (int s = 0; s < m; ++s)
  {
    for (int i = 0; i < n; ++i)
    {
      Node v_i = _T.getNodeByLabel(_F.indexToCharacter(i));
      if (_T.isLeaf(v_i) && _Fmin[s][i] == 0)
      {
        _Umax[s][i] = 0;
      }
      else
      {
        _Umax[s][i] = _Fmax[s][i];
        for (OutArcIt a_ij(T, v_i); a_ij != lemon::INVALID; ++a_ij)
        {
          Node v_j = T.target(a_ij);
          int j = _F.characterToIndex(_T.label(v_j));
          
          _Umax[s][i] -= _Fmin[s][j];
        }
      }
    }
  }
}

void IlpSolverExt::writeDOT(std::ostream& out) const
{
  const int m = _F.getNrSamples();
  
  out << "digraph T {" << std::endl;
  
  for (NodeIt v_i(_G); v_i != lemon::INVALID; ++v_i)
  {
    const int i = (*_pNodeToIndex)[v_i];
    out << "\t" << i << " [label=\"" << _label[v_i];
    if (!isLeaf(v_i))
    {
      int mapped_i = _F.characterToIndex(_T.label(_GtoT[v_i]));
      for (int s = 0; s < m; ++s)
      {
        out << "\n" << _F.indexToSample(s)
            << " : [" << _F.min(s, mapped_i) << "," << _F.max(s, mapped_i) << "]"
            << " ; [" << _Fmin[s][mapped_i] << "," << _Fmax[s][mapped_i] << "]"
            << " ; " << _Umax[s][mapped_i];
      }
    }
    out <<  "\"]" << std::endl;
  }
  
  for (ArcIt a_ij(_G); a_ij != lemon::INVALID; ++a_ij)
  {
    Node v_i = _G.source(a_ij);
    Node v_j = _G.target(a_ij);
    const int i = (*_pNodeToIndex)[v_i];
    const int j = (*_pNodeToIndex)[v_j];
    
    out << "\t" << i << " -> " << j << std::endl;
  }
  
  out << "}" << std::endl;
}

void IlpSolverExt::writeDOT(std::ostream& out,
                            const StringToIntMap& colorMap) const
{
  const int m = _F.getNrSamples();
  
  out << "digraph T {" << std::endl;
  
  for (NodeIt v_i(_G); v_i != lemon::INVALID; ++v_i)
  {
    const int i = (*_pNodeToIndex)[v_i];
    out << "\t" << i << " [label=\"" << _label[v_i];
    if (!isLeaf(v_i))
    {
      int mapped_i = _F.characterToIndex(_T.label(_GtoT[v_i]));
      for (int s = 0; s < m; ++s)
      {
        out << "\n" << _F.indexToSample(s)
            << " : [" << _F.min(s, mapped_i) << "," << _F.max(s, mapped_i) << "]"
            << " ; [" << _Fmin[s][mapped_i] << "," << _Fmax[s][mapped_i] << "]"
            << " ; " << _Umax[s][mapped_i] << " ; " << _resF[_GtoT[v_i]][s] << " ; " << _resU[_GtoT[v_i]][s];
      }
    }
    out << "\\n" << _lPlus[v_i] << "\""
        << ",penwidth=3,colorscheme=set19,color="
        << colorMap.find(_lPlus[v_i])->second << "]" << std::endl;
  }
  
  for (ArcIt a_ij(_G); a_ij != lemon::INVALID; ++a_ij)
  {
    Node v_i = _G.source(a_ij);
    Node v_j = _G.target(a_ij);
    const int i = (*_pNodeToIndex)[v_i];
    const int j = (*_pNodeToIndex)[v_j];
    
    out << "\t" << i << " -> " << j << std::endl;
  }
  
  out << "}" << std::endl;
}

void IlpSolverExt::initIndices()
{
  const int m = _F.getNrSamples();
  _indexToSample = StringVector(m);
  
  // construct G
  constructG();
  
  _indexToArc.clear();
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    _pArcToIndex->set(a, _indexToArc.size());
    _indexToArc.push_back(a);
  }

  _primaryIndex = -1;
  _sampleToIndex.clear();
  for (int s = 0; s < m; ++s)
  {
    _sampleToIndex[_F.indexToSample(s)] = s;
    _indexToSample[s] = _F.indexToSample(s);
    if (_primary == _F.indexToSample(s))
    {
      _primaryIndex = s;
    }
  }
  
  assert(_primaryIndex != -1);
}

void IlpSolverExt::initVariables()
{
  const int nrArcs = _indexToArc.size();
  const int nrNodes = _indexToNode.size();
  const int nrSamples = _indexToSample.size();
  const int nrCharacters = _F.getNrCharacters();
  
  /// x[i][s] = 1 iff vertex v_i is labeled by sample s
  char buf[1024];
  _x = VarMatrix(nrNodes);
  for (int i = 0; i < nrNodes; ++i)
  {
    _x[i] = VarArray(nrSamples + 1);
    for (int s = 0; s < nrSamples; ++s)
    {
      snprintf(buf, 1024, "x_%s_%s",
               label(_indexToNode[i]).c_str(),
               _indexToSample[s].c_str());
      _x[i][s] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
    }
    snprintf(buf, 1024, "x_%s_DUMMY", label(_indexToNode[i]).c_str());
    _x[i][nrSamples] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
  }
  
  /// y[ij] = 1 iff edge (v_i, v_j) is a migration edge
  _y = VarArray(nrArcs);
  for (int ij = 0; ij < nrArcs; ++ij)
  {
    Arc a_ij = _indexToArc[ij];
    Node v_i = _G.source(a_ij);
    Node v_j = _G.target(a_ij);
    
    snprintf(buf, 1024, "y_%s_%s",
             label(v_i).c_str(),
             label(v_j).c_str());
    _y[ij] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
  }
  
  /// z[ij][s] = 1 iff v_i and v_j are labeled by sample s
  _z = VarMatrix(nrArcs);
  for (int ij = 0; ij < nrArcs; ++ij)
  {
    Arc a_ij = _indexToArc[ij];
    Node v_i = _G.source(a_ij);
    Node v_j = _G.target(a_ij);
    
    _z[ij] = VarArray(nrSamples+1);
    for (int s = 0; s < nrSamples; ++s)
    {
      snprintf(buf, 1024, "z_%s_%s_%s",
               label(v_i).c_str(),
               label(v_j).c_str(),
               _indexToSample[s].c_str());
      
      _z[ij][s] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
    }
    snprintf(buf, 1024, "z_%s_%s_DUMMY",
             label(v_i).c_str(),
             label(v_j).c_str());
    _z[ij][nrSamples] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
  }
  
  /// c[s][t] = 1 iff there exists a migration edge (v_i, v_j) where l(v_i) = s and l(v_j) = t
  _c = VarMatrix(nrSamples);
  for (int s = 0; s < nrSamples; ++s)
  {
    _c[s] = VarArray(nrSamples);
    for (int t = 0; t < nrSamples; ++t)
    {
      snprintf(buf, 1024, "c_%s_%s",
               _indexToSample[s].c_str(),
               _indexToSample[t].c_str());
      _c[s][t] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
    }
  }
  
  /// d[s] = 1 iff there exists a migration edge (v_i, v_j) where l(v_i) = s
  _d = VarArray(nrSamples);
  for (int s = 0; s < nrSamples; ++s)
  {
    snprintf(buf, 1024, "d_%s",
             _indexToSample[s].c_str());
    _d[s] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
  }
  
  /// f[s][i] = frequency of mutation incoming to clone i in sample s
  _f = VarMatrix(nrSamples);
  for (int s = 0; s < nrSamples; ++s)
  {
    _f[s] = VarArray(nrCharacters);
    for (int c = 0; c < nrCharacters; ++c)
    {
      snprintf(buf, 1024, "f_%s_%s",
               _indexToSample[s].c_str(),
               _F.indexToCharacter(c).c_str());
      _f[s][c] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
    }
  }

  /// u[s][i] = usage of clone i in sample s
  _u = VarMatrix(nrSamples);
  for (int s = 0; s < nrSamples; ++s)
  {
    _u[s] = VarArray(nrCharacters);
    for (int c = 0; c < nrCharacters; ++c)
    {
      snprintf(buf, 1024, "u_%s_%s",
               _indexToSample[s].c_str(),
               _F.indexToCharacter(c).c_str());
      _u[s][c] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
    }
  }
  
  // w[ij] = 1 iff edge (v_i, v_j) is in the tree
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

void IlpSolverExt::initLeafConstraints()
{
  const int nrSamples = _indexToSample.size();
  
  GRBLinExpr sum;
  
  // Leaf color
  for (NodeIt v_i(_G); v_i != lemon::INVALID; ++v_i)
  {
    if (!isLeaf(v_i))
    {
      continue;
    }
    
    int i = (*_pNodeToIndex)[v_i];
    assert(_sampleToIndex.count(_l[v_i]) == 1);
    int s = _sampleToIndex[l(v_i)];
    
    _model.addConstr(_x[i][s] + _x[i][nrSamples] == 1);
    for (int t = 0; t < nrSamples; ++t)
    {
      if (t == s) continue;
      sum += _x[i][t];
    }
    _model.addConstr(sum == 0);
    sum.clear();
  }
}

void IlpSolverExt::initConstraints()
{
  const int nrArcs = _indexToArc.size();
  const int nrCharacters = _F.getNrCharacters();
  const int nrSamples = _indexToSample.size();
  
  GRBLinExpr sum;
  
  // Unique color
  for (NodeIt v_i(_G); v_i != lemon::INVALID; ++v_i)
  {
    int i = (*_pNodeToIndex)[v_i];
    for (int s = 0; s < nrSamples; ++s)
    {
      sum += _x[i][s];
    }
    sum += _x[i][nrSamples];
    _model.addConstr(sum == 1);
    if (!isLeaf(v_i))
    {
      _model.addConstr(_x[i][nrSamples] == 0);
    }
    sum.clear();
  }
  
  // Matching colors
  for (int ij = 0; ij < nrArcs; ++ij)
  {
    Arc a_ij = _indexToArc[ij];
    Node v_i = _G.source(a_ij);
    Node v_j = _G.target(a_ij);
    int i = (*_pNodeToIndex)[v_i];
    int j = (*_pNodeToIndex)[v_j];
    
    for (int s = 0; s < nrSamples + 1; ++s)
    {
      sum += _z[ij][s];
      
      if (s < nrSamples)
      {
        _model.addConstr(_z[ij][s] <= _x[i][s]);
        _model.addConstr(_z[ij][s] <= _x[j][s]);
        //_model.addConstr(_z[ij][s] >= _x[i][s] + _x[j][s] - 1);
      }
      else if (isLeaf(v_j))
      {
        int mapped_i = _F.characterToIndex(_T.label(_GtoT[v_i]));
        assert(0 <= mapped_i && mapped_i < nrCharacters);
        
        int t = _sampleToIndex[_l[v_j]];
        addMatchingColorsAndUsageConstraint(ij, mapped_i, t);
        _model.addConstr(_x[j][t] >= _u[t][mapped_i]);
      }
      else
      {
        _model.addConstr(_z[ij][nrSamples] == 0);
      }
    }
    addMatchingColorsConstraint(sum, ij);
    
    sum.clear();
  }
  
  // Root color
  _model.addConstr(_x[(*_pNodeToIndex)[root()]][_primaryIndex] == 1);

  // Comigration constraints
  for (int s = 0; s < nrSamples; ++s)
  {
    for (int t = 0; t < nrSamples; ++t)
    {
      if (t == s)
      {
        _model.addConstr(_c[s][t] == 0);
      }
      else
      {
        for (int ij = 0; ij < nrArcs; ++ij)
        {
          Arc a_ij = _indexToArc[ij];
          Node v_i = _G.source(a_ij);
          Node v_j = _G.target(a_ij);
          int i = (*_pNodeToIndex)[v_i];
          int j = (*_pNodeToIndex)[v_j];
          
          addComigrationConstraint(s, t, ij, i, j);
        }
      }
    }
  }
  
  if (nrSamples > 1)
  {
    _model.addConstr(_d[_primaryIndex] == 1);
  }
  
  for (int s = 0; s < nrSamples; ++s)
  {
    for (int t = 0; t < nrSamples; ++t)
    {
      _model.addConstr(_d[s] >= _c[s][t]);
    }
  }
  
  for (NodeIt v_i(_T.tree()); v_i != lemon::INVALID; ++v_i)
  {
    int mapped_i = _F.characterToIndex(_T.label(v_i));
    for (int s = 0; s < nrSamples; ++s)
    {
      for (OutArcIt a_ij(_T.tree(), v_i); a_ij != lemon::INVALID; ++a_ij)
      {
        Node v_j = _T.tree().target(a_ij);
        int mapped_j = _F.characterToIndex(_T.label(v_j));
        
        sum += _f[s][mapped_j];
      }
      _model.addConstr(_f[s][mapped_i] >= _Fmin[s][mapped_i]);
      _model.addConstr(_f[s][mapped_i] <= _Fmax[s][mapped_i]);
      _model.addConstr(_u[s][mapped_i] <= _Umax[s][mapped_i]);
      _model.addConstr(_f[s][mapped_i] >= sum);
      _model.addConstr(_u[s][mapped_i] == _f[s][mapped_i] - sum);
      sum.clear();
    }
  }
  
  for (int ij = 0; ij < nrArcs; ++ij)
  {
    _model.addConstr(_y[ij] - _w[ij] <= 0);
  }
  
  // unique parent: in-degree 1
  for (NodeIt v_j(_G); v_j != lemon::INVALID; ++v_j)
  {
    if (v_j == root())
    {
      continue;
    }
    
    for (InArcIt a_ij(_G, v_j); a_ij != lemon::INVALID; ++a_ij)
    {
      int ij = (*_pArcToIndex)[a_ij];
      sum += _w[ij];
    }
    
    _model.addConstr(sum == 1);
    sum.clear();
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
      int ij = (*_pArcToIndex)[a_ij];
      sum += _w[ij];
    }
    
    _model.addConstr(sum == 2);
    sum.clear();
  }
  
  if (_pattern != MigrationGraph::R)
  {
    for (OutArcIt a_ij(_G, root()); a_ij != lemon::INVALID; ++a_ij)
    {
      Node v_j = _G.target(a_ij);
      int j = (*_pNodeToIndex)[v_j];
      sum += _x[j][_primaryIndex];
    }
    _model.addConstr(sum >= 1);
    sum.clear();
    
    for (NodeIt v_j(_G); v_j != lemon::INVALID; ++v_j)
    {
      Node v_i = _G.source(InArcIt(_G, v_j));
      if (isLeaf(v_j))
      {
        int mapped_i = _F.characterToIndex(_T.label(_GtoT[v_i]));
        assert(_sampleToIndex.count(_l[v_j]) == 1);
        int s = _sampleToIndex[l(v_j)];
        if (s == _primaryIndex
            && OutArcIt(_T.tree(), _GtoT[v_i]) == lemon::INVALID
            && _Fmin[s][mapped_i] > 0)
        {
          initPrimaryConstraint(v_j);
        }
      }
    }
  }
  
  if (_pattern == MigrationGraph::S || _pattern == MigrationGraph::PS)
  {
    initSingleSourceSeedingConstraints();
    if (_pattern == MigrationGraph::PS)
    {
      initParallelSingleSourceSeedingConstraints();
    }
  }
  else if (_pattern == MigrationGraph::M)
  {
    initMultiSourceSeedingConstraints();
  }
}

void IlpSolverExt::initObjective(double upperBound)
{
  const int nrArcs = _indexToArc.size();
  const int nrSamples = _indexToSample.size();
  
  GRBLinExpr obj, sum_migrations, sum_comigrations;
  // migrations
  for (int ij = 0; ij < nrArcs; ++ij)
  {
    obj += _y[ij];
    sum_migrations += _y[ij];
  }
  _model.addConstr(sum_migrations >= nrSamples - 1);
  
  if (_forcedComigrations.empty())
  {
    // comigrations
    {
      const double f = 1. / (nrArcs + 1);
      for (int s = 0; s < nrSamples; ++s)
      {
        for (int t = 0; t < nrSamples; ++t)
        {
          obj += f * _c[s][t];
          sum_comigrations += _c[s][t];
        }
      }
    }
    _model.addConstr(sum_comigrations >= nrSamples - 1);
    _model.addConstr(sum_migrations >= sum_comigrations);
    
    // seeding sites
    {
      const double g = (1. / (nrArcs + 1)) * (1. / (nrSamples + 1));
      for (int s = 0; s < nrSamples; ++s)
      {
        obj += g * _d[s];
      }
    }
  }
  else
  {
    // minimize used #leaves
    const double f = 1. / (nrArcs + 1);
    for (NodeIt v_j(_G); v_j != lemon::INVALID; ++v_j)
    {
      if (OutArcIt(_G, v_j) == lemon::INVALID)
      {
        Node v_i = _G.source(InArcIt(_G, v_j));
        assert(v_i != lemon::INVALID);
        
        obj += f * (1 - _x[(*_pNodeToIndex)[v_j]][nrSamples]);
        
        int t = _sampleToIndex[_l[v_j]];
        int mapped_i = _F.characterToIndex(_T.label(_GtoT[v_i]));
        
        obj += f * f * (_Umax[t][mapped_i] - _u[t][mapped_i]);
      }
    }
  }
  
  if (upperBound >= 0)
  {
    _model.addConstr(obj <= upperBound);
  }
  
  _model.setObjective(obj, GRB_MINIMIZE);
  _model.update();
}

void IlpSolverExt::processSolution()
{
  const int nrSamples = _indexToSample.size();
  const int nrCharacters = _F.getNrCharacters();
  
  for (NodeIt v_i(_G); v_i != lemon::INVALID; ++v_i)
  {
    int i = (*_pNodeToIndex)[v_i];
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
    int ij = (*_pArcToIndex)[a_ij];
    
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
    .nodeMap(*_pNodeToIndex, resNodeToIndex)
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

void IlpSolverExt::labelEdges(const Digraph& resT,
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

void IlpSolverExt::run(const NonBinaryCloneTree& T,
                       const FrequencyMatrix& F,
                       const std::string& primary,
                       const std::string& outputDirectory,
                       const StringToIntMap& colorMap,
                       MigrationGraph::Pattern pattern,
                       int nrThreads,
                       bool outputILP,
                       bool outputSearchGraph,
                       int timeLimit,
                       double UB,
                       const StringPairList& forcedComigrations)
{
  char buf[1024];
  std::string filenameGurobiLog;
  std::string filenameSearchGraph;
  
  if (!outputDirectory.empty())
  {
    snprintf(buf, 1024, "%s/log-%s-%s-binarized.txt",
             outputDirectory.c_str(),
             primary.c_str(),
             MigrationGraph::getPatternString(pattern).c_str());
    
    filenameGurobiLog = buf;
    
    snprintf(buf, 1024, "%s/searchG-%s-%s-binarized.dot",
             outputDirectory.c_str(),
             primary.c_str(),
             MigrationGraph::getPatternString(pattern).c_str());
    
    filenameSearchGraph = buf;
  }
  
  {
    IlpSolverExt solver(T,
                        F,
                        primary,
                        pattern,
                        filenameGurobiLog,
                        forcedComigrations);
    solver.init(UB);
    
    if (!outputDirectory.empty() && outputILP)
    {
      snprintf(buf, 1024, "%s/ilp-%s-%s-binarized.lp",
               outputDirectory.c_str(),
               primary.c_str(),
               MigrationGraph::getPatternString(pattern).c_str());
      solver.exportModel(buf);
    }
    
    if (outputSearchGraph)
    {
      std::ofstream outSearchG(filenameSearchGraph.c_str());
      solver.writeDOT(outSearchG, colorMap);
      outSearchG.close();
    }
    
    lemon::Timer timer;
    std::cerr << "With primary '" << primary << "', "
      << MigrationGraph::getPatternLongString(pattern)
      << " and binarization: " << std::flush;
    if (!solver.solve(nrThreads, timeLimit))
    {
      std::cerr << "No solution found (" << outputDirectory << ")" << std::endl;
//        continue;
      return;
    }
    
    MigrationGraph G(solver.T(), solver.lPlus());
    
    std::cerr << G.getNrMigrations() << " migrations, "
      << G.getNrComigrations(solver.T(), solver.lPlus()) << " comigrations, "
      << G.getNrNonUniqueParentageSamples() << " non-unique parentage sites and "
      << G.getNrSeedingSamples() << " seeding sites";
    if (G.hasReseeding())
    {
      std::cerr << " including reseeding";
    }
    std::cerr << ". [LB, UB] = [" << solver.LB() << ", " << solver.UB() << "]. " << timer.realTime() << " seconds (" << outputDirectory << ")" << std::endl;
    
    if (!outputDirectory.empty())
    {
      snprintf(buf, 1024, "%s/T-%s-%s-binarized.dot",
               outputDirectory.c_str(),
               primary.c_str(),
               MigrationGraph::getPatternString(pattern).c_str());
      std::ofstream outT(buf);
      solver.T().writeDOT(outT,
                          solver.lPlus(),
                          colorMap,
                          solver.getF(),
                          solver.getU());
      outT.close();
      
      snprintf(buf, 1024, "%s/T-%s-%s-binarized-condensed.dot",
               outputDirectory.c_str(),
               primary.c_str(),
               MigrationGraph::getPatternString(pattern).c_str());
      std::ofstream outCondensedT(buf);
      solver.T().writeDOT(outCondensedT,
                          solver.lPlus(),
                          colorMap,
                          solver.getU(),
                          solver.getCharacterLabel());
      outCondensedT.close();
      
      snprintf(buf, 1024, "%s/G-%s-%s-binarized.dot",
               outputDirectory.c_str(),
               primary.c_str(),
               MigrationGraph::getPatternString(pattern).c_str());
      
      std::ofstream outG(buf);
      G.writeDOT(outG, colorMap);
      outG.close();
    }
  }
}
