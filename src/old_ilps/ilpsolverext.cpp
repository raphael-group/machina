/*
 * ilpsolverext.cpp
 *
 *  Created on: 13-apr-2017
 *      Author: M. El-Kebir
 */

#include "ilpsolverext.h"
#include <lemon/time_measure.h>

IlpSolverExt::IlpSolverExt(const CloneTree& T,
                           const FrequencyMatrix& F,
                           const std::string& primary,
                           MigrationGraph::Pattern pattern,
                           const std::string& gurobiLogFilename,
                           const StringPairList& forcedComigrations)
  : IlpSolver(T, primary, pattern, gurobiLogFilename, forcedComigrations)
  , _F(F)
  , _Fmin(DoubleMatrix(F.getNrSamples(),
                       DoubleVector(F.getNrCharacters(), 0)))
  , _Fmax(DoubleMatrix(F.getNrSamples(),
                       DoubleVector(F.getNrCharacters(), 0)))
  , _Umax(DoubleMatrix(F.getNrSamples(),
                       DoubleVector(F.getNrCharacters(), 0)))
  , _sampleToIndex()
  , _indexToSample()
  , _sampleIndexToAnatomicalSiteIndex()
  , _anatomicalSiteIndexToSampleIndices()
  , _G()
  , _label(_G)
  , _l(_G)
  , _lPlus(_G)
  , _TtoG(_T.tree(), lemon::INVALID)
  , _GtoT(_G, lemon::INVALID)
  , _isAnatomicalSiteNode(_G, false)
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

void IlpSolverExt::init(const IntTriple& bounds)
{
  computeFmin(_T.root());
  computeFmax(_T.root());
  computeUmax();
  
  IlpSolver::init(bounds);
}

void IlpSolverExt::constructG()
{
  constructG(_T.root());
}

void IlpSolverExt::computeFmin(Node v_i)
{
  const int k = _F.getNrSamples();
  const int i = _F.characterToIndex(_T.label(v_i));
  
  const Digraph& T = _T.tree();
  
  if (OutArcIt(T, v_i) == lemon::INVALID)
  {
    for (int p = 0; p < k; ++p)
    {
      _Fmin[p][i] = _F.min(p,i);
    }
  }
  else
  {
    for (OutArcIt a_ij(T, v_i); a_ij != lemon::INVALID; ++a_ij)
    {
      Node v_j = T.target(a_ij);
      computeFmin(v_j);
    }

    for (int p = 0; p < k; ++p)
    {
      double sum_p = 0;
      for (OutArcIt a_ij(T, v_i); a_ij != lemon::INVALID; ++a_ij)
      {
        Node v_j = T.target(a_ij);
        int j = _F.characterToIndex(_T.label(v_j));
        sum_p += _Fmin[p][j];
      }
      _Fmin[p][i] = std::max(_F.min(p, i), sum_p);
    }
  }
}

void IlpSolverExt::constructG(Node v)
{
  const int m = _F.getNrAnatomicalSites();
  const int k = _F.getNrSamples();
  
  const Digraph& T = _T.tree();
  
  char buf[1024];
  const int mapped_i = _F.characterToIndex(_T.label(v));

  int out_deg = lemon::countOutArcs(T, v);
  
  NodeVector anatomicalSiteNodes;
  // attach anatomical site s to vertex i if there exists a p in sigma(s) such that Umax[p][i] > 0
  for (int s = 0; s < m; ++s)
  {
    const std::string& sStr = _F.indexToAnatomicalSite(s);
    for (int p : _anatomicalSiteIndexToSampleIndices[s])
    {
      assert(0 <= p && p < k);
      if ((_T.isLeaf(v) && _Fmin[p][mapped_i] > 0) ||
          (!_T.isLeaf(v) && _Umax[p][mapped_i] > 0))
      {
        Node v_is = _G.addNode();
        _isAnatomicalSiteNode[v_is] = true;
        _GtoT[v_is] = lemon::INVALID;
        _label[v_is] = _T.label(v) + "_" + sStr;
        _l[v_is] = sStr;
        _pNodeToIndex->set(v_is, _indexToNode.size());
        _indexToNode.push_back(v_is);
        anatomicalSiteNodes.push_back(v_is);
        break;
      }
    }
  }
  
  out_deg += anatomicalSiteNodes.size();
  if (out_deg > 2)
  {
    NodeVector nodes;
    for (int i = 0; i < out_deg - 1; ++i)
    {
      Node vv = _G.addNode();
      _isAnatomicalSiteNode[vv] = false;
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
    
    // attach anatomical sites
    int is = 0;
    for (Node w : anatomicalSiteNodes)
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
    _isAnatomicalSiteNode[vv] = false;
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
    
    // attach anatomical sites
    for (Node w : anatomicalSiteNodes)
    {
      _G.addArc(vv, w);
    }
  }
}

void IlpSolverExt::computeFmax(Node v_i)
{
  const int k = _F.getNrSamples();
  const int i = _F.characterToIndex(_T.label(v_i));
  const Digraph& T = _T.tree();

  if (v_i == _T.root())
  {
    for (int p = 0; p < k; ++p)
    {
      _Fmax[p][i] = _F.max(p, i);
    }
  }
  else
  {
    Node v_pi_i = T.source(InArcIt(T, v_i));
    int pi_i = _F.characterToIndex(_T.label(v_pi_i));
    
    for (int p = 0; p < k; ++p)
    {
      double sum_p = _Fmax[p][pi_i];
      for (OutArcIt a_pi_i_j(T, v_pi_i); a_pi_i_j != lemon::INVALID; ++a_pi_i_j)
      {
        Node v_j = T.target(a_pi_i_j);
        if (v_j == v_i) continue;

        int j = _F.characterToIndex(_T.label(v_j));
        sum_p -= _Fmin[p][j];
      }
      
      _Fmax[p][i] = std::min(_F.max(p, i), sum_p);
      assert(_Fmax[p][i] >= 0);
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
  const int k = _F.getNrSamples();
  const int n = _F.getNrCharacters();
  const Digraph& T = _T.tree();
  
  // compute Umax
  for (int p = 0; p < k; ++p)
  {
    for (int i = 0; i < n; ++i)
    {
      Node v_i = _T.getNodeByLabel(_F.indexToCharacter(i));
      if (v_i == lemon::INVALID || (_T.isLeaf(v_i) && _Fmin[p][i] == 0))
      {
        _Umax[p][i] = 0;
      }
      else
      {
        _Umax[p][i] = _Fmax[p][i];
        for (OutArcIt a_ij(T, v_i); a_ij != lemon::INVALID; ++a_ij)
        {
          Node v_j = T.target(a_ij);
          int j = _F.characterToIndex(_T.label(v_j));
          
          _Umax[p][i] -= _Fmin[p][j];
        }
      }
    }
  }
}

void IlpSolverExt::writeDOT(std::ostream& out) const
{
  const int k = _F.getNrSamples();
  
  out << "digraph T {" << std::endl;
  
  for (NodeIt v_i(_G); v_i != lemon::INVALID; ++v_i)
  {
    const int i = (*_pNodeToIndex)[v_i];
    out << "\t" << i << " [label=\"" << _label[v_i];
    if (!isLeaf(v_i))
    {
      int mapped_i = _F.characterToIndex(_T.label(_GtoT[v_i]));
      for (int p = 0; p < k; ++p)
      {
        out << "\n" << _F.indexToSample(p)
            << " : [" << _F.min(p, mapped_i) << "," << _F.max(p, mapped_i) << "]"
            << " ; [" << _Fmin[p][mapped_i] << "," << _Fmax[p][mapped_i] << "]"
            << " ; " << _Umax[p][mapped_i];
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
  const int k = _F.getNrSamples();
  
  out << "digraph T {" << std::endl;
  
  for (NodeIt v_i(_G); v_i != lemon::INVALID; ++v_i)
  {
    const int i = (*_pNodeToIndex)[v_i];
    out << "\t" << i << " [label=\"" << _label[v_i];
    if (!isLeaf(v_i))
    {
      int mapped_i = _F.characterToIndex(_T.label(_GtoT[v_i]));
      for (int p = 0; p < k; ++p)
      {
        out << "\n" << _F.indexToSample(p)
            << " : [" << _F.min(p, mapped_i) << "," << _F.max(p, mapped_i) << "]"
            << " ; [" << _Fmin[p][mapped_i] << "," << _Fmax[p][mapped_i] << "]"
            << " ; " << _Umax[p][mapped_i] << " ; " << _resF[_GtoT[v_i]][p] << " ; " << _resU[_GtoT[v_i]][p];
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
  const int m = _F.getNrAnatomicalSites();
  const int k = _F.getNrSamples();
  
  _sampleToIndex.clear();
  _indexToSample = StringVector(k);
  _sampleIndexToAnatomicalSiteIndex = IntVector(k, -1);
  _anatomicalSiteIndexToSampleIndices = IntSetVector(m);
  _sampleIndexToAnatomicalSiteIndex.clear();
  for (int p = 0; p < k; ++p)
  {
    const std::string& pStr = _F.indexToSample(p);
    const int s = _F.sampleIndexToAnatomicalSiteIndex(p);
    
    _sampleToIndex[pStr] = p;
    _indexToSample[p] = pStr;
    
    _sampleIndexToAnatomicalSiteIndex[p] = s;
    _anatomicalSiteIndexToSampleIndices[s].insert(p);
  }
  
  _indexToAnatomicalSite = StringVector(m);
  
  // construct G
  constructG();
  
  _indexToArc.clear();
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    _pArcToIndex->set(a, _indexToArc.size());
    _indexToArc.push_back(a);
  }

  _primaryIndex = -1;
  _anatomicalSiteToIndex.clear();
  for (int s = 0; s < m; ++s)
  {
    _anatomicalSiteToIndex[_F.indexToAnatomicalSite(s)] = s;
    _indexToAnatomicalSite[s] = _F.indexToAnatomicalSite(s);
    if (_primary == _F.indexToAnatomicalSite(s))
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
  const int nrAnatomicalSites = _indexToAnatomicalSite.size();
  const int nrSamples = _indexToSample.size();
  const int nrCharacters = _F.getNrCharacters();
  
  /// x[i][s] = 1 iff vertex v_i is labeled by anatomical site s
  char buf[1024];
  _x = VarMatrix(nrNodes);
  for (int i = 0; i < nrNodes; ++i)
  {
    _x[i] = VarArray(nrAnatomicalSites + 1);
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      snprintf(buf, 255, "x_%s_%s",
               label(_indexToNode[i]).c_str(),
               _indexToAnatomicalSite[s].c_str());
      _x[i][s] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
    }
    snprintf(buf, 255, "x_%s_DUMMY", label(_indexToNode[i]).c_str());
    _x[i][nrAnatomicalSites] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
  }
  
  /// y[ij] = 1 iff edge (v_i, v_j) is a migration edge
  _y = VarArray(nrArcs);
  for (int ij = 0; ij < nrArcs; ++ij)
  {
    Arc a_ij = _indexToArc[ij];
    Node v_i = _G.source(a_ij);
    Node v_j = _G.target(a_ij);
    
    snprintf(buf, 255, "y_%s_%s",
             label(v_i).c_str(),
             label(v_j).c_str());
    _y[ij] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
  }
  
  /// z[ij][s] = 1 iff v_i and v_j are labeled by anatomical site s
  _z = VarMatrix(nrArcs);
  for (int ij = 0; ij < nrArcs; ++ij)
  {
    Arc a_ij = _indexToArc[ij];
    Node v_i = _G.source(a_ij);
    Node v_j = _G.target(a_ij);
    
    _z[ij] = VarArray(nrAnatomicalSites+1);
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      snprintf(buf, 255, "z_%s_%s_%s",
               label(v_i).c_str(),
               label(v_j).c_str(),
               _indexToAnatomicalSite[s].c_str());
      
      _z[ij][s] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
    }
    snprintf(buf, 255, "z_%s_%s_DUMMY",
             label(v_i).c_str(),
             label(v_j).c_str());
    _z[ij][nrAnatomicalSites] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
  }
  
  /// c[s][t] = 1 iff there exists a migration edge (v_i, v_j) where l(v_i) = s and l(v_j) = t
  _c = VarMatrix(nrAnatomicalSites);
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    _c[s] = VarArray(nrAnatomicalSites);
    for (int t = 0; t < nrAnatomicalSites; ++t)
    {
      snprintf(buf, 255, "c_%s_%s",
               _indexToAnatomicalSite[s].c_str(),
               _indexToAnatomicalSite[t].c_str());
      _c[s][t] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
    }
  }
  
  /// d[s] = 1 iff there exists a migration edge (v_i, v_j) where l(v_i) = s
  _d = VarArray(nrAnatomicalSites);
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    snprintf(buf, 255, "d_%s",
             _indexToAnatomicalSite[s].c_str());
    _d[s] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
  }
  
  /// f[p][i] = frequency of mutation incoming to clone i in sample p
  _f = VarMatrix(nrSamples);
  for (int p = 0; p < nrSamples; ++p)
  {
    _f[p] = VarArray(nrCharacters);
    for (int c = 0; c < nrCharacters; ++c)
    {
      snprintf(buf, 255, "f_%s_%s",
               _indexToSample[p].c_str(),
               _F.indexToCharacter(c).c_str());
      _f[p][c] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
    }
  }

  /// u[p][i] = usage of clone i in sample p
  _u = VarMatrix(nrSamples);
  for (int p = 0; p < nrSamples; ++p)
  {
    _u[p] = VarArray(nrCharacters);
    for (int c = 0; c < nrCharacters; ++c)
    {
      snprintf(buf, 255, "u_%s_%s",
               _indexToSample[p].c_str(),
               _F.indexToCharacter(c).c_str());
      _u[p][c] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
    }
  }
  
  // w[ij] = 1 iff edge (v_i, v_j) is in the tree
  _w = VarArray(nrArcs);
  for (int ij = 0; ij < nrArcs; ++ij)
  {
    Arc a_ij = _indexToArc[ij];
    Node v_i = _G.source(a_ij);
    Node v_j = _G.target(a_ij);
    
    snprintf(buf, 255, "w_%s_%s",
             label(v_i).c_str(),
             label(v_j).c_str());
    _w[ij] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
  }

  _model.update();
}

void IlpSolverExt::initLeafConstraints()
{
  const int nrAnatomicalSites = _indexToAnatomicalSite.size();
  
  GRBLinExpr sum;
  
  // Leaf color
  for (NodeIt v_i(_G); v_i != lemon::INVALID; ++v_i)
  {
    if (!isLeaf(v_i))
    {
      continue;
    }
    
    // each leaf is an anatomical site node v_is
    int i = (*_pNodeToIndex)[v_i];
    assert(_anatomicalSiteToIndex.count(_l[v_i]) == 1);
    int s = _anatomicalSiteToIndex[l(v_i)];
    
    _model.addConstr(_x[i][s] + _x[i][nrAnatomicalSites] == 1);
    for (int t = 0; t < nrAnatomicalSites; ++t)
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
  const int nrAnatomicalSites = _indexToAnatomicalSite.size();
  const int nrSamples = _indexToSample.size();
  
  GRBLinExpr sum;
  
  // Unique color
  for (NodeIt v_i(_G); v_i != lemon::INVALID; ++v_i)
  {
    int i = (*_pNodeToIndex)[v_i];
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      sum += _x[i][s];
    }
    sum += _x[i][nrAnatomicalSites];
    _model.addConstr(sum == 1);
    if (!isLeaf(v_i))
    {
      _model.addConstr(_x[i][nrAnatomicalSites] == 0);
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
    
    for (int s = 0; s < nrAnatomicalSites + 1; ++s)
    {
      sum += _z[ij][s];
      
      if (s < nrAnatomicalSites)
      {
        _model.addConstr(_z[ij][s] <= _x[i][s]);
        _model.addConstr(_z[ij][s] <= _x[j][s]);
        //_model.addConstr(_z[ij][s] >= _x[i][s] + _x[j][s] - 1);
      }
      else if (isLeaf(v_j))
      {
        int mapped_i = _F.characterToIndex(_T.label(_GtoT[v_i]));
        assert(0 <= mapped_i && mapped_i < nrCharacters);
        
        int t = _anatomicalSiteToIndex[_l[v_j]];
        assert(0 <= t && t < nrAnatomicalSites);
        assert(!_anatomicalSiteIndexToSampleIndices[t].empty());
        for (int p : _anatomicalSiteIndexToSampleIndices[t])
        {
          assert(0 <= p && p < nrSamples);
          addMatchingColorsAndUsageConstraint(ij, mapped_i, p);
          _model.addConstr(_x[j][t] >= _u[p][mapped_i]);
        }
      }
      else
      {
        _model.addConstr(_z[ij][nrAnatomicalSites] == 0);
      }
    }
    addMatchingColorsConstraint(sum, ij);
    
    sum.clear();
  }
  
  // Root color
  _model.addConstr(_x[(*_pNodeToIndex)[root()]][_primaryIndex] == 1);

  // Comigration constraints
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    for (int t = 0; t < nrAnatomicalSites; ++t)
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
  
  if (nrAnatomicalSites > 1)
  {
    _model.addConstr(_d[_primaryIndex] == 1);
  }
  
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    for (int t = 0; t < nrAnatomicalSites; ++t)
    {
      _model.addConstr(_d[s] >= _c[s][t]);
    }
  }
  
  for (NodeIt v_i(_T.tree()); v_i != lemon::INVALID; ++v_i)
  {
    int mapped_i = _F.characterToIndex(_T.label(v_i));
    for (int p = 0; p < nrSamples; ++p)
    {
      for (OutArcIt a_ij(_T.tree(), v_i); a_ij != lemon::INVALID; ++a_ij)
      {
        Node v_j = _T.tree().target(a_ij);
        int mapped_j = _F.characterToIndex(_T.label(v_j));
        
        sum += _f[p][mapped_j];
      }
      _model.addConstr(_f[p][mapped_i] >= _Fmin[p][mapped_i]);
      _model.addConstr(_f[p][mapped_i] <= _Fmax[p][mapped_i]);
      _model.addConstr(_u[p][mapped_i] <= _Umax[p][mapped_i]);
      _model.addConstr(_f[p][mapped_i] >= sum);
      _model.addConstr(_u[p][mapped_i] == _f[p][mapped_i] - sum);
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
        assert(_anatomicalSiteToIndex.count(_l[v_j]) == 1);
        int s = _anatomicalSiteToIndex[l(v_j)];
        if (s == _primaryIndex
            && OutArcIt(_T.tree(), _GtoT[v_i]) == lemon::INVALID)
        {
          for (int p : _anatomicalSiteIndexToSampleIndices[s])
          {
            if (_Fmin[p][mapped_i] > 0)
            {
              initPrimaryConstraint(v_j);
              break;
            }
          }
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

void IlpSolverExt::initObjective(const IntTriple& bounds)
{
  const int nrArcs = _indexToArc.size();
  const int nrAnatomicalSites = _indexToAnatomicalSite.size();
  
  GRBLinExpr obj, sum_migrations, sum_comigrations, sum_seeding_sites;
  // migrations
  for (int ij = 0; ij < nrArcs; ++ij)
  {
    obj += _y[ij];
    sum_migrations += _y[ij];
  }
  _model.addConstr(sum_migrations >= nrAnatomicalSites - 1);
  
  if (bounds.first != -1)
  {
    _model.addConstr(sum_migrations <= bounds.first);
  }
  
//  if (_forcedComigrations.empty())
  {
    // comigrations
    {
      const double f = 1. / (nrArcs + 1);
      for (int s = 0; s < nrAnatomicalSites; ++s)
      {
        for (int t = 0; t < nrAnatomicalSites; ++t)
        {
          obj += f * _c[s][t];
          sum_comigrations += _c[s][t];
        }
      }
    }
    _model.addConstr(sum_comigrations >= nrAnatomicalSites - 1);
    _model.addConstr(sum_migrations >= sum_comigrations);
    
    if (bounds.second.first != -1)
    {
      _model.addConstr(sum_comigrations <= bounds.second.first);
    }
    
    // seeding sites
    {
      const double g = (1. / (nrArcs + 1)) * (1. / (nrAnatomicalSites + 1));
      for (int s = 0; s < nrAnatomicalSites; ++s)
      {
        obj += g * _d[s];
        sum_seeding_sites += _d[s];
      }
    }
    
    if (bounds.second.second != -1)
    {
      _model.addConstr(sum_seeding_sites <= bounds.second.second);
    }
  }
//  else
  {
    // minimize used #leaves
     const double g = (1. / (nrArcs + 1)) * (1. / (nrAnatomicalSites + 1)) * 1. / (nrArcs + 1);
    for (NodeIt v_j(_G); v_j != lemon::INVALID; ++v_j)
    {
      if (OutArcIt(_G, v_j) == lemon::INVALID)
      {
        Node v_i = _G.source(InArcIt(_G, v_j));
        assert(v_i != lemon::INVALID);
        
        obj += g * (1 - _x[(*_pNodeToIndex)[v_j]][nrAnatomicalSites]);
        
        int t = _anatomicalSiteToIndex[_l[v_j]];
        int mapped_i = _F.characterToIndex(_T.label(_GtoT[v_i]));
        
        for (int p : _anatomicalSiteIndexToSampleIndices[t])
        {
          obj += g * g * (_Umax[p][mapped_i] - _u[p][mapped_i]);
        }
      }
    }
  }
  
  _model.setObjective(obj, GRB_MINIMIZE);
  _model.update();
}

void IlpSolverExt::processSolution()
{
  const int nrAnatomicalSites = _indexToAnatomicalSite.size();
  const int nrSamples = _indexToSample.size();
  const int nrCharacters = _F.getNrCharacters();
  
  for (NodeIt v_i(_G); v_i != lemon::INVALID; ++v_i)
  {
    int i = (*_pNodeToIndex)[v_i];
    bool found = false;
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      if (_x[i][s].get(GRB_DoubleAttr_X) >= 0.4)
      {
        _lPlus[v_i] = _indexToAnatomicalSite[s];
        found = true;
        break;
      }
    }
  }
  
  for (NodeIt v_i(_T.tree()); v_i != lemon::INVALID; ++v_i)
  {
    int mapped_i = _F.characterToIndex(_T.label(v_i));
    assert(0 <= mapped_i && mapped_i < nrCharacters);

    for (int p = 0; p < nrSamples; ++p)
    {
      // set U and F
      _resU[v_i][p] = _u[p][mapped_i].get(GRB_DoubleAttr_X);
      _resF[v_i][p] = _f[p][mapped_i].get(GRB_DoubleAttr_X);
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
      int s = _anatomicalSiteToIndex[_l[v_j]];
      bool remove = true;
      for (int p : _anatomicalSiteIndexToSampleIndices[s])
      {
        if (g_tol.nonZero(_resU[_GtoT[v_i]][p]))
        {
          remove = false;
          break;
        }
      }
      
      if (remove)
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
  
  // PRE: unused anatomical site nodes have been removed
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
        // out deg 0 and not a (used) anatomical site node
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
  StringNodeMap resCharacterLabel(resT);
  BoolVector resVisited(_F.getNrCharacters(), false);
  labelEdges(resT,
             resNodeToIndex,
             resRoot,
             resCharacterLabel,
             resVisited);
  
  _pResCloneTree = new CloneTree(resT, resRoot, resLabel, resLPlus);
  _pResLPlus = new StringNodeMap(_pResCloneTree->tree());
  _pResF = new DoubleVectorNodeMap(_pResCloneTree->tree());
  _pResU = new DoubleVectorNodeMap(_pResCloneTree->tree());
  _pResCharacterLabel = new StringNodeMap(_pResCloneTree->tree());
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
    
    if (_isAnatomicalSiteNode[org_v])
    {
      Node parent_org_v = _G.source(InArcIt(_G, org_v));
//      int i = _nodeToIndex[parent_org_v];
      int s = _anatomicalSiteToIndex[_l[org_v]];
      
      (*_pResU)[new_v].clear();
      for (int p : _anatomicalSiteIndexToSampleIndices[s])
      {
        (*_pResU)[new_v].push_back(_resU[_GtoT[parent_org_v]][p]);
      }
    }
  }
  
  assert(_lPlus[root()] == _primary);
}

void IlpSolverExt::labelEdges(const Digraph& resT,
                              const IntNodeMap& resNodeToIndex,
                              Node v_i,
                              StringNodeMap& characterLabel,
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
    characterLabel[v_i] = _F.indexToCharacter(characterIndex);
  }
  
  for (OutArcIt a(resT, v_i); a != lemon::INVALID; ++a)
  {
    Node v_j = resT.target(a);
    labelEdges(resT, resNodeToIndex, v_j, characterLabel, visited);
  }
}

IntTriple IlpSolverExt::run(const CloneTree& T,
                            const FrequencyMatrix& F,
                            const std::string& primary,
                            const std::string& outputDirectory,
                            const std::string& fileNamePrefix,
                            const StringToIntMap& colorMap,
                            MigrationGraph::Pattern pattern,
                            int nrThreads,
                            bool outputILP,
                            bool outputSearchGraph,
                            int timeLimit,
                            const IntTriple& bounds,
                            const StringPairList& forcedComigrations)
{
  char buf[1024];
  std::string filenameGurobiLog;
  std::string filenameSearchGraph;
  
  if (!outputDirectory.empty())
  {
    snprintf(buf, 1024, "%s/%s-log-%s-%s-binarized.txt",
             outputDirectory.c_str(),
             fileNamePrefix.c_str(),
             primary.c_str(),
             MigrationGraph::getPatternString(pattern).c_str());
    
    filenameGurobiLog = buf;
    
    snprintf(buf, 1024, "%s/%s-searchG-%s-%s-binarized.dot",
             outputDirectory.c_str(),
             fileNamePrefix.c_str(),
             primary.c_str(),
             MigrationGraph::getPatternString(pattern).c_str());
    
    filenameSearchGraph = buf;
  }
  
  IlpSolverExt solver(T,
                      F,
                      primary,
                      pattern,
                      filenameGurobiLog,
                      forcedComigrations);
  solver.init(bounds);
  
  if (!outputDirectory.empty() && outputILP)
  {
    snprintf(buf, 1024, "%s/%s-ilp-%s-%s-binarized.lp",
             outputDirectory.c_str(),
             fileNamePrefix.c_str(),
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
  
  bool solved = solver.solve(nrThreads, timeLimit);
  if (!solved)
  {
    std::cout << fileNamePrefix << "\t"
              << "(" << MigrationGraph::getAllowedPatternsString(pattern) << ")\t"
              << "-" << "\t"
              << "-" << "\t"
              << "-" << "\t"
              << "-" << "\t"
              << "-" << "\t"
              << "-" << "\t"
              << timer.realTime()
              << std::endl;
    return bounds;
  }
  
  MigrationGraph G(solver.T(), solver.lPlus());
  
  int mu = G.getNrMigrations();
  int gamma = G.getNrComigrations(solver.T(), solver.lPlus());
  int sigma = G.getNrSeedingSites();
  
  std::cout << fileNamePrefix << "\t"
            << "(" << MigrationGraph::getAllowedPatternsString(pattern) << ")\t"
            << mu << "\t"
            << gamma << "\t"
            << sigma << "\t"
            << G.getPatternString(G.getPattern(), G.isMonoclonal()) << "\t"
            << solver.LB() << "\t"
            << solver.UB() << "\t"
            << timer.realTime()
            << std::endl;

  
  if (!outputDirectory.empty())
  {
    snprintf(buf, 1024, "%s/%s-T-%s-%s-binarized.dot",
             outputDirectory.c_str(),
             fileNamePrefix.c_str(),
             primary.c_str(),
             MigrationGraph::getPatternString(pattern).c_str());
    std::ofstream outT(buf);
    solver.T().writeDOT(outT,
                        solver.lPlus(),
                        colorMap,
                        solver.getU(),
                        solver.getCharacterLabel());
    outT.close();
          
    snprintf(buf, 1024, "%s/%s-G-%s-%s-binarized.dot",
             outputDirectory.c_str(),
             fileNamePrefix.c_str(),
             primary.c_str(),
             MigrationGraph::getPatternString(pattern).c_str());
    std::ofstream outG(buf);
    G.writeDOT(outG, colorMap);
    outG.close();
    
    snprintf(buf, 1024, "%s/%s-T-%s-%s-binarized.tree",
             outputDirectory.c_str(),
             fileNamePrefix.c_str(),
             primary.c_str(),
             MigrationGraph::getPatternString(pattern).c_str());
    std::ofstream outTree(buf);
    solver.T().write(outTree);
    outTree.close();
    
    snprintf(buf, 1024, "%s/%s-T-%s-%s-binarized.labeling",
             outputDirectory.c_str(),
             fileNamePrefix.c_str(),
             primary.c_str(),
             MigrationGraph::getPatternString(pattern).c_str());
    std::ofstream outLabeling(buf);
    solver.T().writeVertexLabeling(outLabeling, solver.lPlus());
    outLabeling.close();
  }

  return std::make_pair(mu, std::make_pair(gamma, sigma));
}
