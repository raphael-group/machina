/*
 * ilpsolverext.cpp
 *
 *  Created on: 13-apr-2017
 *      Author: M. El-Kebir
 */

#include "ilpsolverext.h"

IlpSolverExt::IlpSolverExt(const NonBinaryCloneTree& T,
                           const FrequencyMatrix& F,
                           const std::string& primary,
                           Mode mode,
                           const std::string& gurobiLogFilename,
                           const StringPairList& forcedComigrations)
  : _T(T)
  , _F(F)
  , _Fmin(DoubleMatrix(F.getNrSamples(), DoubleVector(F.getNrCharacters(), 0)))
  , _Fmax(DoubleMatrix(F.getNrSamples(), DoubleVector(F.getNrCharacters(), 0)))
  , _Umax(DoubleMatrix(F.getNrSamples(), DoubleVector(F.getNrCharacters(), 0)))
  , _primary(primary)
  , _mode(mode)
  , _forcedComigrations(forcedComigrations)
  , _G()
  , _label(_G)
  , _l(_G)
  , _lPlus(_G)
  , _TtoG(_T.tree(), lemon::INVALID)
  , _GtoT(_G, lemon::INVALID)
  , _indexToArc()
  , _arcToIndex(_G)
  , _indexToNode()
  , _nodeToIndex(_G)
  , _sampleToIndex()
  , _indexToSample()
  , _primaryIndex(-1)
  , _env(gurobiLogFilename)
  , _model(_env)
  , _x()
  , _y()
  , _z()
  , _c()
  , _d()
  , _f()
  , _LB(-1)
  , _UB(-1)
  , _resF(_T.tree(), DoubleVector(F.getNrSamples(), 0))
  , _resU(_T.tree(), DoubleVector(F.getNrSamples(), 0))
  , _pResCloneTree(NULL)
  , _pResLPlus(NULL)
  , _pResF(NULL)
  , _pResU(NULL)
  , _pResCharacterLabel(NULL)
{
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

  initIndices();
  initVariables();
  initConstraints();
  initLeafConstraints();
  initForcedComigrations();
  initObjective(upperBound);
  _model.update();
}

void IlpSolverExt::constructG()
{
  const int m = _F.getNrSamples();
  const int n = _F.getNrCharacters();
  
  const Digraph& T = _T.tree();
  
  lemon::digraphCopy(T, _G)
    .nodeMap(_T.getLabelMap(), _label)
    .nodeRef(_TtoG)
    .nodeCrossRef(_GtoT)
    .node(_T.root(), _root)
    .run();
  
  _indexToNode = NodeVector(_F.getNrCharacters(), lemon::INVALID);
  for (NodeIt v_i(_G); v_i != lemon::INVALID; ++v_i)
  {
    int i = _F.characterToIndex(_T.label(_GtoT[v_i]));
    _indexToNode[i] = v_i;
    _nodeToIndex[v_i] = i;
  }
  
  // attach sample s to vertex i if Umax[s][i] > 0
  for (int s = 0; s < m; ++s)
  {
    const std::string& sStr = _F.indexToSample(s);
    for (int i = 0; i < n; ++i)
    {
      Node v_i = _indexToNode[i];
      Node v_i_inT = _GtoT[v_i];
      
      if ((_T.isLeaf(v_i_inT) && _Fmin[s][i] > 0) ||
          (!_T.isLeaf(v_i_inT) && _Umax[s][i] > 0))
      {
        Node v_i = _indexToNode[i];
        Node v_is = _G.addNode();
        _label[v_is] = _label[v_i] + "_" + sStr;
        _l[v_is] = sStr;
        _nodeToIndex[v_is] = _indexToNode.size();
        _indexToNode.push_back(v_is);
        _G.addArc(v_i, v_is);
      }
    }
  }
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
    const int i = _nodeToIndex[v_i];
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
    const int i = _nodeToIndex[v_i];
    const int j = _nodeToIndex[v_j];
    
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
    const int i = _nodeToIndex[v_i];
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
    const int i = _nodeToIndex[v_i];
    const int j = _nodeToIndex[v_j];
    
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
    _arcToIndex[a] = _indexToArc.size();
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
    
    int i = _nodeToIndex[v_i];
    assert(_sampleToIndex.count(l(v_i)) == 1);
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
  const int nrNodes = _indexToNode.size();
  const int nrCharacters = _F.getNrCharacters();
  const int nrSamples = _indexToSample.size();
  
  GRBLinExpr sum;
  
  // Unique color
  for (NodeIt v_i(_G); v_i != lemon::INVALID; ++v_i)
  {
//    if (isLeaf(v_i)) continue;

    int i = _nodeToIndex[v_i];
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
    int i = _nodeToIndex[v_i];
    int j = _nodeToIndex[v_j];
    
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
  _model.addConstr(_x[_nodeToIndex[root()]][_primaryIndex] == 1);

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
          int i = _nodeToIndex[v_i];
          int j = _nodeToIndex[v_j];
          
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
  
  if (_mode != RESEEDING)
  {
    for (OutArcIt a_ij(_G, root()); a_ij != lemon::INVALID; ++a_ij)
    {
      Node v_j = _G.target(a_ij);
      int j = _nodeToIndex[v_j];
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
        assert(_sampleToIndex.count(l(v_j)) == 1);
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
  
  if (_mode == SINGLE_SOURCE_SEEDING || _mode == PARALLEL_SINGLE_SOURCE_SEEDING)
  {
    initSingleSourceSeedingConstraints();
    if (_mode == PARALLEL_SINGLE_SOURCE_SEEDING)
    {
      initParallelSingleSourceSeedingConstraints();
    }
  }
  else if (_mode == MULTI_SOURCE_SEEDING)
  {
    initMultiSourceSeedingConstraints();
  }
}

void IlpSolverExt::initPrimaryConstraint(Node v_j)
{
  GRBLinExpr sum;
  
  for (InArcIt a_ij(_G, v_j); a_ij != lemon::INVALID; ++a_ij)
  {
    Node v_i = _G.source(a_ij);
    int i = _nodeToIndex[v_i];
    sum += _x[i][_primaryIndex];
  }
  
  if (sum.size() > 0)
  {
    _model.addConstr(sum >= 1);
    if (sum.size() == 1)
    {
      Node v_i = _G.source(InArcIt(_G, v_j));
      initPrimaryConstraint(v_i);
    }
  }
}

void IlpSolverExt::initMultiSourceSeedingConstraints()
{
  const int nrSamples = _indexToSample.size();
  assert(0 <= _primaryIndex && _primaryIndex < nrSamples);
  
  IntVector mets;
  for (int s = 0; s < nrSamples; ++s)
  {
    if (s != _primaryIndex)
    {
      mets.push_back(s);
    }
  }
  
  GRBLinExpr sum;
  
  assert(mets.size() == nrSamples - 1);
  BoolVector subset(mets.size(), false);
  while (next(subset))
  {
    IntVector selectedMets;
    for (int i = 0; i < nrSamples - 1; ++i)
    {
      if (subset[i])
        selectedMets.push_back(mets[i]);
    }
    
    if (selectedMets.size() < 2)
      continue;
    
    // make all permutations
    do
    {
      // add cycle inequality
      for (int i = 1; i < selectedMets.size(); ++i)
      {
        sum += _c[selectedMets[i-1]][selectedMets[i]];
      }
      _model.addConstr(sum <= selectedMets.size() - 1);
      sum.clear();
    } while(std::next_permutation(selectedMets.begin(), selectedMets.end()));
  }
  
  // primary has no incoming edges
  for (int t = 0; t < nrSamples; ++t)
  {
    sum += _c[t][_primaryIndex];
  }
  _model.addConstr(sum == 0);
  sum.clear();
}

bool IlpSolverExt::next(BoolVector& subset)
{
  int n = subset.size();
  for (int i = 0; i < n; ++i)
  {
    if (subset[i] == false)
    {
      subset[i] = true;
      for (int j = 0; j < i; ++j)
      {
        subset[j] = false;
      }
      return true;
    }
  }
  
  return false;
}

void IlpSolverExt::initSingleSourceSeedingConstraints()
{
  const int nrSamples = _indexToSample.size();
  
  GRBLinExpr sum;
  
  // Comigration constraints
  for (int s = 0; s < nrSamples; ++s)
  {
    if (s == _primaryIndex)
    {
      for (int t = 0; t < nrSamples; ++t)
      {
        sum += _c[t][s];
      }
      _model.addConstr(sum == 0);
      sum.clear();
    }
    else
    {
      for (int t = 0; t < nrSamples; ++t)
      {
        sum += _c[t][s];
      }
      _model.addConstr(sum == 1);
      sum.clear();
    }
  }
  
  for (int s = 0; s < nrSamples; ++s)
  {
//    _model.addConstr(_e[s] == 0);
    for (int t = 0; t < nrSamples; ++t)
    {
      sum += _c[s][t];
    }
  }
  
  _model.addConstr(sum == nrSamples - 1);
  sum.clear();
}

void IlpSolverExt::initParallelSingleSourceSeedingConstraints()
{
  const int nrSamples = _indexToSample.size();
  
  for (int s = 0; s < nrSamples; ++s)
  {
    if (s != _primaryIndex)
    {
      _model.addConstr(_c[_primaryIndex][s] == 1);
    }
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
        
        obj += f * (1 - _x[_nodeToIndex[v_j]][nrSamples]);
        
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

void IlpSolverExt::initForcedComigrations()
{
  for (const StringPair& st : _forcedComigrations)
  {
    assert(_sampleToIndex.count(st.first) == 1);
    assert(_sampleToIndex.count(st.second) == 1);
    
    int s = _sampleToIndex[st.first];
    int t = _sampleToIndex[st.second];
  
    _model.addConstr(_c[s][t] == 1);
  }
}

void IlpSolverExt::processSolution()
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

    if (isLeaf(v_j))
    {
      int s = _sampleToIndex[_l[v_j]];
      if (_resU[_GtoT[v_i]][s] == 0)
      {
        filterNodeMap[v_j] = false;
        filterArcMap[a_ij] = false;
      }
    }
  }
  
  lemon::SubDigraph<Digraph> TT(_G, filterNodeMap, filterArcMap);
  Digraph resT;
  StringNodeMap resLPlus(resT);
  StringNodeMap resLabel(resT);
  DoubleVectorNodeMap resF(resT);
  Node resRoot;
  lemon::digraphCopy(TT, resT)
    .nodeMap(_resF, resF)
    .nodeMap(_lPlus, resLPlus)
    .nodeMap(_label, resLabel)
    .node(root(), resRoot)
    .run();
  
  _pResCloneTree = new NonBinaryCloneTree(resT, resRoot, resLabel, resLPlus);
  
  _pResLPlus = new StringNodeMap(_pResCloneTree->tree());
  _pResF = new DoubleVectorNodeMap(_pResCloneTree->tree());
  for (NodeIt v(resT); v != lemon::INVALID; ++v)
  {
    Node new_v = _pResCloneTree->getNodeByLabel(resLabel[v]);
    _pResLPlus->set(new_v, resLPlus[v]);
    _pResF->set(new_v, resF[v]);
  }
  
  assert(_lPlus[root()] == _primary);
}

void IlpSolverExt::exportModel(const std::string& filename)
{
  _model.write(filename);
}

bool IlpSolverExt::solve(int nrThreads, int timeLimit)
{
  try
  {
    if (nrThreads > 0)
    {
      _model.getEnv().set(GRB_IntParam_Threads, nrThreads);
    }
    if (timeLimit > 0)
    {
      _model.getEnv().set(GRB_DoubleParam_TimeLimit, timeLimit);
    }
    _model.getEnv().set(GRB_IntParam_LogToConsole, 0);
    _model.optimize();
    int status = _model.get(GRB_IntAttr_Status);
    if (status == GRB_OPTIMAL || status == GRB_TIME_LIMIT)
    {
      _LB = _model.get(GRB_DoubleAttr_ObjBound);
      _UB = _model.get(GRB_DoubleAttr_ObjVal);
      processSolution();
      return true;
    }
    else if (status == GRB_INF_OR_UNBD)
    {
//      std::cerr << "Model is infeasible or unbounded" << std::endl;
      return false;
    }
    else if (status == GRB_INFEASIBLE)
    {
//      std::cerr << "Model is infeasible" << std::endl;
//      _model.computeIIS();
//      _model.write("/tmp/model_IIS.ilp");
      return false;
    }
    else if (status == GRB_UNBOUNDED)
    {
      std::cerr << "Model is unbounded" << std::endl;
      return false;
    }
  }
  catch (const GRBException& e)
  {
    std::cerr << "Error code = " << e.getErrorCode() << std::endl;
    std::cerr << e.getMessage() << std::endl;
    
    return false;
  }
  catch (...)
  {
    return false;
  }
  
  return true;
}


