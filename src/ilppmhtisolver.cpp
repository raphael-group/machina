/*
 * ilppmhtisolver.cpp
 *
 *  Created on: 13-sep-2017
 *      Author: M. El-Kebir
 */

#include "ilppmhtisolver.h"

IlpPmhTiSolver::IlpPmhTiSolver(const CloneTree& T,
                               const FrequencyMatrix& F,
                               const std::string& primary,
                               MigrationGraph::Pattern pattern,
                               const std::string& gurobiLogFilename,
                               const StringPairList& forcedComigrations,
                               bool disablePolytomyResolution)
  : IlpPmhTrSolver(T, primary, pattern, gurobiLogFilename, forcedComigrations)
  , _F(F)
  , _disablePolytomyResolution(disablePolytomyResolution)
  , _Fmin(DoubleMatrix(F.getNrSamples(),
                       DoubleVector(F.getNrCharacters(), 0)))
  , _Fmax(DoubleMatrix(F.getNrSamples(),
                       DoubleVector(F.getNrCharacters(), 0)))
  , _Umax(DoubleMatrix(F.getNrSamples(),
                       DoubleVector(F.getNrCharacters(), 0)))
  , _extT()
  , _rootExtT(lemon::INVALID)
  , _labelExtT( _extT)
  , _leafAnatomicalSiteLabelExtT(_extT)
  , _TtoExtT(IlpPmhTrSolver::getTree(), lemon::INVALID)
  , _extTtoT(_extT, lemon::INVALID)
  , _isAnatomicalSiteNodeExtT(_extT)
  , _f()
  , _u()
  , _pResF(NULL)
  , _pResU(NULL)
{
}

IlpPmhTiSolver::~IlpPmhTiSolver()
{
  delete _pResF;
  delete _pResU;
}

void IlpPmhTiSolver::initLeafVariables()
{
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  
  char buf[1024];
  
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    if (!isLeaf(v_i)) continue;
    
    const int i = (*_pNodeToIndex)[v_i];
    _x[i] = VarMatrix(nrAnatomicalSites + 1);
    for (int s = 0; s < nrAnatomicalSites + 1; ++s)
    {
      int size_L_s = s < nrAnatomicalSites ? _L[s].size() : 1;
      _x[i][s] = VarArray(size_L_s);
      for (int c = 0; c < size_L_s; ++c)
      {
        if (s == nrAnatomicalSites)
        {
          snprintf(buf, 1024, "x;%s;DUMMY;%d",
                   getLabel(v_i).c_str(), c);
        }
        else
        {
          snprintf(buf, 1024, "x;%s;%s;%d",
                   getLabel(v_i).c_str(),
                   _indexToAnatomicalSite[s].c_str(), c);
        }
        _x[i][s][c] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
      }
    }
  }
  
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    if (!isLeaf(v_i)) continue;
    const int i = (*_pNodeToIndex)[v_i];
    _r[i] = VarMatrix(nrAnatomicalSites + 1);
    for (int s = 0; s < nrAnatomicalSites + 1; ++s)
    {
      int size_L_s = s < nrAnatomicalSites ? _L[s].size() : 1;
      _r[i][s] = VarArray(size_L_s);
      for (int c = 0; c < size_L_s; ++c)
      {
        if (s == nrAnatomicalSites)
        {
          snprintf(buf, 1024, "r;%s;DUMMY;%d",
                   getLabel(v_i).c_str(), c);
        }
        else
        {
          snprintf(buf, 1024, "r;%s;%s;%d",
                   getLabel(v_i).c_str(),
                   _indexToAnatomicalSite[s].c_str(), c);
        }
        _r[i][s][c] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
      }
    }
  }
  
  _model.update();
}

void IlpPmhTiSolver::processSolution()
{
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  const int nrSamples = _F.getNrSamples();
  const int nrNodes = _indexToNode.size();
//  const int nrArcs = _indexToArc.size();
  
//  for (int ij = 0; ij < nrArcs; ++ij)
//  {
//    for (int s = 0; s < nrAnatomicalSites; ++s)
//    {
//      int size_L_s = _L[s].size();
//      for (int c = 0; c < size_L_s; ++c)
//      {
//        for (int t = 0; t < nrAnatomicalSites; ++t)
//        {
//          int size_L_t = _L[t].size();
//          for (int d = 0; d < size_L_t; ++d)
//          {
//            if (_xx[ij][s][c][t][d].get(GRB_DoubleAttr_X) >= 0.4)
//            {
//              std::cout << _xx[ij][s][c][t][d].get(GRB_StringAttr_VarName)
//              << " = " << _xx[ij][s][c][t][d].get(GRB_DoubleAttr_X) << std::endl;
//            }
//          }
//        }
//      }
//    }
//  }
//  
//  std::cout << std::endl;
//  
//  for (int i = 0; i < nrNodes; ++i)
//  {
//    for (int s = 0; s < nrAnatomicalSites; ++s)
//    {
//      int size_L_s = _L[s].size();
//      for (int c = 0; c < size_L_s; ++c)
//      {
//        for (int t = 0; t < nrAnatomicalSites; ++t)
//        {
//          int size_L_t = _L[t].size();
//          for (int d = 0; d < size_L_t; ++d)
//          {
//            if (_zz[i][s][c][t][d].get(GRB_DoubleAttr_X) >= 0.4)
//            {
//              std::cout << _zz[i][s][c][t][d].get(GRB_StringAttr_VarName)
//              << " = " << _zz[i][s][c][t][d].get(GRB_DoubleAttr_X) << std::endl;
//            }
//          }
//        }
//      }
//    }
//  }
//
//  std::cout << std::endl;
//  
//  for (int i = 0; i < nrNodes; ++i)
//  {
//    for (int s = 0; s < nrAnatomicalSites; ++s)
//    {
//      int size_L_s = _L[s].size();
//      for (int c = 0; c < size_L_s; ++c)
//      {
//        if (_x[i][s][c].get(GRB_DoubleAttr_X) >= 0.4)
//        {
//          std::cout << _x[i][s][c].get(GRB_StringAttr_VarName)
//                    << " = " << _x[i][s][c].get(GRB_DoubleAttr_X) << std::endl;
//          std::cout << _r[i][s][c].get(GRB_StringAttr_VarName)
//                    << " = " << _r[i][s][c].get(GRB_DoubleAttr_X) << std::endl;
//        }
//      }
//    }
//  }
//
//  std::cout << std::endl;
//  
//  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
//  {
//    if (!isLeaf(v_i)) continue;
//    const int i = (*_pNodeToIndex)[v_i];
//    if (_x[i][nrAnatomicalSites][0].get(GRB_DoubleAttr_X) >= 0.4)
//    {
//      std::cout << _x[i][nrAnatomicalSites][0].get(GRB_StringAttr_VarName)
//                << " = " << _x[i][nrAnatomicalSites][0].get(GRB_DoubleAttr_X) << std::endl;
//    }
//  }
//
//  std::cout << std::endl;
//
//  for (int s = 0; s < nrAnatomicalSites; ++s)
//  {
//    int size_L_s = _L[s].size();
//    for (int c = 0; c < size_L_s; ++c)
//    {
//      if (_y[s][c].get(GRB_DoubleAttr_X) >= 0.4)
//      {
//        std::cout << _y[s][c].get(GRB_StringAttr_VarName)
//                  << " = " << _y[s][c].get(GRB_DoubleAttr_X) << std::endl;
//      }
//    }
//  }
//  
//  std::cout << std::endl;
//  
//  for (NodeIt v_i(getOrgT().tree()); v_i != lemon::INVALID; ++v_i)
//  {
//    int mapped_i = _F.characterToIndex(getOrgT().label(v_i));
//    for (int p = 0; p < nrSamples; ++p)
//    {
//      std::cout << _f[p][mapped_i].get(GRB_StringAttr_VarName)
//                << " = " << _f[p][mapped_i].get(GRB_DoubleAttr_X) << std::endl;
//      std::cout << _u[p][mapped_i].get(GRB_StringAttr_VarName)
//                << " = " << _u[p][mapped_i].get(GRB_DoubleAttr_X) << std::endl;
//    }
//  }
  
  constructGraph();
  delete _pNodeToStateSet;
  delete _pNodeToRootState;
  _pNodeToStateSet = new IntPairSetNodeMap(getTree());
  _pNodeToRootState = new IntPairNodeMap(getTree(), std::make_pair(-1, -1));
  BoolNodeMap leafPresence(getTree(), true);

  for (int i = 0; i < nrNodes; ++i)
  {
    Node v_i = _indexToNode[i];
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      int size_L_s = _L[s].size();
      for (int c = 0; c < size_L_s; ++c)
      {
        if (_r[i][s][c].get(GRB_DoubleAttr_X) >= 0.4)
        {
          //          std::cout << _r[i][s][c].get(GRB_StringAttr_VarName)
          //                    << " = " << _r[i][s][c].get(GRB_DoubleAttr_X) << std::endl;
          (*_pNodeToRootState)[v_i] = std::make_pair(s, c);
        }
      }
    }
  }
  
  for (int i = 0; i < nrNodes; ++i)
  {
    Node v_i = _indexToNode[i];
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      int size_L_s = _L[s].size();
      for (int c = 0; c < size_L_s; ++c)
      {
        if (_x[i][s][c].get(GRB_DoubleAttr_X) >= 0.4)
        {
//          std::cout << _x[i][s][c].get(GRB_StringAttr_VarName)
//                    << " = " << _x[i][s][c].get(GRB_DoubleAttr_X) << std::endl;
          (*_pNodeToStateSet)[v_i].insert(std::make_pair(s, c));
        }
      }
    }
  }
  
  for (int i = 0; i < nrNodes; ++i)
  {
    Node v_i = _indexToNode[i];
    if (isLeaf(v_i))
    {
      if (_x[i][nrAnatomicalSites][0].get(GRB_DoubleAttr_X) >= 0.4)
      {
//        std::cout << _x[i][nrAnatomicalSites][0].get(GRB_StringAttr_VarName)
//                  << " = " << _x[i][nrAnatomicalSites][0].get(GRB_DoubleAttr_X) << std::endl;
        leafPresence[v_i] = false;
      }
    }
  }
  
  StringToStringMap toMutLabel;
  refine(leafPresence, toMutLabel);
  
  // infer resU and resF
  _pResU = new DoubleVectorNodeMap(_pTprime->tree());
  _pResF = new DoubleVectorNodeMap(_pTprime->tree());
  
  const int nrCharacters = _F.getNrCharacters();
  for (int i = 0; i < nrCharacters; ++i)
  {
    const std::string& label_i = _F.indexToCharacter(i);
    Node v_i = _pTprime->getNodeByLabel(label_i);
    if (v_i == lemon::INVALID)
    {
      // mutation not in mutation tree
      continue;
    }
    assert(!_pTprime->isLeaf(v_i));
    (*_pResF)[v_i] = DoubleVector(nrSamples, 0);

    for (int p = 0; p < nrSamples; ++p)
    {
      double val = _f[p][i].get(GRB_DoubleAttr_X);
      if (g_tol.nonZero(val))
      {
        (*_pResF)[v_i][p] = val;
      }
    }
  }
  
  for (Node v_is : _pTprime->leafSet())
  {
    assert(_pTprime->isLeaf(v_is));
    
    Node v_parent_is = _pTprime->parent(v_is);

    const std::string& sStr = _pTprime->l(v_is);
    const int mapped_s = _F.anatomicalSiteToIndex(sStr);
    const IntSet& samples = _F.anatomicalSiteIndexToSampleIndices(mapped_s);
    
    const std::string& label_v_parent_is = _pTprime->label(v_parent_is);
    assert(toMutLabel.count(label_v_parent_is) == 1);
    const std::string& mut_label = toMutLabel[label_v_parent_is];
    const int mapped_i = _F.characterToIndex(mut_label);
    
    (*_pResU)[v_is] = DoubleVector(samples.size(), 0);
    int idx = 0;
    for (int p : samples)
    {
      double val = _u[p][mapped_i].get(GRB_DoubleAttr_X);
      if (g_tol.nonZero(val))
      {
        (*_pResU)[v_is][idx] = val;
      }
      ++idx;
    }
  }
}

void IlpPmhTiSolver::initIndices()
{
  // construct _extT
  computeFmin(getOrgT().root());
  computeFmax(getOrgT().root());
  computeUmax();
  
  _pNodeToIndex = new IntNodeMap(_extT, -1);
  _pArcToIndex = new IntArcMap(_extT, -1);
  
  constructExtT(getOrgT().root());
  
  IlpPmhTrSolver::initIndices();
}

void IlpPmhTiSolver::computeFmin(Node v_i)
{
  const int k = _F.getNrSamples();
  const int i = _F.characterToIndex(getOrgT().label(v_i));
  
  const Digraph& T = getOrgT().tree();
  
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
        int j = _F.characterToIndex(getOrgT().label(v_j));
        sum_p += _Fmin[p][j];
      }
      _Fmin[p][i] = std::max(_F.min(p, i), sum_p);
    }
  }
}

void IlpPmhTiSolver::computeFmax(Node v_i)
{
  const int k = _F.getNrSamples();
  const int i = _F.characterToIndex(getOrgT().label(v_i));
  const Digraph& T = getOrgT().tree();
  
  if (v_i == getOrgT().root())
  {
    for (int p = 0; p < k; ++p)
    {
      _Fmax[p][i] = _F.max(p, i);
    }
  }
  else
  {
    Node v_pi_i = T.source(InArcIt(T, v_i));
    int pi_i = _F.characterToIndex(getOrgT().label(v_pi_i));
    
    for (int p = 0; p < k; ++p)
    {
      double sum_p = _Fmax[p][pi_i];
      for (OutArcIt a_pi_i_j(T, v_pi_i); a_pi_i_j != lemon::INVALID; ++a_pi_i_j)
      {
        Node v_j = T.target(a_pi_i_j);
        if (v_j == v_i) continue;
        
        int j = _F.characterToIndex(getOrgT().label(v_j));
        sum_p -= _Fmin[p][j];
      }
      
      if (!g_tol.nonZero(sum_p))
      {
        sum_p = 0;
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

void IlpPmhTiSolver::computeUmax()
{
  const int k = _F.getNrSamples();
  const int n = _F.getNrCharacters();
  const Digraph& T = getOrgT().tree();
  
  // compute Umax
  for (int p = 0; p < k; ++p)
  {
    for (int i = 0; i < n; ++i)
    {
      Node v_i = getOrgT().getNodeByLabel(_F.indexToCharacter(i));
      if (v_i == lemon::INVALID || (getOrgT().isLeaf(v_i) && _Fmin[p][i] == 0))
      {
        _Umax[p][i] = 0;
      }
      else
      {
        _Umax[p][i] = _Fmax[p][i];
        for (OutArcIt a_ij(T, v_i); a_ij != lemon::INVALID; ++a_ij)
        {
          Node v_j = T.target(a_ij);
          int j = _F.characterToIndex(getOrgT().label(v_j));
          
          _Umax[p][i] -= _Fmin[p][j];
        }
      }
    }
  }
}

void IlpPmhTiSolver::constructExtT(Node v)
{
  const int m = _F.getNrAnatomicalSites();
  const int k = _F.getNrSamples();
  
  const Digraph& T = getOrgT().tree();
  
  const int mapped_i = _F.characterToIndex(getOrgT().label(v));
  
  NodeVector anatomicalSiteNodes;
  for (int s = 0; s < m; ++s)
  {
    const std::string& sStr = _F.indexToAnatomicalSite(s);
    for (int p : _F.anatomicalSiteIndexToSampleIndices(s))
    {
      assert(0 <= p && p < k);
      if ((getOrgT().isLeaf(v) && _Fmin[p][mapped_i] > 0) ||
          (!getOrgT().isLeaf(v) && _Umax[p][mapped_i] > 0))
      {
        Node v_is = _extT.addNode();
        _isAnatomicalSiteNodeExtT[v_is] = true;
        _extTtoT[v_is] = lemon::INVALID;
        _labelExtT[v_is] = getOrgT().label(v) + "_" + sStr;
        _leafAnatomicalSiteLabelExtT[v_is] = sStr;
        _pNodeToIndex->set(v_is, _indexToNode.size());
        _indexToNode.push_back(v_is);
        anatomicalSiteNodes.push_back(v_is);
        break;
      }
    }
  }

  Node vv = _extT.addNode();
  _isAnatomicalSiteNodeExtT[vv] = false;
  _labelExtT[vv] = getOrgT().label(v);
  _TtoExtT[v] = vv;
  _extTtoT[vv] = v;
  _pNodeToIndex->set(vv, _indexToNode.size());
  _indexToNode.push_back(vv);
  if (v == getOrgT().root())
  {
    _rootExtT = vv;
  }
  
  for (OutArcIt a(T, v); a != lemon::INVALID; ++a)
  {
    Node w = T.target(a);
    constructExtT(w);
    
    _extT.addArc(vv, _TtoExtT[w]);
  }
  
  // attach anatomical sites
  for (Node w : anatomicalSiteNodes)
  {
    _extT.addArc(vv, w);
  }
}

void IlpPmhTiSolver::writeSearchGraphDOT(std::ostream& out) const
{
  const int k = _F.getNrSamples();
  
  out << "digraph extT {" << std::endl;
  
  for (NodeIt v_i(_extT); v_i != lemon::INVALID; ++v_i)
  {
    const int i = (*_pNodeToIndex)[v_i];
    out << "\t" << i << " [label=\"" << _labelExtT[v_i];
    if (!isLeaf(v_i))
    {
      int mapped_i = _F.characterToIndex(getOrgT().label(_extTtoT[v_i]));
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
  
  for (ArcIt a_ij(_extT); a_ij != lemon::INVALID; ++a_ij)
  {
    Node v_i = _extT.source(a_ij);
    Node v_j = _extT.target(a_ij);
    const int i = (*_pNodeToIndex)[v_i];
    const int j = (*_pNodeToIndex)[v_j];
    
    out << "\t" << i << " -> " << j << std::endl;
  }
  
  out << "}" << std::endl;
}

void IlpPmhTiSolver::initConstraints()
{
  IlpPmhTrSolver::initConstraints();
  
  const int nrSamples = _F.getNrSamples();
//  const int nrCharacters = _F.getNrCharacters();
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  
  GRBLinExpr sum, sum2;
  
  for (NodeIt v_i(getOrgT().tree()); v_i != lemon::INVALID; ++v_i)
  {
    int mapped_i = _F.characterToIndex(getOrgT().label(v_i));
    for (int p = 0; p < nrSamples; ++p)
    {
      for (OutArcIt a_ij(getOrgT().tree(), v_i); a_ij != lemon::INVALID; ++a_ij)
      {
        Node v_j = getOrgT().tree().target(a_ij);
        const std::string& label_j = getOrgT().label(v_j);
        int mapped_j = _F.characterToIndex(label_j);
        
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
  
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    if (!isLeaf(v_i)) continue;
    
    const int i = (*_pNodeToIndex)[v_i];
    const std::string& label_i = getLabel(getParent(v_i));
    const int mapped_i = _F.characterToIndex(label_i);
    
    const std::string& sStr = getLeafAnatomicalSiteLabel(v_i);
    assert(_anatomicalSiteToIndex.count(sStr) == 1);

    const int s = _anatomicalSiteToIndex[sStr];
    const int mapped_s = _F.anatomicalSiteToIndex(sStr);
    const int size_L_s = _L[s].size();
    for (int c = 0; c < size_L_s; ++c)
    {
      sum += _x[i][s][c];
    }

    for (int p : _F.anatomicalSiteIndexToSampleIndices(mapped_s))
    {
      _model.addConstr(sum >= _u[p][mapped_i]);
    }
    
    sum.clear();
  }
  
  if (_disablePolytomyResolution)
  {
    // Every inner vertex is labeled by exactly one (s,c)
    for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
    {
      if (isLeaf(v_i)) continue;
      
      const int i = (*_pNodeToIndex)[v_i];
      for (int s = 0; s < nrAnatomicalSites; ++s)
      {
        const int size_L_s = _L[s].size();
        for (int c = 0; c < size_L_s; ++c)
        {
          sum += _x[i][s][c];
        }
      }
      
      _model.addConstr(sum == 1);
      sum.clear();
    }
  }
  
  _model.update();
  
//  for (int i = 0; i < _indexToNode.size(); ++i)
//  {
//    std::cout << i << " : " << getLabel(_indexToNode[i]) << std::endl;
//  }
//  _model.addConstr(_x[21][0][0] == 1);
//  _model.addConstr(_x[21][2][0] == 1);
//  _model.addConstr(_x[21][3][0] == 1);
}

void IlpPmhTiSolver::initSingleSourceSeedingConstraints()
{
  const int nrAnatomicalSites = _indexToAnatomicalSite.size();
  GRBLinExpr sum;
  
  for (Node v_i : _L[_primaryIndex])
  {
    if (v_i == lemon::INVALID) continue;
    assert(isLeaf(v_i));
    
    const int i = (*_pNodeToIndex)[v_i];
    _model.addConstr(_x[i][_primaryIndex][0] + _x[i][nrAnatomicalSites][0] == 1);
  }
  
  for (int t = 0; t < nrAnatomicalSites; ++t)
  {
    if (t == _primaryIndex) continue;
    
    const int size_L_t = _L[t].size();
    for (int d = 0; d < size_L_t; ++d)
    {
      for (int e = 0; e < size_L_t; ++e)
      {
        if (d == e) continue;
        for (int s = 0; s < nrAnatomicalSites; ++s)
        {
          sum += _w[s][t][d][e];
        }
        
        _model.addConstr(sum >= _y[t][d] + _y[t][e] - 1);
        _model.addConstr(sum <= _y[t][d]);
        _model.addConstr(sum <= _y[t][e]);
        sum.clear();
      }
    }
  }
}

void IlpPmhTiSolver::initLeafConstraints()
{
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  
  GRBLinExpr sum, sum2;
  
  // Ensure that each leaf v_i receives a unique label
  // satisfying the definition of a sublabeling
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    if (!isLeaf(v_i)) continue;
    
    const int i = (*_pNodeToIndex)[v_i];
    const std::string& sStr = getLeafAnatomicalSiteLabel(v_i);
    
    assert(_anatomicalSiteToIndex.count(sStr) == 1);
    const int s = _anatomicalSiteToIndex[sStr];
    
    const int size_L_s = _L[s].size();
    for (int c = 0; c < size_L_s; ++c)
    {
      sum += _x[i][s][c];
    }
    _model.addConstr(sum + _x[i][nrAnatomicalSites][0] == 1);
    sum.clear();
    
    for (int t = 0; t < nrAnatomicalSites; ++t)
    {
      if (t == s) continue;
      
      const int size_L_t = _L[t].size();
      for (int c = 0; c < size_L_t; ++c)
      {
        sum += _x[i][t][c];
      }
    }
    _model.addConstr(sum == 0);
    sum.clear();
  }
  
  // one root color
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    if (!isLeaf(v_i)) continue;
    const int i = (*_pNodeToIndex)[v_i];
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      const int size_L_s = _L[s].size();
      for (int c = 0; c < size_L_s; ++c)
      {
        _model.addConstr(_r[i][s][c] <= _x[i][s][c]);
        sum += _r[i][s][c];
      }
    }
    _model.addConstr(sum + _r[i][nrAnatomicalSites][0] == 1);
    _model.addConstr(_r[i][nrAnatomicalSites][0] == _x[i][nrAnatomicalSites][0]);
    sum.clear();
  }
  
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    if (!isLeaf(v_i)) continue;
    const int i = (*_pNodeToIndex)[v_i];
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      const int size_L_s = _L[s].size();
      for (int c = 0; c < size_L_s; ++c)
      {
        sum2 += _x[i][s][c];
        for (int t = 0; t < nrAnatomicalSites; ++t)
        {
          const int size_L_t = _L[t].size();
          for (int d = 0; d < size_L_t; ++d)
          {
            sum += _zz[i][s][c][t][d];
            _model.addConstr(_zz[i][s][c][t][d] <= _x[i][s][c]);
            _model.addConstr(_zz[i][s][c][t][d] <= _x[i][t][d]);
          }
        }
      }
    }
    _model.addConstr(sum2 + _x[i][nrAnatomicalSites][0] - 1 == sum);
    sum.clear();
    sum2.clear();
  }
}

void IlpPmhTiSolver::initVariables()
{
  IlpPmhTrSolver::initVariables();
  
  const int nrSamples = _F.getNrSamples();
  const int nrCharacters = _F.getNrCharacters();
  
  char buf[1024];
  
  /// f[p][i] = frequency of mutation incoming to clone i in sample p
  _f = VarMatrix(nrSamples);
  for (int p = 0; p < nrSamples; ++p)
  {
    _f[p] = VarArray(nrCharacters);
    for (int c = 0; c < nrCharacters; ++c)
    {
      snprintf(buf, 1024, "f;%s;%s",
               _F.indexToSample(p).c_str(),
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
      snprintf(buf, 1024, "u;%s;%s",
               _F.indexToSample(p).c_str(),
               _F.indexToCharacter(c).c_str());
      _u[p][c] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
    }
  }
  
  _model.update();
}

void IlpPmhTiSolver::run(const CloneTree& T,
                         const FrequencyMatrix& F,
                         const std::string& primary,
                         const std::string& outputDirectory,
                         const std::string& outputPrefix,
                         const StringToIntMap& colorMap,
                         MigrationGraph::Pattern pattern,
                         int nrThreads,
                         bool outputILP,
                         bool outputSearchGraph,
                         int timeLimit,
                         const IntTriple& bounds,
                         const StringPairList& forcedComigrations,
                         bool disablePolytomyResolution)
{
  std::string filenameGurobiLog;
  if (!outputDirectory.empty())
  {
    char buf[1024];
    snprintf(buf, 1024, "%s/%slog-%s-%s.txt",
             outputDirectory.c_str(),
             outputPrefix.c_str(),
             primary.c_str(),
             MigrationGraph::getPatternString(pattern).c_str());
    
    filenameGurobiLog = buf;
  }
  
  IlpPmhTiSolver solver(T,
                        F,
                        primary,
                        pattern,
                        filenameGurobiLog,
                        forcedComigrations,
                        disablePolytomyResolution);
  
  IlpPmhSolver::run(solver,
                    T,
                    primary,
                    outputDirectory,
                    outputPrefix,
                    colorMap,
                    pattern,
                    nrThreads,
                    outputILP,
                    outputSearchGraph,
                    timeLimit,
                    bounds,
                    forcedComigrations);
}

GRBLinExpr IlpPmhTiSolver::initObjective(const IntTriple& bounds)
{
  GRBLinExpr obj = IlpPmhTrSolver::initObjective(bounds);
  
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  const int nrNodes = _indexToNode.size();
  const int nrCharacters = _F.getNrCharacters();
  const int nrSamples = _F.getNrSamples();

  const double g = 1000 * (1. / (nrAnatomicalSites * nrAnatomicalSites))
                 * (1. / (nrNodes * nrNodes))
                 * (1. / (nrAnatomicalSites + 1))
                 * (1. / nrNodes);

  
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    const int i = (*_pNodeToIndex)[v_i];
    if (isLeaf(v_i))
    {
      obj += g * (1 - _x[i][nrAnatomicalSites][0]);
    }
  }
  
  // minimize distance
  const double gg = 1000 * (1. / (nrAnatomicalSites * nrAnatomicalSites))
                 * (1. / (nrNodes * nrNodes))
                 * (1. / (nrAnatomicalSites + 1))
                 * (1. / nrNodes)
                 * (1. / (nrCharacters * nrSamples + 1));
  
  for (int i = 0; i < nrCharacters; ++i)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      obj += gg * (_Umax[p][i] - _u[p][i]);
    }
  }
  
  _model.setObjective(obj, GRB_MINIMIZE);
  _model.update();
  
  return obj;
}

void IlpPmhTiSolver::writeCloneTree(std::ostream& out,
                                    const StringToIntMap& colorMap) const
{
  T().writeDOT(out, lPlus(), colorMap, getU(), T().getIdMap());
}
