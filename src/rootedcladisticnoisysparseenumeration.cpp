/*
 * rootedcladisticnoisysparseenumeration.cpp
 *
 *  Created on: 20-aug-2017
 *      Author: M. El-Kebir
 */

#include "rootedcladisticnoisysparseenumeration.h"

namespace gm {

RootedCladisticNoisySparseEnumeration::RootedCladisticNoisySparseEnumeration(const RootedCladisticNoisyAncestryGraph& G,
                                      int limit,
                                      int timeLimit,
                                      int threads,
                                      int lowerbound,
                                      bool monoclonal,
                                      bool fixTrunk,
                                      const IntSet& whiteList)
  : RootedCladisticNoisyEnumeration(G, limit, timeLimit, threads, lowerbound, monoclonal, fixTrunk, whiteList)
  , _env()
  , _model(_env)
  , _f()
  , _u()
  , _z()
  , _vertexToIndex(G.G(), -1)
{
};
  
void RootedCladisticNoisySparseEnumeration::initF(int solIdx,
                                                  RealTensor& F) const
{
  const int m = _G.F().m();
  const int n = _G.F().n();
  const int k = _G.F().k();
  
  const Digraph& G = _G.G();
  
  // 1. Initialize solution tree
  BoolNodeMap filterNodesT(G, false);
  BoolArcMap filterArcsT(G, false);
  SubDigraph T(G, filterNodesT, filterArcsT);
  
  ArcListList::const_iterator it = _result.begin();
  int tmp = solIdx;
  for (; tmp > 0; --tmp, ++it);
  
  const ArcList& arcs = *it;
  for (ArcListIt it = arcs.begin(); it != arcs.end(); ++it)
  {
    Arc a_cidj = *it;
    Node v_ci = G.source(a_cidj);
    Node v_dj = G.target(a_cidj);
    T.enable(v_ci);
    T.enable(v_dj);
    T.enable(a_cidj);
  }
  
  // 2. Formulate ILP
  initVariables(T);
  initConstraints(T);
  initObjective(T);
  
  // 3. Solve
  if (g_verbosity != VERBOSE_DEBUG)
  {
    _model.getEnv().set(GRB_IntParam_LogToConsole, 0);
  }
  _model.optimize();
  
  int status = _model.get(GRB_IntAttr_Status);
  if (status == GRB_OPTIMAL)
  {
    for (int i = 0; i < k; ++i)
    {
      for (int p = 0; p < m; ++p)
      {
        for (int c = 0; c < n; ++c)
        {
          double val = _f[i][p][c].get(GRB_DoubleAttr_X);
          assert(0 <= val && val <= 1);
          F.set(i, p, c, val);
        }
      }
    }
  }
  else
  {
    isValid(T, _G.root(), F);
  }
  
  for (int p = 0; p < F.m(); ++p)
  {
    F.setRowLabel(p, _G.F().getRowLabel(p));
  }
  for (int c = 0; c < F.n(); ++c)
  {
    F.setColLabel(c, _G.F().getColLabel(c));
  }
}
  
void RootedCladisticNoisySparseEnumeration::initVariables(const SubDigraph& T) const
{
  const int m = _G.F().m();
  const int n = _G.F().n();
  const int k = _G.F().k();
  const int nrVertices = lemon::countNodes(T);
  
  char buf[1024];
  
  /// f[i][p][c]
  _f = Var3Matrix(k);
  for (int i = 0; i < k; ++i)
  {
    _f[i] = VarMatrix(m);
    for (int p = 0; p < m; ++p)
    {
      _f[i][p] = VarArray(n);
      for (int c = 0; c < n; ++c)
      {
        snprintf(buf, 1024, "f_%d_%s_%s",
                 i,
                 _G.F().getRowLabel(p).c_str(),
                 _G.F().getColLabel(c).c_str());
        _f[i][p][c] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
      }
    }
  }

  int idx = 0;
  for (SubNodeIt v(T); v != lemon::INVALID; ++v)
  {
    _vertexToIndex[v] = idx;
    ++idx;
  }
  
  /// u[p][v] = usage of vertex v in sample p
  _u = VarMatrix(m);
  _z = VarMatrix(m);
  for (int p = 0; p < m; ++p)
  {
    _u[p] = VarArray(nrVertices);
    _z[p] = VarArray(nrVertices);
    for (SubNodeIt v(T); v != lemon::INVALID; ++v)
    {
      int idx = _vertexToIndex[v];
      assert(idx < nrVertices);
      snprintf(buf, 1024, "u_%s_%d",
               _G.F().getRowLabel(p).c_str(),
               idx);
      _u[p][idx] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
      
      snprintf(buf, 1024, "z_%s_%d",
               _G.F().getRowLabel(p).c_str(),
               idx);
      _z[p][idx] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
    }
  }
  
  _model.update();
}
  
void RootedCladisticNoisySparseEnumeration::initConstraints(const SubDigraph& T) const
{
  const int m = _G.F().m();
  const int n = _G.F().n();
  const int k = _G.F().k();
  
  for (int i = 0; i < k; ++i)
  {
    for (int p = 0; p < m; ++p)
    {
      for (int c = 0; c < n; ++c)
      {
        double lb = _noisyG.F_lb()(i, p, c);
        double ub = _noisyG.F_ub()(i, p, c);
        _model.addConstr(_f[i][p][c] >= std::min(lb, ub));
        _model.addConstr(_f[i][p][c] <= std::max(lb, ub));
      }
    }
  }
  
  for (SubNodeIt v(T); v != lemon::INVALID; ++v)
  {
    int idx = _vertexToIndex[v];
    for (int p = 0; p < m; ++p)
    {
      _model.addConstr(_u[p][idx] >= 0);
      _model.addConstr(_u[p][idx] <= 1);
    }
  }
  
  // sum condition
  for (SubNodeIt v_ci(T); v_ci != lemon::INVALID; ++v_ci)
  {
    int idx = _vertexToIndex[v_ci];
    for (int p = 0; p < m; ++p)
    {
      IntPair ci = *(_G.nodeToCharState(v_ci).begin());
      
      GRBLinExpr lhs;
      for (int l : _G.S(ci.first).D(ci.second))
      {
        lhs += _f[l][p][ci.first];
      }
      
      GRBLinExpr rhs;
      for (SubOutArcIt a(T, v_ci); a != lemon::INVALID; ++a)
      {
        Node v_dj = T.target(a);
        IntPair dj = *(_G.nodeToCharState(v_dj).begin());
        
        for (int l : _G.S(dj.first).D(dj.second))
        {
          rhs += _f[l][p][dj.first];
        }
      }
      
//      _model.addConstr(lhs >= rhs);
      _model.addConstr(_u[p][idx] == lhs - rhs);
      _model.addConstr(_z[p][idx] >= _u[p][idx]);
    }
  }
  
  _model.update();
}
  
void RootedCladisticNoisySparseEnumeration::initObjective(const SubDigraph& T) const
{
  const int m = _G.F().m();
  const int nrVertices = lemon::countNodes(T);
  
  GRBLinExpr obj;
  GRBLinExpr usages;
  
  for (SubNodeIt v_ci(T); v_ci != lemon::INVALID; ++v_ci)
  {
    int idx = _vertexToIndex[v_ci];
    for (int p = 0; p < m; ++p)
    {
      obj += _z[p][idx];
      obj -= 1. / (nrVertices * m) * _u[p][idx];
    }
  }
  
  _model.setObjective(obj, GRB_MINIMIZE);
  _model.update();
}
  
} // namespace gm
