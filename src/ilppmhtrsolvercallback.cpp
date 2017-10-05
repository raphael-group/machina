/*
 * ilppmhtrsolvercallback.cpp
 *
 *  Created on: 15-sep-2017
 *      Author: M. El-Kebir
 */

#include "ilppmhtrsolvercallback.h"
#include <lemon/connectivity.h>

IlpPmhPrSolverCycleElimination::IlpPmhPrSolverCycleElimination(const StringVector& indexToAnatomicalSite,
                                                               int primaryIndex,
                                                               const VarMatrix& y,
                                                               const Var4Matrix& z)
  : _indexToAnatomicalSite(indexToAnatomicalSite)
  , _primaryIndex(primaryIndex)
  , _y(y)
  , _z(z)
  , _G()
  , _filterNodes(_G)
  , _filterArcs(_G)
  , _subG(_G, _filterNodes, _filterArcs)
  , _rootG(lemon::INVALID)
  , _stateToNodeG()
  , _nodeToState(_G)
  , _sccMap(_G)
  , _arcLookUp(_G)
  , _costMap(_G)
  , _mmc(_subG, _costMap)
{
  // initialize graph
  const int nrAnatomicalSites = _indexToAnatomicalSite.size();
  
  _stateToNodeG = NodeMatrix(nrAnatomicalSites);
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    const int size_L_s = _z[s].size();
    _stateToNodeG[s] = NodeVector(size_L_s, lemon::INVALID);
    for (int c = 0; c < size_L_s; ++c)
    {
      Node v_sc = _G.addNode();
      _nodeToState[v_sc] = std::make_pair(s, c);
      _stateToNodeG[s][c] = v_sc;
      _filterNodes[v_sc] = false;
      
      if (s == _primaryIndex && c == 0)
      {
        _rootG = v_sc;
      }
    }
  }
}

void IlpPmhPrSolverCycleElimination::updateG()
{
  const int nrAnatomicalSites = _indexToAnatomicalSite.size();
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    const int size_L_s = _z[s].size();
    for (int c = 0; c < size_L_s; ++c)
    {
      Node v_sc = _stateToNodeG[s][c];
      double sol = getSolution(_y[s][c]);
      _filterNodes[v_sc] = (sol >= 0.4);
    }
  }
  
  lemon::mapFill(_G, _filterArcs, false);
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    const int size_L_s = _z[s].size();
    for (int c = 0; c < size_L_s; ++c)
    {
      Node v_sc = _stateToNodeG[s][c];

      for (int t = 0; t < nrAnatomicalSites; ++t)
      {
        const int size_L_t = _z[t].size();
        for (int d = 0; d < size_L_t; ++d)
        {
          Node v_td = _stateToNodeG[t][d];
          double sol = getSolution(_z[s][c][t][d]);
          if (sol >= 0.4)
          {
            Arc a_sctd = _arcLookUp(v_sc, v_td);
            if (a_sctd == lemon::INVALID)
            {
              a_sctd = _G.addArc(v_sc, v_td);
              _costMap[a_sctd] = 1;
            }
            _filterArcs[a_sctd] = true;
//            std::cout << "" << _indexToAnatomicalSite[s] << "_" << c <<" -> " << _indexToAnatomicalSite[t] << "_" << d << "" << std::endl;
          }
        }
      }
    }
  }
}

void IlpPmhPrSolverCycleElimination::callback()
{
  if (where == GRB_CB_MIPSOL)
  {
    updateG();
    
    if (_mmc.run())
    {
      GRBLinExpr sum;
      
      const lemon::Path<SubDigraph>& cycle = _mmc.cycle();
      const int n = cycle.length();
      if (n > 1)
      {
        for (int i = 0; i < n; ++i)
        {
          Arc a_sctd = cycle.nth(i);
          Node v_sc = _G.source(a_sctd);
          Node v_td = _G.target(a_sctd);
          const IntPair& sc = _nodeToState[v_sc];
          const IntPair& td = _nodeToState[v_td];
          
//          std::cout << "" << _indexToAnatomicalSite[sc.first] << "_" << sc.second <<" -> " << _indexToAnatomicalSite[td.first] << "_" << td.second << "" << "\t";
          sum += _z[sc.first][sc.second][td.first][td.second];
        }
//        std::cout << std::endl;
      }
      addLazy(sum <= n - 1);
      sum.clear();
    }
  }
}
