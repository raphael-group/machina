/*
 * ilppmhprsolvercallback.h
 *
 *  Created on: 15-sep-2017
 *      Author: M. El-Kebir
 */

#ifndef ILPPMHPRSOLVERCYCLEELIMNATION
#define ILPPMHPRSOLVERCYCLEELIMNATION

#include "gurobi_c++.h"
#include "utils.h"
#include <lemon/hartmann_orlin_mmc.h>

class IlpPmhPrSolverCycleElimination : public GRBCallback
{
public:
  /// Vector of node sets
  typedef std::vector<NodeSet> NodeSetVector;
  /// Gurobi variable array
  typedef std::vector<GRBVar> VarArray;
  /// Gurobi variable matrix
  typedef std::vector<VarArray> VarMatrix;
  /// Gurobi variable 3D matrix
  typedef std::vector<VarMatrix> Var3Matrix;
  /// Gurobi variable 4D matrix
  typedef std::vector<Var3Matrix> Var4Matrix;
  /// Gurobi variable 5D matrix
  typedef std::vector<Var4Matrix> Var5Matrix;
  /// Node matrix
  typedef std::vector<NodeVector> NodeMatrix;
  
public:
  IlpPmhPrSolverCycleElimination(const StringVector& indexToAnatomicalSite,
                                 int primaryIndex,
                                 const VarMatrix& y,
                                 const Var4Matrix& z);
  
protected:
  void callback();
  
  void updateG();
  
private:
  const StringVector& _indexToAnatomicalSite;
  const int _primaryIndex;
  const VarMatrix& _y;
  const Var4Matrix& _z;
  Digraph _G;
  BoolNodeMap _filterNodes;
  BoolArcMap _filterArcs;
  SubDigraph _subG;
  Node _rootG;
  NodeMatrix _stateToNodeG;
  IntPairNodeMap _nodeToState;
  IntNodeMap _sccMap;
  lemon::DynArcLookUp<Digraph> _arcLookUp;
  IntArcMap _costMap;
  lemon::HartmannOrlinMmc<SubDigraph, IntArcMap> _mmc;
};

#endif // ILPPMHPRSOLVERCYCLEELIMNATION
