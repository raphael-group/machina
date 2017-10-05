/*
 * rootedcladisticnoisysparseenumeration.h
 *
 *  Created on: 20-aug-2017
 *      Author: M. El-Kebir
 */

#ifndef ROOTEDCLADISTICNOISYSPARSEENUMERATION_H
#define ROOTEDCLADISTICNOISYSPARSEENUMERATION_H

#include "spruce/rootedcladisticnoisyenumeration.h"
#include <gurobi_c++.h>

namespace gm {
  
class RootedCladisticNoisySparseEnumeration : public RootedCladisticNoisyEnumeration
{
public:
  RootedCladisticNoisySparseEnumeration(const RootedCladisticNoisyAncestryGraph& G,
                                  int limit,
                                  int timeLimit,
                                  int threads,
                                  int lowerbound,
                                  bool monoclonal,
                                  bool fixTrunk,
                                        const IntSet& whiteList);
private:
  virtual void initF(int solIdx, RealTensor& F) const;
  /// Gurobi variable array
  typedef std::vector<GRBVar> VarArray;
  /// Gurobi variable matrix
  typedef std::vector<VarArray> VarMatrix;
  /// Gurobi variable 3D matrix
  typedef std::vector<VarMatrix> Var3Matrix;
  
  void initVariables(const SubDigraph& T) const;
  
  void initConstraints(const SubDigraph& T) const;
  
  void initObjective(const SubDigraph& T) const;

private:
  /// Gurobi environment
  mutable GRBEnv _env;
  /// Gurobi model
  mutable GRBModel _model;
  /// Variable f[p][c] is frequency of mutation incoming to clone i in sample p
  mutable Var3Matrix _f;
  /// Variable u[p][i] is usage of clone i in sample p
  mutable VarMatrix _u;
  /// Variable z[p][i] indicates whether u[p][i] > 0
  mutable VarMatrix _z;
  /// Vertex to index map
  mutable IntNodeMap _vertexToIndex;
};

} // namespace gm

#endif // #ROOTEDCLADISTICNOISYSPARSEENUMERATION_H
