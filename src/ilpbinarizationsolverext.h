/*
 * ilpbinarizationsolverext.h
 *
 *  Created on: 15-apr-2017
 *      Author: M. El-Kebir
 */

#ifndef ILPSOLVERBINARIZATIONEXT_H
#define ILPSOLVERBINARIZATIONEXT_H

#include "utils.h"
#include "ilpsolverext.h"

class IlpBinarizationSolverExt : public IlpSolverExt
{
public:
  IlpBinarizationSolverExt(const NonBinaryCloneTree& T,
                           const FrequencyMatrix& F,
                           const std::string& primary,
                           Mode mode,
                           const std::string& gurobiLogFilename,
                           const StringPairList& forcedComigrations);
  
protected:
  virtual void initVariables();
  virtual void initConstraints();
  virtual void processSolution();
  
  virtual void addMatchingColorsConstraint(GRBLinExpr sum_z, int ij)
  {
    _model.addConstr(sum_z + _y[ij] == _w[ij]);
  }
  
  virtual void addMatchingColorsAndUsageConstraint(int ij,
                                                   int mapped_i,
                                                   int t)
  {
    const int nrSamples = _indexToSample.size();
    _model.addConstr(_z[ij][nrSamples] <= 1 - _u[t][mapped_i]);
  }
  
  virtual void addComigrationConstraint(int s, int t,
                                        int ij, int i, int j)
  {
    _model.addConstr(_c[s][t] - _x[i][s] - _x[j][t] >= -1 - (1-_w[ij]));
  }
  
  virtual void constructG();
  
  void constructG(Node v);
  
  int getCharacterIndex(const Digraph& resT,
                        const IntNodeMap& resNodeToIndex,
                        Node v_i) const
  {
    int i = -1;
    if (OutArcIt(resT, v_i) == lemon::INVALID)
    {
      Node parent = resT.source(InArcIt(resT, v_i));
      i = resNodeToIndex[parent];
    }
    else
    {
      i = resNodeToIndex[v_i];
    }
    
    Node org_v_i = _indexToNode[i];
    Node tree_v_i  = _GtoT[org_v_i];
    
    int mapped_i = _F.characterToIndex(_T.label(tree_v_i));
    return mapped_i;
  }
  
  void labelEdges(const Digraph& resT,
                  const IntNodeMap& resNodeToIndex,
                  Node v_i,
                  IntNodeMap& characterLabel,
                  BoolVector& visited);
  
  typedef Digraph::NodeMap<Arc> ArcNodeMap;
  
protected:
  BoolNodeMap _isSampleNode;
  ArcNodeMap _sampleNodeToArcMap;
  /// w[ij] = 1 iff edge (v_i, v_j) is in the tree
  VarArray _w;
};

#endif // ILPSOLVERBINARIZATIONEXT_H
