/*
 * ilpsolverext.h
 *
 *  Created on: 13-apr-2017
 *      Author: M. El-Kebir
 */

#ifndef ILPSOLVEREXT_H
#define ILPSOLVEREXT_H

#include "utils.h"
#include "frequencymatrix.h"
#include "nonbinaryclonetree.h"
#include <gurobi_c++.h>

/// This class implements an ILP for solving PMH-CTI under topological constraints
class IlpSolverExt
{
public:
  const static int _nrModes = 4;
  
  enum Mode {
    PARALLEL_SINGLE_SOURCE_SEEDING,
    SINGLE_SOURCE_SEEDING,
    MULTI_SOURCE_SEEDING,
    RESEEDING
  };
  
  IlpSolverExt(const NonBinaryCloneTree& T,
               const FrequencyMatrix& F,
               const std::string& primary,
               Mode mode,
               const std::string& gurobiLogFilename,
               const StringPairList& forcedComigrations);
  
  virtual ~IlpSolverExt();
  
  bool solve(int nrThreads, int timeLimit);
  
  void exportModel(const std::string& filename);
  
  virtual const StringNodeMap& lPlus() const
  {
    return *_pResLPlus;
  }
  
  virtual const NonBinaryCloneTree& T() const
  {
    return *_pResCloneTree;
  }
  
  virtual void init(double upperBound);
  
  double LB() const
  {
    return _LB;
  }
  
  double UB() const
  {
    return _UB;
  }
  
  const DoubleVectorNodeMap& getF() const
  {
    return *_pResF;
  }
  
  const DoubleNodeMap& getU() const
  {
    return *_pResU;
  }
  
  const IntNodeMap& getCharacterLabel() const
  {
    return *_pResCharacterLabel;
  }
  
  void writeDOT(std::ostream& out) const;
  
  void writeDOT(std::ostream& out,
                const StringToIntMap& colorMap) const;
  
protected:
  virtual void constructG();
  virtual void initIndices();
  virtual void initVariables();
  virtual void initLeafConstraints();
  virtual void initConstraints();
  virtual void initSingleSourceSeedingConstraints();
  virtual void initParallelSingleSourceSeedingConstraints();
  virtual void initMultiSourceSeedingConstraints();
  virtual void initObjective(double upperBound);
  virtual void initPrimaryConstraint(Node v_j);
  virtual void initForcedComigrations();
  virtual void processSolution();

  virtual int getNrBackBoneVertices() const
  {
    return _F.getNrCharacters();
  }
  
  virtual void addMatchingColorsConstraint(GRBLinExpr sum_z, int ij)
  {
    _model.addConstr(sum_z + _y[ij] == 1);
  }
  
  virtual void addMatchingColorsAndUsageConstraint(int ij, int mapped_i, int t)
  {
    const int nrSamples = _indexToSample.size();
    _model.addConstr(_z[ij][nrSamples] <= 1 - _u[t][mapped_i]);
  }
  
  virtual void addComigrationConstraint(int s, int t,
                                        int ij, int i, int j)
  {
    _model.addConstr(_c[s][t] - _x[i][s] - _x[j][t] >= -1);
  }
  
  virtual const std::string& l(Node v) const
  {
    assert(isLeaf(v));
    return _l[v];
  }
  
  virtual const std::string& label(Node v) const
  {
    return _label[v];
  }
   
  virtual Node root() const
  {
    return _root;
  }
  
  virtual bool isLeaf(Node v) const
  {
    return OutArcIt(_G, v) == lemon::INVALID;
  }
  
  typedef std::vector<Arc> ArcVector;
  
  typedef std::vector<GRBVar> VarArray;
  typedef std::vector<VarArray> VarMatrix;
  typedef std::vector<VarMatrix> Var3Matrix;
  typedef std::vector<Var3Matrix> Var4Matrix;
  typedef std::vector<Var4Matrix> Var5Matrix;
  typedef std::vector<Var5Matrix> Var6Matrix;

  static bool next(BoolVector& mets);
  
  void computeFmin(Node v_i);
  void computeFmax(Node v_i);
  void computeUmax();
  
protected:
  const NonBinaryCloneTree& _T;
  const FrequencyMatrix& _F;
  DoubleMatrix _Fmin;
  DoubleMatrix _Fmax;
  DoubleMatrix _Umax;
  const std::string& _primary;
  const Mode _mode;
  const StringPairList& _forcedComigrations;
  
  Digraph _G;
  Node _root;
  StringNodeMap _label;
  StringNodeMap _l;
  StringNodeMap _lPlus;
  NodeNodeMap _TtoG;
  NodeNodeMap _GtoT;
  
  ArcVector _indexToArc;
  IntArcMap _arcToIndex;
  NodeVector _indexToNode;
  IntNodeMap _nodeToIndex;
  StringToIntMap _sampleToIndex;
  StringVector _indexToSample;
  int _primaryIndex;
  
  GRBEnv _env;
  GRBModel _model;
  /// x[i][s] = 1 iff vertex v_i is labeled by sample s
  VarMatrix _x;
  /// y[ij] = 1 iff edge (v_i, v_j) is a migration edge
  VarArray _y;
  /// z[ij][s] = 1 iff v_i and v_j are labeled by sample s
  VarMatrix _z;
  /// c[s][t] = 1 iff there exists a migration edge (v_i, v_j) where l(v_i) = s and l(v_j) = t
  VarMatrix _c;
  /// d[s] = 1 iff there exists a migration edge (v_i, v_j) where l(v_i) = s
  VarArray _d;
  /// f[s][i] = frequency of mutation incoming to clone i in sample s
  VarMatrix _f;
  /// u[s][i] = usage of clone i in sample s
  VarMatrix _u;
  
  double _LB;
  double _UB;
  
  // Result
  DoubleVectorNodeMap _resF;
  DoubleVectorNodeMap _resU;
  NonBinaryCloneTree* _pResCloneTree;
  StringNodeMap* _pResLPlus;
  DoubleVectorNodeMap* _pResF;
  DoubleNodeMap* _pResU;
  IntNodeMap* _pResCharacterLabel;
};

#endif // ILPSOLVEREXT_H
