/*
 * ilpsolver.h
 *
 *  Created on: 2-feb-2017
 *      Author: M. El-Kebir
 */

#ifndef ILPSOLVER_H
#define ILPSOLVER_H

#include "utils.h"
#include "nonbinaryclonetree.h"
#include "migrationgraph.h"
#include <gurobi_c++.h>

/// This class implements an ILP for solving the Parsimonious Migration History
/// (PMH) problem under various topological constraints (PS, S, M and R)
///
/// \brief This class implements an ILP for solving PMH under topological
/// constraints
class IlpSolver
{
public:
  /// Constructor
  ///
  /// @param T Non-binary clone tree
  /// @param primary Primary tumor
  /// @param pattern Topological constraint
  /// @param gurobiLogFilename Gurobi logging filename
  /// @param forcedComigrations List of ordered pairs of anatomical sites
  /// that must be present
  IlpSolver(const NonBinaryCloneTree& T,
            const std::string& primary,
            MigrationGraph::Pattern pattern,
            const std::string& gurobiLogFilename,
            const StringPairList& forcedComigrations);
  
  /// Destructor
  virtual ~IlpSolver()
  {
    delete _pLPlus;
    delete _pArcToIndex;
    delete _pNodeToIndex;
  };
  
  /// Solve ILP
  ///
  /// @param nrThreads Number of threads
  /// @param timeLimit Time limit in seconds (-1: no time limit)
  bool solve(int nrThreads, int timeLimit);
  
  /// Export ILP model
  void exportModel(const std::string& filename);
  
  /// Return vertex labeling
  virtual const StringNodeMap& lPlus() const
  {
    return *_pLPlus;
  }
  
  /// Initialize ILP
  virtual void init(double upperBound);
  
  /// Return lower bound
  double LB() const
  {
    return _LB;
  }
  
  /// Return upper bound
  double UB() const
  {
    return _UB;
  }
  
protected:
  /// Initialize indices and mappings
  virtual void initIndices();
  
  /// Initialize ILP variables
  virtual void initVariables();
  
  /// Initialize ILP constraints on the leaves
  virtual void initLeafConstraints();
  
  /// Initialize ILP constraints
  virtual void initConstraints();
  
  /// Initialize ILP single source seeding constraints
  virtual void initSingleSourceSeedingConstraints();
  
  /// Initialize ILP parallel single source seeding constraints
  virtual void initParallelSingleSourceSeedingConstraints();
  
  /// Initialize ILP multi source seeding constraints
  virtual void initMultiSourceSeedingConstraints();
  
  /// Initialize ILP objective
  ///
  /// @param upperBound Upper bound on objective function
  virtual void initObjective(double upperBound);
  
  /// Initialize ILP constraints regarding primary tumor anatomical site
  /// (only applicable for PS and S)
  ///
  /// @param v_j Node in T
  virtual void initPrimaryConstraint(Node v_j);
  
  /// Initialized ILP constraints regarding forced comigrations
  virtual void initForcedComigrations();
  
  /// Process ILP solution
  virtual void processSolution();
  
  /// Add ILP contraint involving matching anatomical sites of adjacent vertices
  ///
  /// @param sum_z Sum of z variables
  /// @param ij Index of arc (v_i, v_j)
  virtual void addMatchingColorsConstraint(GRBLinExpr sum_z, int ij)
  {
    _model.addConstr(sum_z + _y[ij] == 1);
  }
  
  /// Add ILP comigration constraint
  ///
  /// @param s Source anatomical site
  /// @param t Target anatomical site
  /// @param ij Arc (v_i, v_j) index
  /// @param i Node v_i index
  /// @param j Node v_j index
  virtual void addComigrationConstraint(int s, int t, int ij, int i, int j)
  {
    _model.addConstr(_c[s][t] - _x[i][s] - _x[j][t] >= -1);
  }
  
  /// Return underlying LEMON tree
  virtual const Digraph& tree() const
  {
    return _T.tree();
  }
  
  /// Return node identifier
  ///
  /// @param v Node in T
  virtual const std::string& label(Node v) const
  {
    return _T.label(v);
  }
  
  /// Return leaf anatomical site label
  ///
  /// @param v Node in T
  virtual const std::string& l(Node v) const
  {
    assert(_T.isLeaf(v));
    return _T.l(v);
  }
  
  /// Return root node
  virtual Node root() const
  {
    return _T.root();
  }
  
  /// Return whether given node is a leaf
  ///
  /// @param v Node in T
  virtual bool isLeaf(Node v) const
  {
    return _T.isLeaf(v);
  }
  
  /// Vector of arcs
  typedef std::vector<Arc> ArcVector;
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
  /// Gurobi variable 6D matrix
  typedef std::vector<Var5Matrix> Var6Matrix;
  
  /// Generate next combination of metastases
  static bool next(BoolVector& mets);
  
protected:
  /// Non-binary clone tree
  const NonBinaryCloneTree& _T;
  /// Label of primary tumor anatomical site
  const std::string& _primary;
  /// Topological constraint
  const MigrationGraph::Pattern _pattern;
  /// List of ordered pairs of anatomical sites that must be present
  const StringPairList& _forcedComigrations;
  /// Index to arc map
  ArcVector _indexToArc;
  /// Map to arc index
  IntArcMap* _pArcToIndex;
  /// Index to node map
  NodeVector _indexToNode;
  /// Node to index map
  IntNodeMap* _pNodeToIndex;
  /// Sample label to index map
  StringToIntMap _sampleToIndex;
  /// Index to sample label map
  StringVector _indexToSample;
  /// Primary tumor index
  int _primaryIndex;
  /// Gurobi environment
  GRBEnv _env;
  /// Gurobi model
  GRBModel _model;
  /// Variable x[i][s] = 1 iff vertex v_i is labeled by sample s
  VarMatrix _x;
  /// Variable y[ij] = 1 iff edge (v_i, v_j) is a migration edge
  VarArray _y;
  /// Variable z[ij][s] = 1 iff v_i and v_j are labeled by sample s
  VarMatrix _z;
  /// Variable c[s][t] = 1 iff there exists a migration edge (v_i, v_j)
  /// where l(v_i) = s and l(v_j) = t
  VarMatrix _c;
  /// Variable d[s] = 1 iff there exists a migration edge (v_i, v_j)
  /// where l(v_i) = s
  VarArray _d;
  /// Resulting vertex labeling
  StringNodeMap* _pLPlus;
  /// Lowerbound on the optimal solution
  double _LB;
  /// Upperbound on the optimal solution
  double _UB;
};

#endif // ILPSOLVER_H
