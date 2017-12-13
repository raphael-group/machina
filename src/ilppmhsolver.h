/*
 * ilppmhsolver.h
 *
 *  Created on: 25-aug-2017
 *      Author: M. El-Kebir
 */

#ifndef ILPPMHSOLVER_H
#define ILPPMHSOLVER_H

#include "utils.h"
#include "clonetree.h"
#include "migrationgraph.h"
#include <gurobi_c++.h>

/// This class implements an ILP for solving the Parsimonious Migration History
/// (PMH) problem under various topological constraints (PS, S, M and R)
///
/// \brief This class implements an ILP for solving PMH under topological
/// constraints
class IlpPmhSolver
{
public:
  /// Constructor
  ///
  /// @param T Clone tree
  /// @param primary Primary tumor
  /// @param pattern Topological constraint
  /// @param gurobiLogFilename Gurobi logging filename
  /// @param forcedComigrations List of ordered pairs of anatomical sites
  /// that must be present
  IlpPmhSolver(const CloneTree& T,
               const std::string& primary,
               MigrationGraph::Pattern pattern,
               const std::string& gurobiLogFilename,
               const StringPairList& forcedComigrations);
  
  /// Destructor
  virtual ~IlpPmhSolver();
  
  /// Return clone tree
  virtual const CloneTree& T() const
  {
    return _T;
  }
  
  /// Return vertex labeling
  virtual const StringNodeMap& lPlus() const
  {
    return *_pLPlus;
  }
  
  /// Solve ILP
  ///
  /// @param nrThreads Number of threads
  /// @param timeLimit Time limit in seconds (-1: no time limit)
  bool solve(int nrThreads, int timeLimit);
  
  /// Export ILP model
  void exportModel(const std::string& filename);
  
  /// Initialize ILP
  ///
  /// @param bounds Upper bounds on mu, gamma and sigma
  virtual void init(const IntTriple& bounds);
  
  /// Return lower bound
  double LB() const
  {
    return _LB;
  }
  
  /// Write solution graph with sublabels in DOT format
  ///
  /// @param out Output stream
  /// @param colorMap Color map
  void writeSolutionGraphDOT(std::ostream& out,
                             const StringToIntMap& colorMap) const;
  
  /// Write clone tree in DOT format
  ///
  /// @param out Output stream
  /// @param colorMap Color map
  virtual void writeCloneTree(std::ostream& out,
                              const StringToIntMap& colorMap) const;
  
  /// Return upper bound
  double UB() const
  {
    return _UB;
  }
  
  /// Solve PMH under a topological constraint
  ///
  /// @param T Non-binary clone tree
  /// @param primary Primary tumor
  /// @param outputDirectory Output directory
  /// @param outputPrefix Prefix prepended to every output filename
  /// @param colorMap Color map
  /// @param pattern Topological constraint
  /// @param nrThreads Number of threads
  /// @param outputILP Output ILP model
  /// @param outputSearchGraph Output search graph
  /// @param timeLimit Time limit in seconds
  /// @param bounds Upper bounds on mu, gamma and sigma
  /// @param forcedComigrations List of ordered pairs of anatomical sites
  /// that must be present
  static IntTriple run(const CloneTree& T,
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
                       const StringPairList& forcedComigrations);
  
protected:
  /// Initialize indices and mappings
  virtual void initIndices();
  
  /// Initialize ILP variables
  virtual void initVariables();
  
  /// Initialize ILP leaf variables
  virtual void initLeafVariables();
  
  /// Initialize ILP callbacks
  virtual void initCallbacks();
  
  /// Initialize ILP constraints regarding vertex labelings
  virtual void initVertexLabelingConstraints();
  
  /// Initialized ILP constraints regarding forced comigrations
  virtual void initForcedComigrations();
  
  /// Initialize ILP leaf constraints
  virtual void initLeafConstraints();
  
  /// Initialize ILP constraints regarding _G
  virtual void initConstraintsG();
  
  /// Initialize ILP constraints regarding nonedges of _G
  virtual void initConstraintsNonEdgesG();
  
  /// Initialize ILP constraints
  virtual void initConstraints();
  
  /// Intialize ILP M constraints
  virtual void initMultiSourceSeedingConstraints();
  
  /// Initialize ILP S constraints
  virtual void initSingleSourceSeedingConstraints();
  
  /// Initialize initial solution
  virtual void initWarmStart();
  
  /// Initialize ILP objective function
  ///
  /// @param bounds Upper bounds on mu, gamma and sigma
  virtual GRBLinExpr initObjective(const IntTriple& bounds);
  
  /// Process ILP solution
  virtual void processSolution();
  
  /// Return underlying tree
  virtual const Digraph& getTree() const
  {
    return _T.tree();
  }
  
  /// Return root vertex of underlying tree
  virtual Node getRoot() const
  {
    return _T.root();
  }
  
  /// Return whether given node of underlying tree is a leaf
  virtual bool isLeaf(Node v) const
  {
    return _T.isLeaf(v);
  }
  
  /// Return anatomical site label of given leaf of underlying tree
  virtual const std::string& getLeafAnatomicalSiteLabel(Node v) const
  {
    assert(_T.isLeaf(v));
    return _T.l(v);
  }
  
  /// Return label of given vertex of underlying tree
  virtual const std::string& getLabel(Node v) const
  {
    return _T.label(v);
  }
  
  /// Return lca vertex of given leaf labels of underlying tree
  virtual Node getLCA(const NodeSet& vertices) const
  {
    return _T.getLCA(vertices);
  }
  
  /// Return parent of given vertex of underlying tree
  virtual Node getParent(Node v) const
  {
    return _T.parent(v);
  }
  
  /// Return whether vertex u is an ancestor of vertex v in the underlying tree
  virtual bool isAncestor(Node u, Node v) const
  {
    return _T.isAncestor(u, v);
  }

  /// Return anatomical sites
  virtual StringSet getAnatomicalSites() const
  {
    return _T.getAnatomicalSites();
  }
  
  /// Return original clone tree
  const CloneTree& getOrgT() const
  {
    return _T;
  }
  
  /// Construct graph
  void constructGraph();
  
  void initConstraintsLCA(int s,
                          Node v_k,
                          Node v_i,
                          Node v_j);
  
  /// Solve PMH under a topological constraint
  ///
  /// @param solver Solver
  /// @param T Non-binary clone tree
  /// @param primary Primary tumor
  /// @param outputDirectory Output directory
  /// @param outputPrefix Prefix prepended to every output filename
  /// @param colorMap Color map
  /// @param pattern Topological constraint
  /// @param nrThreads Number of threads
  /// @param outputILP Output ILP model
  /// @param outputSearchGraph Output search graph
  /// @param timeLimit Time limit in seconds
  /// @param bounds Upper bounds on mu, gamma and sigma
  /// @param forcedComigrations List of ordered pairs of anatomical sites
  /// that must be present
  static IntTriple run(IlpPmhSolver& solver,
                       const CloneTree& T,
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
                       const StringPairList& forcedComigrations);
  
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

private:
  /// Clone tree
  const CloneTree& _T;

protected:
  /// Label of primary tumor anatomical site
  const std::string& _primary;
  /// Topological constraint
  const MigrationGraph::Pattern _pattern;
  /// List of ordered pairs of anatomical sites that must be present
  const StringPairList& _forcedComigrations;
  
  /// lca(s)
  NodeVector _lca;
  /// Anatomical site label to index map
  StringToIntMap _anatomicalSiteToIndex;
  /// Index to anatomical site label map
  StringVector _indexToAnatomicalSite;
  /// Node to index
  NodeVector _indexToNode;
  /// Index to node
  IntNodeMap* _pNodeToIndex;
  /// Arc to index
  ArcVector _indexToArc;
  /// Index to arc
  IntArcMap* _pArcToIndex;
  /// Primary tumor index
  int _primaryIndex;
  /// Leaves per anatomical site index
  NodeSetVector _L;
  /// Maximum number of edges in G
  int _maxNrEdgesInG;

  /// Gurobi environment
  GRBEnv _env;
  /// Gurobi model
  GRBModel _model;
  /// x[i][s][c] = 1 iff vertex v_i occurs in X_(s,c)
  Var3Matrix _x;
  /// xx[ij][s][c][t][d] = 1 iff vertex (v_i, v_j) in E(T), x[i][s][c] = 1 and x[j][t][d] = 1
  Var5Matrix _xx;
  /// y[s][c] = 1 iff there exists a leaf labeled by (s,c)
  VarMatrix _y;
  /// z[s][c][t][d] = 1 if and only if there exists an edge from (s,c) to (t,d)
  Var4Matrix _z;
  /// w[s][t][d][e] = 1 if and only if anatomical sites (t,d) and (t,e)
  /// both have an incoming edge from an anatomical site $s$
  Var4Matrix _w;
  /// gamma[s][t] = 1 if and only if there exists c,d s.t. z[s][c][t][d] = 1
  VarMatrix _gamma;
  /// sigma[s] = 1 if and only if there exists t s.t. gamma[s][t] = 1
  VarArray _sigma;
  
  /// Lowerbound on the optimal solution (inferred by Gurobi)
  double _LB;
  /// Upperbound on the optimal solution (inferred by Gurobi)
  double _UB;
  
  /// Graph
  Digraph _G;
  /// Root node
  Node _rootG;
  /// Node to sub label
  IntPairNodeMap _nodeGToSubLabel;
  /// Sub label to node
  NodeMatrix _subLabelToNodeG;
  /// Resulting vertex labeling
  StringNodeMap* _pLPlus;
};

#endif // ILPPMHSOLVER_H
