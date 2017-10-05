/*
 * ilppmhtisolver.h
 *
 *  Created on: 13-sep-2017
 *      Author: M. El-Kebir
 */

#ifndef ILPPMHTISOLVER_H
#define ILPPMHTISOLVER_H

#include "ilppmhtrsolver.h"
#include "frequencymatrix.h"

/// This class solves the Parsimonious Migration History with Tree Inference (PMH-TI) problem
class IlpPmhTiSolver : public IlpPmhTrSolver
{
public:
  /// Constructor
  ///
  /// @param T Mutation tree
  /// @param F Frequency matrices F- and F+
  /// @param primary Primary tumor
  /// @param pattern Topological constraint
  /// @param gurobiLogFilename Gurobi logging filename
  /// @param forcedComigrations List of ordered pairs of anatomical sites
  /// that must be present
  /// @param disablePolytomyResolution No polytomy resolution
  IlpPmhTiSolver(const CloneTree& T,
                 const FrequencyMatrix& F,
                 const std::string& primary,
                 MigrationGraph::Pattern pattern,
                 const std::string& gurobiLogFilename,
                 const StringPairList& forcedComigrations,
                 bool disablePolytomyResolution);
  
  /// Destructor
  virtual ~IlpPmhTiSolver();
  
  /// Write clone tree in DOT format
  ///
  /// @param out Output stream
  /// @param colorMap Color map
  virtual void writeCloneTree(std::ostream& out,
                              const StringToIntMap& colorMap) const;
  
  /// Solve PMH-TI under a topological constraint
  ///
  /// @param T Mutation tree
  /// @param F Frequency matrix
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
  /// @param disablePolytomyResolution No polytomy resolution
  static void run(const CloneTree& T,
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
                  bool disablePolytomyResolution);
  
  /// Write search graph
  void writeSearchGraphDOT(std::ostream& out) const;
  
  /// Return resulting proportion matrix
  const DoubleVectorNodeMap& getU() const
  {
    return *_pResU;
  }
  
  /// Return resulting frequency matrix
  const DoubleVectorNodeMap& getF() const
  {
    return *_pResF;
  }
  
protected:
  /// Initialize indices and mappings
  virtual void initIndices();
  
  /// Initialize ILP variables
  virtual void initVariables();
  
  /// Initialize ILP leaf variables
  virtual void initLeafVariables();
  
  /// Initialize ILP constraints
  virtual void initConstraints();
  
  /// Initialize ILP leaf constraints
  virtual void initLeafConstraints();
  
  /// Initialize ILP S constraints
  virtual void initSingleSourceSeedingConstraints();
  
  /// Initialize ILP objective function
  virtual GRBLinExpr initObjective(const IntTriple& bounds);
  
  /// Return underlying tree
  virtual const Digraph& getTree() const
  {
    return _extT;
  }
  
  /// Return root vertex of underlying tree
  virtual Node getRoot() const
  {
    return _rootExtT;
  }
  
  /// Return whether given node of underlying tree is a leaf
  virtual bool isLeaf(Node v) const
  {
    return OutArcIt(_extT, v) == lemon::INVALID;
  }
  
  /// Return anatomical site label of given leaf of underlying tree
  virtual const std::string& getLeafAnatomicalSiteLabel(Node v) const
  {
    assert(isLeaf(v));
    return _leafAnatomicalSiteLabelExtT[v];
  }
  
  /// Return label of given vertex of underlying tree
  virtual const std::string& getLabel(Node v) const
  {
    return _labelExtT[v];
  }
  
  /// Return lca vertex of given leaves of underlying tree
  virtual Node getLCA(const NodeSet& vertices) const
  {
    return BaseTree::getLCA(_extT, vertices);
  }
  
  /// Return parent of given vertex of underlying tree
  virtual Node getParent(Node v) const
  {
    assert(v != _rootExtT);
    return _extT.source(InArcIt(_extT, v));
  }
  
  /// Return whether vertex u is an ancestor of vertex v in the underlying tree
  virtual bool isAncestor(Node u, Node v) const
  {
    while (v != _rootExtT)
    {
      if (u == v)
        return true;
      else
        v = _extT.source(InArcIt(_extT, v));
    }
    
    return u == v;
  }
  
  /// Return anatomical sites
  virtual StringSet getAnatomicalSites() const
  {
    return _F.getAnatomicalSites();
  }
  
  /// Process ILP solution
  virtual void processSolution();
  
private:
  /// Compute Fmin
  ///
  /// @param v_i Node in mutation tree T
  void computeFmin(Node v_i);
  
  /// Compute Fmax
  ///
  /// @param v_i Node in mutation tree T
  void computeFmax(Node v_i);
  
  /// Compute Umax
  void computeUmax();
  
  /// Construct search graph G recursively
  ///
  /// @param v Node in mutation tree
  void constructExtT(Node v);
  
private:
  /// Frequency matrix
  const FrequencyMatrix& _F;
  /// No polytomy resolution
  bool _disablePolytomyResolution;
  /// Inferred frequency lowerbounds given mutation tree
  DoubleMatrix _Fmin;
  /// Inferred frequency upperbounds given mutation tree
  DoubleMatrix _Fmax;
  /// Inferred mixing proportion upper bounds given mutation tree
  DoubleMatrix _Umax;
  
  /// Extended mutation tree
  Digraph _extT;
  /// Root of extended mutation tree
  Node _rootExtT;
  StringNodeMap _labelExtT;
  /// Anatomical site label of leaf of extT
  StringNodeMap _leafAnatomicalSiteLabelExtT;
  /// Node mapping from V(T) to V(extT)
  NodeNodeMap _TtoExtT;
  /// Node mapping from V(extT) to V(T)
  NodeNodeMap _extTtoT;
  /// Denotes whether node is an anatomical site attachment node
  BoolNodeMap _isAnatomicalSiteNodeExtT;
  
  /// Variable f[p][i] is frequency of mutation incoming to clone i in anatomical site s
  VarMatrix _f;
  /// Variable u[p][i] is mixing proportion of clone i in anatomical site s
  VarMatrix _u;
  
  /// Resulting frequencies in resulting clone tree
  DoubleVectorNodeMap* _pResF;
  /// Resulting mixing proportions of leaves of resulting clone
  DoubleVectorNodeMap* _pResU;
};

#endif // ILPPMHTISOLVER_H
