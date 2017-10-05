/*
 * ilppmhtrsolver.h
 *
 *  Created on: 20-aug-2017
 *      Author: M. El-Kebir
 */

#ifndef ILPPMHTRSOLVER_H
#define ILPPMHTRSOLVER_H

#include "utils.h"
#include "clonetree.h"
#include <gurobi_c++.h>
#include "migrationgraph.h"
#include "solutiongraph.h"
#include "ilppmhsolver.h"
#include "ilppmhtrsolvercallback.h"

/// This class solves the Parsimonious Migration History with Tree Resolution Problem
class IlpPmhTrSolver : public IlpPmhSolver
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
  IlpPmhTrSolver(const CloneTree& T,
                 const std::string& primary,
                 MigrationGraph::Pattern pattern,
                 const std::string& gurobiLogFilename,
                 const StringPairList& forcedComigrations);
  
  virtual ~IlpPmhTrSolver();
  
  /// Solve PMH-TR under a topological constraint
  ///
  /// @param T Clone tree
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
  static void run(const CloneTree& T,
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
  
  /// Return refined clone tree of provided solution index
  virtual const CloneTree& T() const
  {
    return *_pTprime;
  }
  
protected:
  typedef std::map<std::string, std::string> StringToStringMap;
  
  /// Initialize ILP variables
  virtual void initVariables();
  
  /// Initialize ILP leaf variables
  virtual void initLeafVariables();
  
  /// Initialize ILP constraints regarding vertex labelings
  virtual void initVertexLabelingConstraints();
  
  /// Initialize ILP leaf constraints
  virtual void initLeafConstraints();
  
  /// Initialize ILP callbacks
  virtual void initCallbacks();
  
  /// Initialize ILP constraints
  virtual void initConstraints();
  
  /// Initialize ILP constraints regarding G
  virtual void initConstraintsG();
  
  /// Initialize ILP constraints regarding non-edges of G
  virtual void initConstraintsNonEdgesG();
  
  /// Intialize ILP M constraints
  virtual void initMultiSourceSeedingConstraints();
  
  /// Initialize ILP objective function
  virtual GRBLinExpr initObjective(const IntTriple& bounds);
  
  /// Process ILP solution
  virtual void processSolution();
  
  /// Refine clone tree according to identified _G
  void refine(const BoolNodeMap& leafPresence,
              StringToStringMap& toMutLabel);
  
  /// Refine clone tree
  void refine(const BoolNodeMap& leafPresence,
              StringToStringMap& toMutLabel,
              Node v_inT,
              Digraph& Tprime,
              Node v_inTprime,
              StringNodeMap& label,
              StringNodeMap& lPlus);
  
  /// Check whether Sigma_u induces a connected subgraph of G
  ///
  /// @param Sigma_u Set of states
  /// @param root_sc Root of the connected subgraph
  bool isConnected(const IntPairSet& Sigma_u,
                   IntPair& root_sc) const;
  
  /// Check whether subgraph of G induced by Sigma_u are siblings
  ///
  /// @param Sigma_u Set of states
  /// @param parent_sc Parent of siblings
  bool areSiblings(const IntPairSet& Sigma_u,
                   IntPair& parent_sc) const;
  
  /// Compute next combination of anatomical sites
  ///
  /// @param subset Current set of anatomical sites
  bool nextCombinationAnatomicalSites(BoolVector& subset) const;
  
  typedef std::vector<IntPair> IntPairVector;
  
  /// Compute next combination of states
  ///
  /// @param states Current combination of states
  bool nextCombinationStates(IntPairVector& states) const;
  
protected:
  /// Vertex of T to set of states in color
  IntPairSetNodeMap* _pNodeToStateSet;
  /// Vertex of T to root color
  IntPairNodeMap* _pNodeToRootState;
  /// Refined clone tree
  CloneTree* _pTprime;
  /// _z[i][s][c][t][d]
  Var5Matrix _zz;
  /// _r[i][s][c], root color
  Var3Matrix _r;
  /// Callback class
  IlpPmhPrSolverCycleElimination* _pCallback;
};

#endif // ILPPMHTRSOLVER_H
