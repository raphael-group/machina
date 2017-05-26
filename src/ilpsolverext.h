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
#include "ilpsolver.h"
#include <gurobi_c++.h>

/// This class implements an ILP for solving PMH-CTI under topological constraints
class IlpSolverExt : public IlpSolver
{
public:
  /// Constructor
  ///
  /// @param T Non-binary clone tree
  /// @param F Frequency matrices F- and F+
  /// @param primary Primary tumor
  /// @param pattern Topological constraint
  /// @param gurobiLogFilename Gurobi logging filename
  /// @param forcedComigrations List of ordered pairs of anatomical sites
  /// that must be present
  IlpSolverExt(const NonBinaryCloneTree& T,
               const FrequencyMatrix& F,
               const std::string& primary,
               MigrationGraph::Pattern pattern,
               const std::string& gurobiLogFilename,
               const StringPairList& forcedComigrations);
  
  /// Destructor
  virtual ~IlpSolverExt();
  
  /// Return vertex labeling
  virtual const StringNodeMap& lPlus() const
  {
    return *_pResLPlus;
  }
  
  /// Return resulting clone tree
  virtual const NonBinaryCloneTree& T() const
  {
    return *_pResCloneTree;
  }
  
  /// Initialize ILP
  virtual void init(double upperBound);
  
  /// Return resulting frequencies
  const DoubleVectorNodeMap& getF() const
  {
    return *_pResF;
  }
  
  /// Return resulting mixing proportions
  const DoubleNodeMap& getU() const
  {
    return *_pResU;
  }
  
  /// Return resulting character labeling
  const IntNodeMap& getCharacterLabel() const
  {
    return *_pResCharacterLabel;
  }
  
  /// Print search graph G
  ///
  /// @param out Output stream
  void writeDOT(std::ostream& out) const;
  
  /// Print search graph G
  ///
  /// @param out Output stream
  /// @param colorMap Color map
  void writeDOT(std::ostream& out,
                const StringToIntMap& colorMap) const;
  
  /// Solve PMH-CTI under a topological constraint
  ///
  /// @param T Non-binary mutation tree
  /// @param F Frequency matrices F- and F+
  /// @param primary Primary tumor
  /// @param outputDirectory Output directory
  /// @param colorMap Color map
  /// @param pattern Topological constraint
  /// @param nrThreads Number of threads
  /// @param outputILP Output ILP model
  /// @param outputSearchGraph Output search graph
  /// @param timeLimit Time limit in seconds
  ///
  /// @param UB Upper bound on objective value
  /// @param forcedComigrations List of ordered pairs of anatomical sites
  /// that must be present
  static void run(const NonBinaryCloneTree& T,
                  const FrequencyMatrix& F,
                  const std::string& primary,
                  const std::string& outputDirectory,
                  const StringToIntMap& colorMap,
                  MigrationGraph::Pattern pattern,
                  int nrThreads,
                  bool outputILP,
                  bool outputSearchGraph,
                  int timeLimit,
                  double UB,
                  const StringPairList& forcedComigrations);
  
protected:
  /// Construct search graph G
  virtual void constructG();
  
  /// Construct search graph G recursively
  ///
  /// @param v Node
  void constructG(Node v);
  
  /// Initialize indices and mappings
  virtual void initIndices();
  
  /// Initialize ILP variables
  virtual void initVariables();
  
  /// Initialize ILP constraints
  virtual void initLeafConstraints();
  
  /// Process ILP solution
  virtual void initConstraints();
  
  /// Initialize ILP objective
  ///
  /// @param upperBound Upper bound on objective function
  virtual void initObjective(double upperBound);
  
  /// Process ILP solution
  virtual void processSolution();

  /// Return the number of mutation tree vertices
  virtual int getNrMutationTreeVertices() const
  {
    return _F.getNrCharacters();
  }
  
  /// Add ILP contraint involving matching anatomical sites of adjacent vertices
  ///
  /// @param sum_z Sum of z variables
  /// @param ij Index of arc (v_i, v_j)
  virtual void addMatchingColorsConstraint(GRBLinExpr sum_z,
                                           int ij)
  {
    _model.addConstr(sum_z + _y[ij] == _w[ij]);
  }
  
  /// Return underlying LEMON search graph
  virtual const Digraph& tree() const
  {
    return _G;
  }
  
  /// Add ILP usage constraint
  ///
  /// @param ij Arc (v_i, v_j) index
  /// @param mapped_i Character index
  /// @param t Sample index
  virtual void addMatchingColorsAndUsageConstraint(int ij,
                                                   int mapped_i,
                                                   int t)
  {
    const int nrSamples = _indexToSample.size();
    _model.addConstr(_z[ij][nrSamples] <= 1 - _u[t][mapped_i]);
  }
  
  /// Add ILP comigration constraint
  ///
  /// @param s Source anatomical site
  /// @param t Target anatomical site
  /// @param ij Arc (v_i, v_j) index
  /// @param i Node v_i index
  /// @param j Node v_j index
  virtual void addComigrationConstraint(int s, int t,
                                        int ij, int i, int j)
  {
    _model.addConstr(_c[s][t] - _x[i][s] - _x[j][t] >= -1 - (1-_w[ij]));
  }
  
  /// Return leaf anatomical site label
  ///
  /// @param v Node in G
  virtual const std::string& l(Node v) const
  {
    assert(isLeaf(v));
    return _l[v];
  }
  
  /// Return node identifier
  ///
  /// @param v Node in G
  virtual const std::string& label(Node v) const
  {
    return _label[v];
  }
  
  /// Return root node
  virtual Node root() const
  {
    return _root;
  }
  
  /// Return whether given node is a leaf
  ///
  /// @param v Node in G
  virtual bool isLeaf(Node v) const
  {
    return OutArcIt(_G, v) == lemon::INVALID;
  }
  
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
  
  /// Label edges by introduced characters recursively
  ///
  /// @param resT LEMON clone tree
  /// @param resNodeToIndex Node to index mapping
  /// @param v_i Current node
  /// @param characterLabel Output node labeling by introduced characters
  /// @param visited Marks characters as introduced
  void labelEdges(const Digraph& resT,
                  const IntNodeMap& resNodeToIndex,
                  Node v_i,
                  IntNodeMap& characterLabel,
                  BoolVector& visited);
  
  /// Return introduced character index for given node
  ///
  /// @param resT LEMON clone tree
  /// @param resNodeToIndex Node to index mapping
  /// @param v_i Node
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
  
protected:
  /// Frequency matrices F- and F+
  const FrequencyMatrix& _F;
  /// Inferred frequency lowerbounds given mutation tree
  DoubleMatrix _Fmin;
  /// Inferred frequency upperbounds given mutation tree
  DoubleMatrix _Fmax;
  /// Inferred mixing proportion upper bounds given mutation tree
  DoubleMatrix _Umax;
  
  /// Search graph
  Digraph _G;
  /// Root node of _G
  Node _root;
  /// Node identifier
  StringNodeMap _label;
  /// Leaf labeling
  StringNodeMap _l;
  /// Vertex labeling
  StringNodeMap _lPlus;
  /// Node mapping from V(T) to V(G)
  NodeNodeMap _TtoG;
  /// Node mapping from V(G) to V(T)
  NodeNodeMap _GtoT;
  /// Denotes whether a node is a sample attachment node
  BoolNodeMap _isSampleNode;
  
  /// Variable f[s][i] is frequency of mutation incoming to clone i in sample s
  VarMatrix _f;
  /// Variable u[s][i] is mixing proportion of clone i in sample s
  VarMatrix _u;
  /// Variable w[ij] = 1 iff edge (v_i, v_j) is in the tree
  VarArray _w;
  
  /// Resulting frequencies
  DoubleVectorNodeMap _resF;
  /// Resulting mixing proportions
  DoubleVectorNodeMap _resU;
  /// Resulting clone tree
  NonBinaryCloneTree* _pResCloneTree;
  /// Resulting vertex labeling
  StringNodeMap* _pResLPlus;
  /// Resulting frequencies in resulting clone tree
  DoubleVectorNodeMap* _pResF;
  /// Resulting mixing proportions of leaves of resulting clone
  DoubleNodeMap* _pResU;
  /// Resulting character labeling of vertices of resulting clone tree
  IntNodeMap* _pResCharacterLabel;
};

#endif // ILPSOLVEREXT_H
