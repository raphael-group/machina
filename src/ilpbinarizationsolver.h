/*
 * ilpbinarizationsolver.h
 *
 *  Created on: 3-feb-2016
 *      Author: M. El-Kebir
 */

#ifndef ILPBINARIZATIONSOLVER_H
#define ILPBINARIZATIONSOLVER_H

#include "ilpsolver.h"
#include "clonetree.h"

/// This class implements an ILP for solving the Parsimonious Migration History
/// with Polytomy Resolution (PMH-PR) problem under various topological
/// constraints (PS, S, M and R)
///
/// \brief This class implements an ILP for solving PMH-PR under topological
/// constraints.
class IlpBinarizationSolver : public IlpSolver
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
  IlpBinarizationSolver(const NonBinaryCloneTree& T,
                        const std::string& primary,
                        MigrationGraph::Pattern pattern,
                        const std::string& gurobiLogFilename,
                        const StringPairList& forcedComigrations);
  
  /// Destructor
  virtual ~IlpBinarizationSolver()
  {
    delete _pResLPlus;
    delete _pResCloneTree;
  }
  
  /// Return vertex labeling
  virtual const StringNodeMap& lPlus() const
  {
    return *_pResLPlus;
  }
  
  /// Write DOT file using given color map
  ///
  /// @param out Output stream
  /// @param colorMap Color map
  void writeDOT(std::ostream& out,
                const StringToIntMap& colorMap) const;
  
  /// Return search graph
  virtual const Digraph& tree() const
  {
    return _G;
  }
  
  /// Return resulting clone tree
  const CloneTree& getCloneTree() const
  {
    return *_pResCloneTree;
  }
  
protected:
  /// Initialize indices and mappings
  virtual void initIndices();
  
  /// Initialize ILP variables
  virtual void initVariables();
  
  /// Initialize ILP constraints
  virtual void initConstraints();
  
  /// Process ILP solution
  virtual void processSolution();
  
  /// Add ILP contraint involving matching anatomical sites of adjacent vertices
  ///
  /// @param sum_z Sum of z variables
  /// @param ij Index of arc (v_i, v_j)
  virtual void addMatchingColorsConstraint(GRBLinExpr sum_z, int ij)
  {
    _model.addConstr(sum_z + _y[ij] == _w[ij]);
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
  
  /// Construct search graph G recursively
  ///
  /// @param T Clone tree
  /// @param v Node
  void constructG(const Digraph& T, Node v);

  /// Return node identifier
  virtual const std::string& label(Node v) const
  {
    return _label[v];
  }
  
  /// Return leaf anatomical site label
  ///
  /// @param v Node in G
  virtual const std::string& l(Node v) const
  {
    Node org_v = _GtoT[v];
    assert(_T.isLeaf(org_v));
    return _T.l(org_v);
  }
  
  /// Return root node
  virtual Node root() const
  {
    return _TtoG[_T.root()];
  }
  
  /// Return whether given node is a leaf
  ///
  /// @param v Node in G
  virtual bool isLeaf(Node v) const
  {
    return _T.isLeaf(_GtoT[v]);
  }
  
private:
  /// Search graph whose full binary spanning trees are valid binarization
  Digraph _G;
  /// Node identifier
  StringNodeMap _label;
  /// Node mapping from V(T) to V(G)
  NodeNodeMap _TtoG;
  /// Node mapping from V(G) to V(T)
  NodeNodeMap _GtoT;
  /// w[ij] = 1 iff edge (v_i, v_j) is in the tree
  VarArray _w;
  /// Resulting clone tree
  CloneTree* _pResCloneTree;
  /// Resulting vertex labeling
  StringNodeMap* _pResLPlus;
};

#endif // ILPBINARIZATIONSOLVER_H
