/*
 * gabowmyers.h
 *
 *  Created on: 18-oct-2016
 *      Author: M. El-Kebir
 */

#ifndef GABOWMYERS_H
#define GABOWMYERS_H

#include "utils.h"
#include <lemon/bfs.h>

/// This class uses the Gabow-Myers algorithm to enumerate
/// all rooted spanning trees in a directed graph.
/// @brief This class implements the Gabow-Myers algorithm
class GabowMyers
{
public:
  DIGRAPH_TYPEDEFS(Digraph);

  /// Constructor
  ///
  /// @param G Directed graph
  /// @param root Root
  /// @param limit Maximum number of trees to enumerate
  GabowMyers(const Digraph& G, Node root, int limit = -1);
  
  /// Run enumeration
  void run();
  
  /// Print the directed graph in DOT format
  ///
  /// @param label Node labels
  /// @param out Output stream
  void writeDOT(const StringNodeMap& label, std::ostream& out) const;
  
  /// Return an enumerated spanning tree
  ///
  /// @param T Placeholder for the enumerated spanning tree
  /// @param i Index
  void result(SubDigraph& T,
              int i) const;
  
  /// Return the number of enumerated spanning trees
  int getNrTrees() const
  {
    return _result.size();
  }
  
protected:
  typedef lemon::Bfs<SubDigraph> SubBfs;
  typedef std::list<Arc> ArcList;
  typedef ArcList::const_iterator ArcListIt;
  typedef ArcList::iterator ArcListNonConstIt;
  typedef ArcList::const_reverse_iterator ArcListRevIt;
  typedef std::vector<ArcList> ArcListVector;
  typedef std::list<ArcList> ArcListList;
  typedef Digraph::NodeMap<IntSet> IntSetNodeMap;
  
  /// Print a subtree in DOT format
  ///
  /// @param T Subtree
  /// @param out Output stream
  void writeDOT(const SubDigraph& T, std::ostream& out) const;
  
private:
  /// Directed graph
  const Digraph& _G;
  /// Limit
  const int _limit;
  /// Root node
  Node _root;
  /// Resulting list of lists of arcs
  ArcListList _result;
  /// Number of vertices
  int _nrVertices;
  
  /// Extend the provided spanning tree using arcs present
  /// in the provided frontier and directed graph
  /// Returns false if limit reached
  ///
  /// @param G Directed graph
  /// @param T Tree to extend
  /// @param L Last reported tree
  /// @param bfsL BFS levels for L used for the bridge test
  /// @param F Frontier
  bool grow(SubDigraph& G,
            SubDigraph& T,
            SubDigraph& L,
            SubBfs& bfsL,
            ArcList& F);

  /// Initialize the enumeration algorithm.
  /// Add root vertex to T, initialize F with incident arc.
  ///
  /// @param G Directed graph
  /// @param T Tree to extend
  /// @param F Frontier
  void init(SubDigraph& G,
            SubDigraph& T,
            ArcList& F);
};

#endif // GABOWMYERS_H
