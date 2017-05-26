/*
 * gabowmyers.h
 *
 *  Created on: 18-oct-2016
 *      Author: M. El-Kebir
 */

#ifndef GABOWMYERS_H
#define GABOWMYERS_H

#include "utils.h"

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
  GabowMyers(const Digraph& G, Node root);
  
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
  /// Root node
  Node _root;
  /// Resulting list of lists of arcs
  ArcListList _result;
  
  /// Extend the provided spanning tree using arcs present
  /// in the provided frontier and directed graph
  ///
  /// @param G Directed graph
  /// @param T Tree to extend
  /// @param F Frontier
  void grow(SubDigraph& G,
            SubDigraph& T,
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
  
  /// Add the given arc to the subtree
  ///
  /// @param T Subtree
  /// @param a_cidj Arc to add
  void addArc(SubDigraph& T,
              Arc a_cidj) const;
  
  /// Remove the given arc from the subtree
  ///
  /// @param T Subtree
  /// @param a_cidj Arc to remove
  void removeArc(SubDigraph& T,
                 Arc a_cidj) const;
  
  /// Decide whether the given subtree is a directed tree
  bool isArborescence(const SubDigraph& T) const;
  
  /// Report the enumerated spanning tree
  void finalize(const SubDigraph& T);
};

#endif // GABOWMYERS_H
