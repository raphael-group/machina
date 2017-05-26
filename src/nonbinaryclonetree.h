/*
 * nonbinaryclonetree.h
 *
 *  Created on: 19-oct-2016
 *      Author: M. El-Kebir
 */

#ifndef NONBINARYCLONETREE_H
#define NONBINARYCLONETREE_H

#include "utils.h"
#include "binarytree.h"
#include "clonetree.h"

/// This class models a clone tree
class NonBinaryCloneTree : public BaseTree
{
public:
  /// Default constructor
  NonBinaryCloneTree();
  
  /// Copy constructor
  NonBinaryCloneTree(const NonBinaryCloneTree& other);
  
  /// Constructor
  ///
  /// @param T Directed graph
  /// @param root Root node
  /// @param id Node identifier
  /// @param l Leaf label
  NonBinaryCloneTree(const Digraph& T,
                     Node root,
                     const StringNodeMap& id,
                     const StringNodeMap& l);
  
  /// Return the leaf label of the given node
  ///
  /// @param u Node
  const std::string& l(Node u) const
  {
    assert(_isLeaf[u]);
    return _l[u];
  }
  
  /// Return the set of leaf labels of subtree rooted at the given node
  ///
  /// @param u Node
  StringSet ll(Node u) const
  {
    StringSet res;
    for (Node v : _leafSubset[u])
    {
      res.insert(l(v));
    }
    return res;
  }
  
  /// Read leaf labeling
  ///
  /// @param in Input stream
  bool readLeafLabeling(std::istream& in);
  
  /// Write leaf labeling
  ///
  /// @param out Output stream
  void writeLeafLabeling(std::ostream& out) const;
  
  /// Print tree in DOT format
  ///
  /// @param out Output stream
  void writeDOT(std::ostream& out) const;
  
  /// Print tree in DOT format using the given color map
  ///
  /// @param out Output stream
  /// @param colorMap Color map
  void writeDOT(std::ostream& out,
                const StringToIntMap& colorMap) const;
  
  /// Print tree in DOT format using the given vertex labeling
  ///
  /// @param out Output stream
  /// @param lPlus Labeling of each node by an anatomical site
  void writeDOT(std::ostream& out,
                const StringNodeMap& lPlus) const;
  
  /// Print tree in DOT format using the given vertex labeling and color map
  ///
  /// @param out Output stream
  /// @param lPlus Labeling of each node by an anatomical site
  /// @param colorMap Color map
  void writeDOT(std::ostream& out,
                const StringNodeMap& lPlus,
                const StringToIntMap& colorMap) const;
  
  /// Print tree in DOT format using the given vertex labeling and color map
  ///
  /// @param out Output stream
  /// @param lPlus Labeling of each node by an anatomical site
  /// @param colorMap Color map
  /// @param F Frequencies
  /// @param U Mixing proportions of the leaves
  void writeDOT(std::ostream& out,
                const StringNodeMap& lPlus,
                const StringToIntMap& colorMap,
                const DoubleVectorNodeMap& F,
                const DoubleNodeMap& U) const;

  /// Print tree in DOT format using the given vertex labeling and color map
  ///
  /// @param out Output stream
  /// @param lPlus Labeling of each node by an anatomical site
  /// @param colorMap Color map
  /// @param U Mixing proportions of the leaves
  /// @param characterLabel Character index of the character
  ///        introduced on the incoming edge of a node
  void writeDOT(std::ostream& out,
                const StringNodeMap& lPlus,
                const StringToIntMap& colorMap,
                const DoubleNodeMap& U,
                const IntNodeMap& characterLabel) const;
  
  /// Return number of anatomical sites
  int getNrSamples() const
  {
    return getSamples().size();
  }
  
  /// Return anatomical site labels
  StringSet getSamples() const
  {
    StringSet res;
    
    for (NodeIt v(_tree); v != lemon::INVALID; ++v)
    {
      if (isLeaf(v))
      {
        res.insert(_l[v]);
      }
    }
    
    return res;
  }
  
  /// Return all migration edges of the provided vertex labeling
  ///
  /// @param lPlus Labeling of each node by an anatomical site
  ArcSet getMigrationEdges(const StringNodeMap& lPlus) const
  {
    ArcSet res;
    
    for (ArcIt a(_tree); a != lemon::INVALID; ++a)
    {
      Node u = _tree.source(a);
      Node v = _tree.target(a);
      if (lPlus[u] != lPlus[v])
      {
        res.insert(a);
      }
    }
    
    return res;
  }
  
  /// Merge sibling leaves labeled by the same anatomical site
  void mergeSameSampleSiblingLeaves();
  
protected:
  /// Leaf labeling L(T) -> Sigma
  StringNodeMap _l;
};

#endif // NONBINARYCLONETREE_H
