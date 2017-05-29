/*
 * clonetree.h
 *
 *  Created on: 19-oct-2016
 *      Author: M. El-Kebir
 */

#ifndef CLONETREE_H
#define CLONETREE_H

#include "utils.h"
#include "binarytree.h"
#include "nonbinaryclonetree.h"

/// This class models a full binary clone tree
class CloneTree : public BinaryTree
{
public:
  /// Default constructor
  CloneTree();
  
  /// Constructor
  ///
  /// @param T Directed graph
  /// @param root Root node
  /// @param id Node identifier
  /// @param l Leaf label
  CloneTree(const Digraph& T,
            Node root,
            const StringNodeMap& id,
            const StringNodeMap& l);
  
  /// Copy constructor
  CloneTree(const CloneTree& other);
  
  /// Initialize a clone tree given a Dyck word and randomly assigns anatomical
  /// sites to its leaves.
  ///
  /// @param dyckWord Dyck word
  /// @param nrSamples Number of anatomical sites
  CloneTree(const BoolVector& dyckWord, int nrSamples);
  
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
  
  /// Return the leaf label of the given node
  ///
  /// @param u Node
  const std::string& l(Node u) const
  {
    assert(_isLeaf[u]);
    return _l[u];
  }
  
  /// Read leaf labeling
  ///
  /// @param in Input stream
  bool readLeafLabeling(std::istream& in);
  
  /// Write leaf labeling
  ///
  /// @param out Output stream
  void writeLeafLabeling(std::ostream& out) const;
  
  /// Write vertex labeling
  ///
  /// @param out Output stream
  /// @param lPlus Vertex labeling
  void writeVertexLabeling(std::ostream& out,
                           const StringNodeMap& lPlus) const;
  
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
  
protected:
  /// Leaf labeling L(T) -> Sigma
  StringNodeMap _l;
};

#endif // CLONETREE_H
