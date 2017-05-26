/*
 * basetree.h
 *
 *  Created on: 13-oct-2016
 *      Author: M. El-Kebir
 */

#ifndef BASETREE_H
#define BASETREE_H

#include "utils.h"

/// This class models a tree whose nodes are labeled by unique identifiers
class BaseTree
{
public:
  /// Default constructor
  BaseTree();
  
  /// Constructor
  ///
  /// @param T Directed graph
  /// @param root Root node
  /// @param id Node identifier
  BaseTree(const Digraph& T,
           Node root,
           const StringNodeMap& id);
  
  /// Copy constructor
  BaseTree(const BaseTree& other);
  
  /// Destructor
  virtual ~BaseTree()
  {
  }
  
  /// Deserialize
  ///
  /// @param in Input stream
  virtual bool read(std::istream& in);
  
  /// Serialize
  ///
  /// @param out Output stream
  void write(std::ostream& out) const;
  
  /// Print tree in DOT format
  ///
  /// @param out Output stream
  virtual void writeDOT(std::ostream& out) const;
  
  /// Return number of anatomical sites
  virtual int getNrSamples() const
  {
    return 0;
  }
  
  /// Return anatomical site labels
  virtual StringSet getSamples() const
  {
    return StringSet();
  }
  
  /// Return underlying LEMON tree
  const Digraph& tree() const
  {
    return _tree;
  }
  
  /// Decide whether two nodes occur in distinct branches
  ///
  /// @param u Node
  /// @param v Node
  bool areIncomparable(Node u, Node v) const
  {
    return !isAncestor(u, v) && !isAncestor(v, u);
  }
  
  /// Decide whether one node is an ancestor (or self) of another
  ///
  /// @param u Node
  /// @param v Node
  bool isAncestor(Node u, Node v) const
  {
    while (v != _root)
    {
      if (u == v)
        return true;
      else
        v = _tree.source(InArcIt(_tree, v));
    }
    
    return u == v;
  }

  /// Return parent node. If u is the root lemon::INVALID is returned.
  ///
  /// @param u Node
  Node parent(Node u) const
  {
    if (u == _root)
    {
      return lemon::INVALID;
    }
    else
    {
      return _tree.source(InArcIt(_tree, u));
    }
  }
  
  /// Return the root node
  Node root() const
  {
    return _root;
  }
  
  /// Return the node identifier map
  const StringNodeMap& getLabelMap() const
  {
    return _nodeToId;
  }
  
  /// Return the identifier of a node
  ///
  /// @param u Node
  const std::string& label(Node u) const
  {
    return _nodeToId[u];
  }
  
  /// Decide whether a node is a leaf
  ///
  /// @param u Node
  bool isLeaf(Node u) const
  {
    return _isLeaf[u];
  }
  
  /// Return the BFS level of a node. The root node has BFS level 0.
  ///
  /// @param u Node
  int level(Node u) const
  {
    return _level[u];
  }
  
  /// Return the set of leaves that occur in the subtree rooted at the provided node.
  ///
  /// @param u Node
  const NodeSet& leafSubset(Node u) const
  {
    return _leafSubset[u];
  }
  
  /// Return the leaf set
  const NodeSet& leafSet() const
  {
    return _leafSet;
  }
  
  /// Return the LCA of the provided set of nodes
  ///
  /// @param nodes Node set
  Node getLCA(const NodeSet& nodes) const;
  
  /// Decides whether the given set of nodes is connected
  ///
  /// @param nodes Node set
  bool isConnected(const NodeSet& nodes) const;
  
  /// Initialize auxilliary data structures (_isLeaf, _leafSet, _leafSubset, _level)
  void init();
  
  /// Return a node corresponding to the given identifier. If there is no node with the given identifier, lemon::INVALID is returned.
  ///
  /// @param lbl Identifier
  Node getNodeByLabel(const std::string& lbl) const
  {
    auto it = _idToNode.find(lbl);
    if (it == _idToNode.end())
    {
      return lemon::INVALID;
    }
    else
    {
      return it->second;
    }
  }
  
  /// Return the unique path from the root to u
  NodeList pathFromRoot(Node u) const;
  
  /// Return the unique path from u to v
  NodeList path(Node u, Node v) const;
  
  /// Read vertex labeling
  ///
  /// @param in Input stream
  /// @param T Tree
  /// @param label Output labeling
  static bool readVertexLabeling(std::istream& in,
                                 const BaseTree& T,
                                 StringNodeMap& label);
  
  /// Read color map
  ///
  /// @param in Input stream
  /// @param colorMap Output color map
  static bool readColorMap(std::istream& in,
                           StringToIntMap& colorMap);
  
  /// Return a generated color map based on the anatomical sites present in the current tree
  StringToIntMap generateColorMap() const
  {
    StringSet samples = getSamples();
    
    StringToIntMap colorMap;
    int idx = 0;
    for (const std::string& s : samples)
    {
      colorMap[s] = ++idx;
    }
    
    return colorMap;
  }
  
  /// Return the edge list
  StringPairList getEdgeList() const
  {
    StringPairList res;
    for (ArcIt a(_tree); a != lemon::INVALID; ++a)
    {
      Node u = _tree.source(a);
      Node v = _tree.target(a);
      
      res.push_back(StringPair(_nodeToId[u], _nodeToId[v]));
    }
    return res;
  }

protected:
  /// Tree
  Digraph _tree;
  /// Root
  Node _root;
  /// Node label
  StringNodeMap _nodeToId;
  /// Get node by label
  StringToNodeMap _idToNode;
  /// Indicates whether a node is a leaf
  BoolNodeMap _isLeaf;
  /// Leaf set
  NodeSet _leafSet;
  /// Leaf set of subtree rooted at a node
  NodeNodeSetMap _leafSubset;
  /// BFS level of a node
  IntNodeMap _level;
  
  /// Initialize _leafSubset
  void initLeafSubset(Node node);
  
  /// Decide whether the tree is valid
  virtual bool isValid() const;
};

#endif // BASETREE_H
