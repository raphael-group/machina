/*
 * binarytree.h
 *
 *  Created on: 19-oct-2016
 *      Author: M. El-Kebir
 */

#ifndef BINARYTREE_H
#define BINARYTREE_H

#include "utils.h"
#include "basetree.h"

/// This class models a full binary tree
class BinaryTree : public BaseTree
{
public:
  /// Default constructor
  BinaryTree();
  
  /// Constructor
  ///
  /// @param T Directed graph
  /// @param root Root node
  /// @param label Node identifier
  BinaryTree(const Digraph& T,
             Node root,
             const StringNodeMap& label);
  
  /// Copy constructor
  BinaryTree(const BinaryTree& other);
  
  /// Constructs a full binary tree given a Dyck work
  ///
  /// @param dyckWord Dyck word
  BinaryTree(const BoolVector& dyckWord);
  
  /// Returns left child of a node
  ///
  /// @param u Node
  Node left(Node u) const
  {
    if (_isLeaf[u])
    {
      return lemon::INVALID;
    }
    else
    {
      return _tree.target(OutArcIt(_tree, u));
    }
  }
  
  /// Returns right child of a node
  ///
  /// @param u Node
  Node right(Node u) const
  {
    if (_isLeaf[u])
    {
      return lemon::INVALID;
    }
    else
    {
      return _tree.target(++OutArcIt(_tree, u));
    }
  }
  
  /// Returns the sibling of a node. In case of the root, lemon::INVALID is returned.
  ///
  /// @param v Node
  Node sibling(Node v) const
  {
    if (v == _root)
    {
      return lemon::INVALID;
    }
    else
    {
      Node u = parent(v);
      if (v == left(u))
      {
        return right(u);
      }
      else
      {
        return left(u);
      }
    }
  }
  
  /// Returns the corresponding Dyck word
  BoolVector dyckWord() const
  {
    return dyckWord(_root);
  }
  
  /// Returns a string that displays the provided Dyck word in brack notation
  ///
  /// @param w Dyck word
  static std::string bracketNotation(const BoolVector& w);
  
  /// Generates a random Dyck word
  ///
  /// @param n Size
  static BoolVector randomDyckWord(int n);
  
protected:
  /// Decides whether the tree is valid
  virtual bool isValid() const
  {
    if (!BaseTree::isValid())
      return false;
    
    for (NodeIt u(_tree); u != lemon::INVALID; ++u)
    {
      if (isLeaf(u))
        continue;
      
      int count = 0;
      for (OutArcIt a(_tree, u); a != lemon::INVALID; ++a, ++count);
      if (count != 2)
        return false;
    }
    return true;
  }
  
  /// Init tree topology given a Dyck word
  Node initTree(const BoolVector& dyckWord);
  
  BoolVector dyckWord(Node v) const;
  
  static BoolVector map(const BoolVector& w);
  
  static bool isWellFormed(const BoolVector& w);
  
  static bool isBalanced(const BoolVector& w);
  
  static int match(const BoolVector& w, int i);
};

#endif // BINARYTREE_H
