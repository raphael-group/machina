/*
 * binarytree.cpp
 *
 *  Created on: 19-oct-2016
 *      Author: M. El-Kebir
 */

#include "binarytree.h"
#include <lemon/bfs.h>

BinaryTree::BinaryTree()
  : BaseTree()
{
}

BinaryTree::BinaryTree(const Digraph& T,
                       Node root,
                       const StringNodeMap& label)
  : BaseTree(T, root, label)
{
}

BinaryTree::BinaryTree(const BinaryTree& other)
  : BaseTree(other)
{
}

BinaryTree::BinaryTree(const BoolVector& dyckWord)
  : BaseTree()
{
  _root = initTree(dyckWord);
  init();
  for (NodeIt v(_tree); v != lemon::INVALID; ++v)
  {
    std::string str = boost::lexical_cast<std::string>(_tree.id(v));
    _nodeToId[v] = str;
    _idToNode[str] = v;
  }
}

Node BinaryTree::initTree(const BoolVector& dyckWord)
{
  if (!dyckWord.empty())
  {
    assert(isWellFormed(dyckWord));
    assert(dyckWord.front());
    
    Node u = _tree.addNode();
    
    int idx = match(dyckWord, 0);
    BoolVector leftDyckWord(dyckWord.begin() + 1, dyckWord.begin() + idx);
    BoolVector rightDyckWork(dyckWord.begin() + idx + 1, dyckWord.end());
    
    Node v = initTree(leftDyckWord);
    Node w = initTree(rightDyckWork);

    _tree.addArc(u, w);
    _tree.addArc(u, v);
    
    return u;
  }
  else
  {
    Node u = _tree.addNode();
    return u;
  }
}

BoolVector BinaryTree::randomDyckWord(int n)
{
  // See "Generating binary trees at random" by M.D. Atkinson and J.-R. Sack
  // Information Processing Letters 41 (1992) 21--23
  IntVector TwoL(2*n, 0);
  for (int i = 1; i < 2*n; ++i)
  {
    TwoL[i] = TwoL[i-1] + 1;
  }
  
  std::shuffle(TwoL.begin(), TwoL.end(), g_rng);
  
  IntSet L(TwoL.begin(), TwoL.begin() + n);
  
  // X[i] is true iff open paren
  BoolVector X(2*n, false);
  for (int i = 0; i < 2*n; ++i)
  {
    X[i] = L.count(i) == 1;
  }
  
  assert(isBalanced(X));
  BoolVector mapped_X = map(X);
  
//  std::cerr << bracketNotation(X) << " => " << bracketNotation(mapped_X) << std::endl;
  
  return mapped_X;
}

BoolVector BinaryTree::map(const BoolVector& w)
{
  if (w.empty())
  {
    return BoolVector();
  }
  
  // find irreducible u (prefix of w)
  int sum = 0;
  int i = 0;
  for (; i < w.size(); ++i)
  {
    if (w[i])
      ++sum;
    else
      --sum;
    
    if (sum == 0)
    {
      ++i;
      break;
    }
  }
  
  BoolVector u(w.begin(), w.begin() + i);
  BoolVector v(w.begin() + i, w.end());
  
  assert(!u.empty());
  if (isWellFormed(u))
  {
    BoolVector phi_w = u;
    BoolVector phi_v = map(v);
    phi_w.insert(phi_w.end(), phi_v.begin(), phi_v.end());
    return phi_w;
  }
  else
  {
    assert(!u.front() && u.back());
    assert(u.size() >= 2);
    BoolVector t(u.begin() + 1, u.end() - 1);
    
    BoolVector t_star;
    std::reverse_copy(t.begin(), t.end(), std::inserter(t_star, t_star.begin()));
    
    BoolVector phi_w;
    phi_w.push_back(true);
    BoolVector phi_v = map(v);
    phi_w.insert(phi_w.end(), phi_v.begin(), phi_v.end());
    phi_w.push_back(false);
    phi_w.insert(phi_w.end(), t_star.begin(), t_star.end());
    
    return phi_w;
  }
  
  return BoolVector();
}

bool BinaryTree::isWellFormed(const BoolVector& w)
{
  int defect = 0;
  
  int sum = 0;
  for (int i = 0; i < w.size(); ++i)
  {
    if (w[i])
      ++sum;
    else
      --sum;
    
    if (i % 2 == 0 && sum < 0)
      ++defect;
  }
  
  return defect == 0;
}

bool BinaryTree::isBalanced(const BoolVector& w)
{
  int sum = 0;
  
  for (int i = 0; i < w.size(); ++i)
  {
    if (w[i])
      ++sum;
    else
      --sum;
  }
  
  return sum == 0;
}

std::string BinaryTree::bracketNotation(const BoolVector& w)
{
  std::string str;
  for (int i = 0; i < w.size(); ++i)
  {
    if (w[i])
      str += "(";
    else
      str += ")";
  }
  return str;
}

BoolVector BinaryTree::dyckWord(Node v) const
{
  if (isLeaf(v))
  {
    return BoolVector();
  }
  else
  {
    BoolVector res;
    
    BoolVector res_left = dyckWord(left(v));
    BoolVector res_right = dyckWord(right(v));
    
    res.push_back(true);
    res.insert(res.end(), res_left.begin(), res_left.end());
    res.push_back(false);
    res.insert(res.end(), res_right.begin(), res_right.end());
    
    return res;
  }
}

int BinaryTree::match(const BoolVector& w, int i)
{
  assert(0 <= i && i < w.size());
  assert(w[i]);
  
  int sum = 0;
  
  for (; i < w.size(); ++i)
  {
    if (w[i])
      ++sum;
    else
      --sum;
    
    if (sum == 0)
      return i;
  }
  
  assert(false);
  return -1;
}
