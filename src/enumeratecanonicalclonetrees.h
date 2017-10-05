/*
 * enumeratecanonicalclonetrees.h
 *
 *  Created on: 07-sep-2017
 *      Author: M. El-Kebir
 */

#ifndef ENUMERATECANONICALCLONETREES_H
#define ENUMERATECANONICALCLONETREES_H

#include "frequencymatrix.h"
#include "clonetree.h"
#include "spruce/solution.h"
#include "spruce/realtensor.h"

/// This class enumerates all canonical clone trees given a frequency matrix
class EnumerateCanonicalCloneTrees
{
public:
  typedef std::vector<CloneTree> TreeVector;
  
  /// Constructor
  ///
  /// @param F Frequency matrix
  EnumerateCanonicalCloneTrees(const FrequencyMatrix& F);
  
  /// Enumerate canonical clone trees
  ///
  /// @param outputDirectory Output directory, may be empty
  /// @param nrThreads Number of threads to use in the enumeration
  /// @param limit Maximum number of trees to enumerate
  /// @param canonicalCloneTrees Vector to store enumerated canonical clone trees
  void enumerate(const std::string& outputDirectory,
                 const int nrThreads,
                 const int limit,
                 TreeVector& canonicalCloneTrees) const;
  
private:
  struct CanonicalCloneTree
  {
  public:
    CanonicalCloneTree(const Digraph& tree,
                       Node root,
                       const StringNodeMap& idMap,
                       const StringNodeMap& leafLabelMap)
    : _T(tree, root, idMap, leafLabelMap)
    , _usage(_T.tree())
    , _frequency(_T.tree())
    , _characterNodeMap(_T.tree())
    {
    }
    
    CloneTree _T;
    DoubleVectorNodeMap _usage;
    DoubleVectorNodeMap _frequency;
    StringNodeMap _characterNodeMap;
  };
  
  CanonicalCloneTree* constructCanonicalCloneTree(const gm::Solution& solution) const;
  
  /// Construct frequency tensor (dimension: 2 * k * n)
  ///
  /// @param F_lb Frequency lower bounds
  /// @param F_ub Frequency upper bounds
  void getFrequencyTensor(RealTensor& F_lb,
                          RealTensor& F_ub) const;
  
private:
  /// Frequency matrix
  const FrequencyMatrix& _F;
};

#endif // ENUMERATECANONICALCLONETREES_H
