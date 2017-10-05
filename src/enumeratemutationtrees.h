/*
 * enumeratemutationtrees.h
 *
 *  Created on: 07-sep-2017
 *      Author: M. El-Kebir
 */

#ifndef ENUMERATEMUTATIONTREES_H
#define ENUMERATEMUTATIONTREES_H

#include "frequencymatrix.h"
#include "clonetree.h"
#include "spruce/realtensor.h"

/// This class enumerates all mutation trees given a frequency matrix
class EnumerateMutationTrees
{
public:
  typedef std::vector<CloneTree> TreeVector;
  
  /// Constructor
  ///
  /// @param F Frequency matrix
  EnumerateMutationTrees(const FrequencyMatrix& F);
  
  /// Enumerate mutation trees
  ///
  /// @param outputDirectory Output directory, may be empty
  /// @param nrThreads Number of threads to use in the enumeration
  /// @param limit Maximum number of trees to enumerate
  /// @param timeLimit Time limit for the enumeration
  /// @param mutationTrees Vector to store enumerated mutation trees
  void enumerate(const std::string& outputDirectory,
                 const int nrThreads,
                 const int limit,
                 const int timeLimit,
                 TreeVector& mutationTrees);
  
  /// Return number of mutations in the enumerated trees
  int getNrCharactersInTrees() const
  {
    return _nrCharactersInTrees;
  }
  
  /// Pick trees uniformly at random
  ///
  /// @param trees Mutation trees
  /// @param nrTrees to pick
  static void pick(TreeVector& trees,
                   int nrTrees)
  {
    IntVector indices(trees.size(), 0);
    for (int i = 1; i < indices.size(); ++i)
    {
      indices[i] = indices[i-1] + 1;
    }
    std::shuffle(indices.begin(), indices.end(), g_rng);
    
    if (nrTrees > trees.size())
    {
      nrTrees = trees.size();
    }
    
    TreeVector pickedTrees;
    pickedTrees.reserve(nrTrees);
    for (int i = 0; i < nrTrees; ++i)
    {
      pickedTrees.push_back(trees[indices[i]]);
    }
    
    trees = pickedTrees;
  }
  
private:
  /// Construct frequency tensor (dimension: 2 * k * n)
  ///
  /// @param F_lb Frequency lower bounds
  /// @param F_ub Frequency upper bounds
  void getFrequencyTensor(RealTensor& F_lb,
                          RealTensor& F_ub) const;
  
private:
  /// Frequency matrix
  const FrequencyMatrix& _F;
  /// Number of mutations in the enumerated trees
  int _nrCharactersInTrees;
};

/// Write mutation trees to the specified output stream
///
/// @param out Output stream
/// @param trees Mutation trees
std::ostream& operator<<(std::ostream& out,
                         const EnumerateMutationTrees::TreeVector& trees);

/// Read mutation trees from the specified input stream
///
/// @param in Input stream
/// @param trees Mutation trees
std::istream& operator>>(std::istream& in,
                         EnumerateMutationTrees::TreeVector& trees);

#endif // ENUMERATEMUTATIONTREES_H
