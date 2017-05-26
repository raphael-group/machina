/*
 * mpbase.h
 *
 *  Created on: 10-jan-2017
 *      Author: M. El-Kebir
 */

#ifndef MPBASE_H
#define MPBASE_H

#include "utils.h"
#include "charactertree.h"

/// This class models maximum parsimony small phylogeny enumeration algorithms
class MPBase
{
public:
  /// Constructor
  ///
  /// @param T Character-based tree T
  MPBase(const CharacterTree& T);
  
  /// Destructor
  virtual ~MPBase()
  {
    const int nrCharacters = _T.getNrCharacters();
    
    for (int c = 0; c < nrCharacters; ++c)
    {
      int nrSolutions = _stateVector[c].size();
      for (int solIdx = 0; solIdx < nrSolutions; ++solIdx)
      {
        delete _stateVector[c][solIdx];
      }
    }
  };
  
  /// Print search tree in DOT format
  void writeDOT(std::ostream& out) const;
  
  /// Return the state of the given node and character
  int state(Node u, int c) const
  {
    assert(0 <= c && c < _T.getNrCharacters());
    
    int solIdx = _solutionIndex[c];
    return (*_stateVector[c][solIdx])[u];
  }
  
  /// Return whether the given character underwent homoplasy
  bool homoplasy(int c) const
  {
    assert(0 <= c && c < _T.getNrCharacters());

    int solIdx = _solutionIndex[c];
    return _homoplasyVector[c][solIdx];
  }
  
  /// Return the set of characters that underwent homoplasy
  IntSet getHomoplasyCharacters() const
  {
    const int nrCharacters = _T.getNrCharacters();

    IntSet result;
    for (int c = 0; c < nrCharacters; ++c)
    {
      if (homoplasy(c))
      {
        result.insert(c);
      }
    }
    return result;
  }
  
  /// Return the number of solutions
  unsigned long long getNrSolutions() const
  {
    const int nrCharacters = _T.getNrCharacters();
    
    int res = 1;
    for (int c = 0; c < nrCharacters; ++c)
    {
      assert(_stateVector[c].size() > 0);
      res = res * _stateVector[c].size();
    }
    return res;
  }
  
  /// Switch to the next solution
  bool nextSolution();
  
  /// Switch to the first solution
  void resetSolutions()
  {
    const int nrCharacters = _T.getNrCharacters();
    
    for (int c = 0; c < nrCharacters; ++c)
    {
      _solutionIndex[c] = 0;
    }
  }
  
  /// Run the enumeration
  ///
  /// @param rootState Denotes the root state, when set to -1 all possible states are considered.
  virtual void run(int rootState = -1) = 0;
  
protected:
  typedef std::pair<Node, int> NodeStatePair;
  typedef std::list<NodeStatePair> NodeStatePairList;
  
  typedef std::vector<IntNodeMap*> IntNodeMapVector;
  typedef std::vector<IntNodeMapVector> IntNodeMapVectorVector;
  
  typedef std::vector<IntVectorNodeMap*> IntVectorNodeMapVector;
  typedef std::vector<BoolVector> BoolVectorVector;
  
  /// Underlying character-based tree whose leaves are labeled
  const CharacterTree& _T;
  
  /// Character-specific solution index
  IntVector _solutionIndex;
  
  /// Character-specific state vector
  // (*_stateVector[c][solIdx])[u]
  IntNodeMapVectorVector _stateVector;
  
  /// Character-specific homoplasy vector
  // _homoplasyVector[c][solIdx]
  BoolVectorVector _homoplasyVector;
  
  /// Update homoplasy vector
  void updateHomoplasy();
};

#endif // MPBASE_H
