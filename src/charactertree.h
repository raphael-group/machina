/*
 * charactertree.h
 *
 *  Created on: 13-oct-2016
 *      Author: M. El-Kebir
 */

#ifndef CHARACTERTREE_H
#define CHARACTERTREE_H

#include "utils.h"
#include "basetree.h"
#include "nonbinaryclonetree.h"

/// This class models a character-based (phylogenetic) tree
class CharacterTree : public BaseTree
{
public:
  typedef std::map<std::string, int> StringToIntMap;

public:
  /// Default construct
  CharacterTree();
  
  /// Copy constructor
  CharacterTree(const CharacterTree& other);
  
  /// Constructor
  ///
  /// @param T Non binary clone tree
  /// @param primary Primary tumor anatomical site label
  /// @param stateToSample State index to anatomical site label
  /// @param sampleToState Anatomical site label to state index
  CharacterTree(const NonBinaryCloneTree& T,
                const std::string& primary,
                StringVector& stateToSample,
                StringToIntMap& sampleToState);
  
  /// Read leaf labeling
  ///
  /// @param in Input stream
  bool readLeafLabeling(std::istream& in);
  
  /// Return number of characters
  int getNrCharacters() const
  {
    return _indexToCharacter.size();
  }
  
  /// Return number of states for given character
  ///
  /// @param characterIndex Character index
  int getNrStates(int characterIndex) const
  {
    assert(0 <= characterIndex && characterIndex < getNrCharacters());
    return _nrStates[characterIndex];
  }
  
  /// Return character index corresponding to the given string
  ///
  /// @param c Character label
  int characterIndex(const std::string& c) const
  {
    if (_characterToIndex.count(c) == 0)
    {
      return -1;
    }
    else
    {
      return _characterToIndex.find(c)->second;
    }
  }
  
  /// Return character label given an index
  ///
  /// @param characterIndex Character index
  const std::string& characterLabel(int characterIndex) const
  {
    assert(0 <= characterIndex && characterIndex < getNrCharacters());
    return _indexToCharacter[characterIndex];
  }
  
  /// Return the state of the given character in the given node
  ///
  /// @param u Node
  /// @param characterIndex Character index
  int state(Node u, int characterIndex) const
  {
    return _stateVector[u][characterIndex];
  }
  
  /// Remove given leaf
  ///
  /// @param v Node
  void eraseLeaf(Node v);
  
protected:
  /// State vector of a given node
  IntVectorNodeMap _stateVector;
  /// Number of states of a given character
  IntVector _nrStates;
  /// Index to character label
  StringVector _indexToCharacter;
  /// Character label to index
  StringToIntMap _characterToIndex;
};

#endif // CHARACTERTREE_H
