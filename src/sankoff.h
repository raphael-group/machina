/*
 * sankoff.h
 *
 *  Created on: 10-jan-2017
 *      Author: M. El-Kebir
 */

#ifndef SANKOFF_H
#define SANKOFF_H

#include "utils.h"
#include "charactertree.h"
#include "mpbase.h"
#include <map>

/// This class implements the Sankoff algorithm with an extension to
/// enumerate all maximum parsimony solutions
/// @brief This class implements the Sankoff algorithm
class Sankoff : public MPBase
{
public:
  /// Constructor
  ///
  /// @param T Character-based tree T
  Sankoff(const CharacterTree& T);
  
  /// Destructor
  virtual ~Sankoff()
  {
    int nrCharacters = _T.getNrCharacters();
    for (int c = 0; c < nrCharacters; ++c)
    {
      delete _M[c];
    }
  }
  
  /// Run the enumeration
  ///
  /// @param rootState Denotes the root state, when set to -1 all
  /// possible states are considered.
  void run(int rootState = -1);
  
  /// Write the solution for the given character
  ///
  /// @param c Character
  /// @param out Output stream
  void write(int c, std::ostream& out) const;
  
protected:
  typedef std::pair<Digraph::Node, int> NodeStatePair;
  typedef std::vector<NodeStatePair> NodeStatePairVector;
  typedef std::map<Node, IntSet> NodeIntSetMap;

  /// Dynamic programming entry
  struct Entry {
    /// Default constructor
    Entry()
      : _cost(std::numeric_limits<int>::max())
      , _previousEntries()
    {
    }
    
    /// Number of state changes
    int _cost;
    /// Previous entries
    NodeIntSetMap _previousEntries;
  };
  
  typedef std::vector<Entry> EntryVector;
  typedef Digraph::NodeMap<EntryVector> EntryMatrix;
  typedef std::vector<EntryMatrix*> CharacterTable;
  
  /// Dynamic programming table
  CharacterTable _M;
  
  /// Solve for the given character and node
  ///
  /// @param characterIndex Character
  /// @param u Node
  void run(int characterIndex, Node u);
  
  /// Return the dynamic programming entry corresponding to (c, u, s)
  ///
  /// @param c Character
  /// @param u Node
  /// @param s State
  Entry& entry(int c, Node u, int s)
  {
    return (*_M[c])[u][s];
  }

  /// Return the dynamic programming entry corresponding to (c, u, s)
  ///
  /// @param c Character
  /// @param u Node
  /// @param s State
  const Entry& entry(int c, Node u, int s) const
  {
    return (*_M[c])[u][s];
  }
  
  /// Compute the cost of reaching (c, v, t) from (c, u, s)
  ///
  /// @param c Character
  /// @param u Node
  /// @param s State
  /// @param v Node
  /// @param t State
  int computeCost(int c, Node u, int s, Node v, int t) const
  {
    int cost_vt = entry(c, v, t)._cost;
    if (s != t)
    {
      if (cost_vt == std::numeric_limits<int>::max())
      {
        return cost_vt;
      }
      else
      {
        return cost_vt + 1;
      }
    }
    else
    {
      return cost_vt;
    }
  }
  
  /// Construct a back trace from (c, u, s)
  ///
  /// @param c Character
  /// @param u Node
  /// @param s State
  void constructBackTrace(int c, Node u, int s);
  
  /// Construct all back traces
  ///
  /// @param c Character
  /// @param stateVector Set of states
  /// @param frontier Frontier of node state pairs leading to optimal solutions
  void constructBackTrace(int c,
                          IntNodeMap& stateVector,
                          NodeStatePairList& frontier);
};

#endif // SANKOFF_H
