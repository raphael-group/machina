/*
 * sankofflabeling.h
 *
 *  Created on: 12-jan-2017
 *      Author: M. El-Kebir
 */

#ifndef SANKOFFLABELING_H
#define SANKOFFLABELING_H

#include "sankoff.h"
#include "clonetree.h"
#include "migrationgraph.h"

/// This class is a wrapper for solving the Parsimonious Migration Problem
/// with a minimum migration objective.
/// @brief This class solves the PMH problem with minimum number of migrations
class SankoffLabeling
{
public:
  typedef Digraph::NodeMap<std::string> StringNodeMap;
  typedef std::pair<int, std::pair<int, std::pair<bool, int> > > IntTriple;
  typedef std::map<IntTriple, int> IntTripleToIntMap;
  
public:
  /// Constructor
  ///
  /// @param T Clone tree
  /// @param primary Primary tumor label
  SankoffLabeling(const CloneTree& T,
                  const std::string& primary)
    : _T(T)
    , _primary(primary)
    , _stateToSample()
    , _sampleToState()
    , _charT(_T, primary, _stateToSample, _sampleToState)
  {
  }
  
  /// Destructor
  ~SankoffLabeling()
  {
    for (StringNodeMap* pLabeling : _mpLabelings)
    {
      delete pLabeling;
    }
  }
  
  /// Run the enumeration
  void run();
  
  /// Return the number of labelings
  int getNrLabelings() const
  {
    return _mpLabelings.size();
  }
  
  /// Classify the results into (gamma, sigma, pattern) triples
  IntTripleToIntMap classify() const
  {
    IntTripleToIntMap result;
    
    const int nrSolutions = getNrLabelings();
    for (int solIdx = 0; solIdx < nrSolutions; ++solIdx)
    {
      MigrationGraph G(_T, getLabeling(solIdx));
      int nrComigrations = G.getNrComigrations(_T, getLabeling(solIdx));
      int nrSeedingSamples = G.getNrSeedingSites();
      int migrationPattern = static_cast<int>(G.getPattern());
      IntTriple triple = std::make_pair(nrComigrations,
                                        std::make_pair(nrSeedingSamples,
                                                       std::make_pair(G.isMonoclonal(),
                                                                      migrationPattern)));
      
      if (result.count(triple) == 0)
      {
        result[triple] = 1;
      }
      else
      {
        ++result[triple];
      }
    }
    
    return result;
  }
  
  /// Return labeling corresponding to the given solution idnex
  ///
  /// @param solIdx Solution index
  const StringNodeMap& getLabeling(int solIdx) const
  {
    assert(0 <= solIdx && solIdx < _mpLabelings.size());
    return *(_mpLabelings[solIdx]);
  }
  
  /// Return migration graph corresponding to the given solution index
  ///
  /// @param solIdx Solution index
  MigrationGraph getMigrationGraph(int solIdx) const
  {
    return MigrationGraph(_T, getLabeling(solIdx));
  }
  
  /// Run Sankoff enumeration algorithm to infer all vertex labelings
  /// with minimum number of migrations
  ///
  /// @param T Clone tree
  /// @param primary Primary tumor label
  /// @param outputDirectory Output directory
  /// @param colorMap Color map
  static void run(const CloneTree& T,
                  const std::string& primary,
                  const std::string& outputDirectory,
                  const StringToIntMap& colorMap);
  
private:
  typedef CharacterTree::StringToIntMap StringToIntMap;
  typedef std::vector<StringNodeMap*> StringNodeMapVector;
  
private:
  /// Clone tree
  const CloneTree& _T;
  /// Primary tumor label
  const std::string& _primary;
  /// Mapping from state to anatomical state label
  StringVector _stateToSample;
  /// Mapping from anatomical state label to state
  StringToIntMap _sampleToState;
  /// Character tree corresponding to _T
  const CharacterTree _charT;
  /// Maximum parsimony (minimum migration) vertex labelings
  StringNodeMapVector _mpLabelings;
};

#endif // SANKOFFLABELING_H
