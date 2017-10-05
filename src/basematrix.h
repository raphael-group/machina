/*
 * basematrix.h
 *
 *  Created on: 5-sep-2017
 *      Author: M. El-Kebir
 */

#ifndef BASEMATRIX_H
#define BASEMATRIX_H

#include "utils.h"

/// This class models an k * m * n matrix
class BaseMatrix
{
public:
  /// Constructor
  BaseMatrix();
  
  /// Return number of anatomical sites
  int getNrAnatomicalSites() const
  {
    return _m;
  }
  
  /// Return number of samples
  int getNrSamples() const
  {
    return _k;
  }
  
  /// Return number of characters
  int getNrCharacters() const
  {
    return _n;
  }
  
  /// Return whether given label corresponds to a sample
  ///
  /// @param pStr Sample label
  bool isSample(const std::string& pStr) const
  {
    return _sampleToIndex.count(pStr) == 1;
  }
  
  /// Return whether given label corresponds to an anatomical site
  ///
  /// @param sStr Anatomical site label
  bool isAnatomicalSite(const std::string& sStr) const
  {
    return _anatomicalSiteToIndex.count(sStr) == 1;
  }
  
  /// Return whether given label corresponds to a character
  ///
  /// @param cStr Character label
  bool isCharacter(const std::string& cStr) const
  {
    return _characterToIndex.count(cStr) == 1;
  }
  
  /// Return label of given anatomical site index
  ///
  /// @param s Anatomical site index
  const std::string& indexToAnatomicalSite(int s) const
  {
    assert(0 <= s && s < _m);
    
    return _indexToAnatomicalSite[s];
  }
  
  /// Return index of given anatomical site label, -1 is returned if provided
  /// anatomical site is undefined
  ///
  /// @param sStr Anatomical site label
  int anatomicalSiteToIndex(const std::string& sStr) const
  {
    if (_anatomicalSiteToIndex.count(sStr) == 1)
      return _anatomicalSiteToIndex.find(sStr)->second;
    else
      return -1;
  }
  
  /// Return label of given sample index
  ///
  /// @param p Sample index
  const std::string& indexToSample(int p) const
  {
    assert(0 <= p && p < _k);
    
    return _indexToSample[p];
  }
  
  /// Return index of given sample label, -1 is returned if provided
  /// sample is undefined
  ///
  /// @param pStr Sample label
  int sampleToIndex(const std::string& pStr) const
  {
    if (_sampleToIndex.count(pStr) == 1)
      return _sampleToIndex.find(pStr)->second;
    else
      return -1;
  }
  
  /// Return label of given character index
  ///
  /// @param c Character index
  const std::string& indexToCharacter(int c) const
  {
    assert(0 <= c && c < _n);
    
    return _indexToCharacter[c];
  }
  
  /// Return index of given character.
  /// If the provided does not correspond to a character, -1 is returned.
  ///
  /// @param cStr Character label
  int characterToIndex(const std::string& cStr) const
  {
    if (_characterToIndex.count(cStr) == 1)
      return _characterToIndex.find(cStr)->second;
    else
      return -1;
  }
  
  /// Return anatomical site labels
  StringSet getAnatomicalSites() const
  {
    StringSet Sigma;
    for (int s = 0; s < _m; ++s)
    {
      Sigma.insert(_indexToAnatomicalSite[s]);
    }
    return Sigma;
  }
  
  /// Return anatomical site index associated with provided sample index
  ///
  /// @param p Sample index
  int sampleIndexToAnatomicalSiteIndex(int p) const
  {
    assert(0 <= p && p < _k);
    
    return _sampleIndexToAnatomicalSiteIndex[p];
  }
  
  /// Return set of sample indices associated with provided anatomical site index
  ///
  /// @param s Anatomical site index
  const IntSet& anatomicalSiteIndexToSampleIndices(int s) const
  {
    assert(0 <= s && s < _m);
    
    return _anatomicalSiteIndexToSampleIndices[s];
  }
  
  /// Return index to anatomical sites vector
  const StringVector& getIndexToAnatomicalSites() const
  {
    return _indexToAnatomicalSite;
  }
  
  /// Return map of sample indices to anatomical site indices
  const IntVector& getSampleIndexToAnatomicalSiteIndex() const
  {
    return _sampleIndexToAnatomicalSiteIndex;
  }
  
  /// Return map of anatomical site indices to sample indices
  const IntSetVector& getAnatomicalSiteIndexToSampleIndices() const
  {
    return _anatomicalSiteIndexToSampleIndices;
  }
  
protected:
  /// Number of anatomical sites
  int _m;
  /// Number of samples
  int _k;
  /// Number of two-state characters
  int _n;
  /// Index to sample label
  StringVector _indexToAnatomicalSite;
  /// Sample label to index
  StringToIntMap _anatomicalSiteToIndex;
  /// Index to sample label
  StringVector _indexToSample;
  /// Sample label to index
  StringToIntMap _sampleToIndex;
  /// Index to character label
  StringVector _indexToCharacter;
  /// Character label to index
  StringToIntMap _characterToIndex;
  /// Sample index to anatomical site index
  IntVector _sampleIndexToAnatomicalSiteIndex;
  /// Anatomical site index to sample indices
  IntSetVector _anatomicalSiteIndexToSampleIndices;
};

#endif // BASEMATRIX_H
