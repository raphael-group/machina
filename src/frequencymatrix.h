/*
 * frequencymatrix.h
 *
 *  Created on: 13-apr-2017
 *      Author: M. El-Kebir
 */

#ifndef FREQUENCYMATRIX_H
#define FREQUENCYMATRIX_H

#include "utils.h"

/// This class models a frequency confidence interval matrix
class FrequencyMatrix
{
public:
  /// Construct
  FrequencyMatrix();
  
  /// Return number of samples
  int getNrSamples() const
  {
    return _m;
  }
  
  /// Return number of characters
  int getNrCharacters() const
  {
    return _n;
  }
  
  /// Return whether given label corresponds to a sample
  ///
  /// @param sStr Sample label
  bool isSample(const std::string& sStr) const
  {
    return _sampleToIndex.count(sStr) == 1;
  }
  
  /// Return whether given label corresponds to a character
  ///
  /// @param cStr Character label
  bool isCharacter(const std::string& cStr) const
  {
    return _characterToIndex.count(cStr) == 1;
  }
  
  /// Return label of given sample index
  ///
  /// @param s Sample index
  const std::string& indexToSample(int s) const
  {
    assert(0 <= s && s < _m);

    return _indexToSample[s];
  }
  
  /// Return index of given sample label
  ///
  /// @param sStr Sample label
  int sampleToIndex(const std::string& sStr) const
  {
    assert(_sampleToIndex.count(sStr) == 1);
    
    return _sampleToIndex.find(sStr)->second;
  }
  
  /// Return label of given character index
  ///
  /// @param c Character index
  const std::string& indexToCharacter(int c) const
  {
    assert(0 <= c && c < _n);
    
    return _indexToCharacter[c];
  }
  
  /// Return index of given character
  ///
  /// @param cStr Character label
  int characterToIndex(const std::string& cStr) const
  {
    assert(_characterToIndex.count(cStr) == 1);
    
    return _characterToIndex.find(cStr)->second;
  }
  
  /// Return frequency lower bound
  ///
  /// @param s Sample index
  /// @param i Character index
  double min(int s, int i) const
  {
    assert(0 <= s && s < _m);
    assert(0 <= i && i < _n);
    
    return _f[s][i].first;
  }

  /// Return frequency upper bound
  ///
  /// @param s Sample index
  /// @param i Character index
  double max(int s, int i) const
  {
    assert(0 <= s && s < _m);
    assert(0 <= i && i < _n);
    
    return _f[s][i].second;
  }

  /// Return sample labels
  StringSet getSamples() const
  {
    StringSet Sigma;
    for (int s = 0; s < _m; ++s)
    {
      Sigma.insert(_indexToSample[s]);
    }
    return Sigma;
  }
  
private:
  typedef std::pair<double, double> DoublePair;
  typedef std::vector<DoublePair> DoublePairVector;
  typedef std::vector<DoublePairVector> DoublePairMatrix;
  
private:
  /// Number of samples
  int _m;
  /// Number of characters
  int _n;
  /// Index to sample label
  StringVector _indexToSample;
  /// Sample label to index
  StringToIntMap _sampleToIndex;
  /// Index to character label
  StringVector _indexToCharacter;
  /// Character label to index
  StringToIntMap _characterToIndex;
  /// Frequencies
  DoublePairMatrix _f;
  
  friend std::ostream& operator<<(std::ostream& out, const FrequencyMatrix& F);
  friend std::istream& operator>>(std::istream& in, FrequencyMatrix& F);
};

std::ostream& operator<<(std::ostream& out, const FrequencyMatrix& F);
std::istream& operator>>(std::istream& in, FrequencyMatrix& F);

#endif // FREQUENCYMATRIX_H
