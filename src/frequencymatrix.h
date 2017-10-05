/*
 * frequencymatrix.h
 *
 *  Created on: 13-apr-2017
 *      Author: M. El-Kebir
 */

#ifndef FREQUENCYMATRIX_H
#define FREQUENCYMATRIX_H

#include "utils.h"
#include "basematrix.h"

/// This class models a frequency confidence interval matrix
class FrequencyMatrix : public BaseMatrix
{
public:
  /// Constructor
  FrequencyMatrix();
  
  /// Constructor
  FrequencyMatrix(const BaseMatrix& other);
  
  /// Return frequency lower bound
  ///
  /// @param p Sample index
  /// @param i Character index
  double min(int p, int i) const
  {
    assert(0 <= p && p < _k);
    assert(0 <= i && i < _n);
    
    return _f[p][i].first;
  }

  /// Return frequency upper bound
  ///
  /// @param p Sample index
  /// @param i Character index
  double max(int p, int i) const
  {
    assert(0 <= p && p < _k);
    assert(0 <= i && i < _n);
    
    return _f[p][i].second;
  }
  
  /// Set frequency
  ///
  /// @param p Sample index
  /// @param i Character index
  /// @param f_lb Frequency lower bound
  /// @param f_ub Frequency upper bound
  void set(int p, int i, double f_lb, double f_ub)
  {
    assert(0 <= p && p < _k);
    assert(0 <= i && i < _n);
    assert(0 <= f_lb && f_lb <= f_ub && f_ub <= 1);
    _f[p][i].first = f_lb;
    _f[p][i].second = f_ub;
  }
  
  /// Return whether a mutation is surely subclonal
  ///
  /// @param s Anatomical site index
  /// @param i Character index
  bool isSurelySubclonal(int s, int i) const;
  
  /// Return whether a mutation is surely present
  ///
  /// @param s Anatomical site index
  /// @param i Character index
  bool isSurelyPresent(int s, int i) const;
  
  /// Return whether a mutation is surely descendant
  ///
  /// @param s Anatomical site index
  /// @param i Character index
  /// @param P Primary index
  bool isSurelyDescendant(int s, int i, int P) const;
  
  /// Return whether (F^-, F+^) satisfy a necessary condition for mS
  ///
  /// @param P Primary index
  bool mS(int P) const;
  
private:
  typedef std::pair<double, double> DoublePair;
  typedef std::vector<DoublePair> DoublePairVector;
  typedef std::vector<DoublePairVector> DoublePairMatrix;
  
private:
  /// Frequencies
  DoublePairMatrix _f;
  
  friend std::ostream& operator<<(std::ostream& out, const FrequencyMatrix& F);
  friend std::istream& operator>>(std::istream& in, FrequencyMatrix& F);
};

/// Writes a frequency matrix from the specified output stream
///
/// @param out Output stream
/// @param F Frequency matrix
std::ostream& operator<<(std::ostream& out, const FrequencyMatrix& F);

/// Reads a frequency matrix from the specified input stream
///
/// @param in Input stream
/// @param F Frequency matrix
std::istream& operator>>(std::istream& in, FrequencyMatrix& F);

#endif // FREQUENCYMATRIX_H
