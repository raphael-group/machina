/*
 * frequencymatrix.cpp
 *
 *  Created on: 13-apr-2017
 *      Author: M. El-Kebir
 */

#include "frequencymatrix.h"

FrequencyMatrix::FrequencyMatrix()
  : _m(0)
  , _n(0)
  , _indexToSample()
  , _sampleToIndex()
  , _indexToCharacter()
  , _characterToIndex()
  , _f()
{
}

std::ostream& operator<<(std::ostream& out, const FrequencyMatrix& F)
{
  out << F._m << " #samples" << std::endl;
  out << F._n << " #characters" << std::endl;
  out << "#sample_index\tsample_label\tcharacter_index\tcharacter_label\tf-\tf+" << std::endl;
  for (int s = 0; s < F._m; ++s)
  {
    const std::string& sStr = F.indexToSample(s);
    for (int c = 0; c < F._n; ++c)
    {
      const std::string& cStr = F.indexToCharacter(c);
      out << s << "\t" << sStr << "\t" << c << "\t"
          << cStr << "\t" << F._f[s][c].first << "\t"
          << F._f[s][c].second << std::endl;
    }
  }
  
  return out;
}

std::istream& operator>>(std::istream& in, FrequencyMatrix& F)
{
  std::string line;
  getline(in, line);
  
  int m = -1;
  int n = -1;
  
  std::stringstream ss(line);
  ss >> m;
  
  if (m <= 0)
  {
    throw std::runtime_error("Error: m should be nonnegative");
  }
  
  getline(in, line);
  ss.clear();
  ss.str(line);
  ss >> n;
  
  if (n <= 0)
  {
    throw std::runtime_error("Error: n should be nonnegative");
  }
  
  F._m = m;
  F._n = n;
  
  F._indexToSample = StringVector(m);
  F._indexToCharacter = StringVector(n);
  F._f = FrequencyMatrix::DoublePairMatrix(m, FrequencyMatrix::DoublePairVector(n, std::make_pair(0., 0.)));
  
  std::vector<std::vector<bool> > present(m, std::vector<bool>(n, false));
  while (in.good())
  {
    getline(in, line);
    if (line == "" || line[0] == '#')
      continue;
    
    ss.clear();
    ss.str(line);
    
    int s = -1;
    std::string sStr;
    int c = -1;
    std::string cStr;
    double f_min = -1;
    double f_max = -1;
    
    ss >> s >> sStr >> c >> cStr >> f_min >> f_max;
    
    if (!(0 <= s && s < m) || !(0 <= c && c < n)
        || !(0 <= f_min && f_min <= 1) || !(0 <= f_max && f_max <= 1))
    {
      throw std::runtime_error("Invalid character ("
                               + boost::lexical_cast<std::string>(s)
                               + ","
                               + boost::lexical_cast<std::string>(c)
                               + ")");
    }
    
    if (present[s][c])
    {
      throw std::runtime_error("Duplicate character ("
                               + boost::lexical_cast<std::string>(s)
                               + ","
                               + boost::lexical_cast<std::string>(c)
                               + ")");
    }
    
    F._indexToSample[s] = sStr;
    F._sampleToIndex[sStr] = s;
    F._indexToCharacter[c] = cStr;
    F._characterToIndex[cStr] = c;
    F._f[s][c] = std::make_pair(f_min, f_max);
    
    present[s][c] = true;
  }
  
  return in;
}
