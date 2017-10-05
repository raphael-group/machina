/*
 * frequencymatrix.cpp
 *
 *  Created on: 13-apr-2017
 *      Author: M. El-Kebir
 */

#include "frequencymatrix.h"

FrequencyMatrix::FrequencyMatrix()
  : BaseMatrix()
  , _f()
{
}

FrequencyMatrix::FrequencyMatrix(const BaseMatrix& other)
  : BaseMatrix(other)
  , _f(DoublePairMatrix(_k,
                        DoublePairVector(_n,
                                         std::make_pair(0., 0.))))
{

}


bool FrequencyMatrix::isSurelyPresent(int s, int i) const
{
  assert(0 <= s && s < _m);
  assert(0 <= i && i < _n);
  
  for (int p : _anatomicalSiteIndexToSampleIndices[s])
  {
    if (_f[p][i].first > 0)
      return true;
  }
  
  return false;
}

bool FrequencyMatrix::isSurelySubclonal(int s, int i) const
{
  assert(0 <= s && s < _m);
  assert(0 <= i && i < _n);
  
  if (!isSurelyPresent(s, i))
    return false;

  for (int q : _anatomicalSiteIndexToSampleIndices[s])
  {
    for (int j = 0; j < _n; ++j)
    {
      if (i == j) continue;
      
      if (_f[q][i].second < _f[q][j].first)
        return true;
    }
  }
  
  return false;
}

bool FrequencyMatrix::isSurelyDescendant(int s, int i, int P) const
{
  if (s == P)
  {
    return isSurelyPresent(s, i);
  }
  else
  {
    return isSurelySubclonal(s, i);
  }
}

bool FrequencyMatrix::mS(int P) const
{
  bool res = true;
  
  IntSetVector surelyDescedantSites(_n);
  for (int i = 0; i < _n; ++i)
  {
    for (int s = 0; s < _m; ++s)
    {
        if (isSurelyDescendant(s, i, P))
        {
          surelyDescedantSites[i].insert(s);
        }
    }
    
    if (surelyDescedantSites[i].size() > 2)
    {
      std::cerr << "Mutation '" << _indexToCharacter[i] << "'"
        << " is surely descendant in anatomical sites";
      for (int s : surelyDescedantSites[i])
      {
        std::cerr << " '" << _indexToAnatomicalSite[s] << "'";
      }
      std::cerr << std::endl;
      res = false;
    }
  }

  return res;
}

std::ostream& operator<<(std::ostream& out, const FrequencyMatrix& F)
{
  out << F._m << " #anatomical sites" << std::endl;
  out << F._k << " #samples" << std::endl;
  out << F._n << " #characters" << std::endl;
  out << "#sample_index\tsample_label\tanatomical_site_index\tanatomical_site_label"\
         "\tcharacter_index\tcharacter_label\tf-\tf+"
      << std::endl;
  for (int p = 0; p < F._k; ++p)
  {
    int s = F.sampleIndexToAnatomicalSiteIndex(p);
    const std::string sStr = F.indexToAnatomicalSite(s);
    
    const std::string& pStr = F.indexToSample(p);
    for (int c = 0; c < F._n; ++c)
    {
      const std::string& cStr = F.indexToCharacter(c);
      out << p << "\t" << pStr << "\t"
          << s << "\t" << sStr << "\t"
          << c << "\t" << cStr << "\t"
          << F._f[p][c].first << "\t" << F._f[p][c].second << std::endl;
    }
  }
  
  return out;
}

std::istream& operator>>(std::istream& in, FrequencyMatrix& F)
{
  std::string line;
  getline(in, line);
  
  int m = -1;
  int k = -1;
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
  ss >> k;

  if (k < m)
  {
    throw std::runtime_error("Error: k should be at least m");
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
  F._k = k;
  F._n = n;

  F._indexToAnatomicalSite = StringVector(m);
  F._indexToSample = StringVector(k);
  F._indexToCharacter = StringVector(n);
  F._f = FrequencyMatrix::DoublePairMatrix(k, FrequencyMatrix::DoublePairVector(n, std::make_pair(0., 0.)));
  F._sampleIndexToAnatomicalSiteIndex = IntVector(k, -1);
  F._anatomicalSiteIndexToSampleIndices = IntSetVector(m);
  
  std::vector<std::vector<bool> > present(k, std::vector<bool>(n, false));
  while (in.good())
  {
    getline(in, line);
    if (line == "" || line[0] == '#')
      continue;
    
    ss.clear();
    ss.str(line);
    
    int p = -1;
    std::string pStr;
    int s = -1;
    std::string sStr;
    int c = -1;
    std::string cStr;
    double f_min = -1;
    double f_max = -1;
    
    ss >> p >> pStr >> s >> sStr >> c >> cStr >> f_min >> f_max;
    
    if (!(0 <= p && p < k) || !(0 <= s && s < m) || !(0 <= c && c < n)
        || !(0 <= f_min) || !(0 <= f_max)
        || !(f_min <= f_max))
    {
      throw std::runtime_error("Invalid character ("
                               + boost::lexical_cast<std::string>(p)
                               + ","
                               + boost::lexical_cast<std::string>(c)
                               + ")");
    }
    
    if (present[p][c])
    {
      throw std::runtime_error("Duplicate character ("
                               + boost::lexical_cast<std::string>(p)
                               + ","
                               + boost::lexical_cast<std::string>(c)
                               + ")");
    }

    F._indexToAnatomicalSite[s] = sStr;
    F._anatomicalSiteToIndex[sStr] = s;
    F._indexToSample[p] = pStr;
    F._sampleToIndex[pStr] = p;
    F._indexToCharacter[c] = cStr;
    F._characterToIndex[cStr] = c;
    F._f[p][c] = std::make_pair(f_min, f_max);
    F._sampleIndexToAnatomicalSiteIndex[p] = s;
    F._anatomicalSiteIndexToSampleIndices[s].insert(p);
    
    present[p][c] = true;
  }
  
  return in;
}
