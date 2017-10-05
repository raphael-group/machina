/*
 *  perfectphylomatrix.cpp
 *
 *   Created on: 25-sep-2015
 *       Author: M. El-Kebir
 */

#include "perfectphylomatrix.h"

namespace gm {
  
PerfectPhyloMatrix::PerfectPhyloMatrix()
  : _n(0)
  , _k(0)
  , _A()
{
}
  
PerfectPhyloMatrix::PerfectPhyloMatrix(int n, int k)
  : _n(n)
  , _k(k)
  , _A(n, IntMatrix(k, IntVector()))
{
}
  
std::ostream& operator<<(std::ostream& out, const PerfectPhyloMatrix& A)
{
  const int n = A.n();
  const int k = A.k();
  
  out << A.n() << " #n" << std::endl;
  out << A.k() << " #k" << std::endl;
  
  for (int i = 0; i < k; ++i)
  {
    for (int c = 0; c < n; ++c)
    {
      if (c > 0 && i == 0)
      {
        // print root row (*,0) only once
        continue;
      }
      
      if (A.defined(c, i))
      {
        bool first = true;
        for (int d = 0; d < A.n(); ++d)
        {
          if (!first)
            out << " ";
          else
            first = false;
          
          out << A(c, i, d);
        }
        out << " ";
      }
      out << "#(" << c << "," << i << ")" << std::endl;
    }
  }
  
  return out;
}
  
std::istream& operator>>(std::istream& in, PerfectPhyloMatrix& A)
{
  A._n = A._k = -1;
  
  std::string line;
  getline(in, line);
  std::stringstream ss(line);
  ss >> A._n;
  
  if (A._n <= 0)
  {
    throw std::runtime_error("Error: n should be nonnegative");
  }
  
  getline(in, line);
  ss.clear();
  ss.str(line);
  ss >> A._k;
  
  if (A._k <= 0)
  {
    throw std::runtime_error("Error: k should be nonnegative");
  }
  
  A._A = PerfectPhyloMatrix::StlInt3Matrix(A._n,
                                           IntMatrix(A._k,
                                                     IntVector(A._n, 0)));
  
  for (int i = 0; i < A._k; ++i)
  {
    for (int c = 0; c < A._n; ++c)
    {
      if (c > 0 && i == 0)
      {
        // read in root row (*,0) only once
        continue;
      }
      
      getline(in, line);
      if (line[0] != '#')
      {
        ss.clear();
        ss.str(line);
        
        for (int d = 0; d < A._n; ++d)
        {
          int j = -1;
          ss >> j;
          
          if (!(0 <= j && j < A._k))
          {
            throw std::runtime_error("Invalid entry");
          }
          
          A.set(c, i, d, j);
        }
      }
      else
      {
        A._A[c][i].clear();
      }
    }
  }
  
  return in;
}
  
} // namespace gm
