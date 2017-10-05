/*
 * utils.cpp
 *
 *  Created on: 19-oct-2016
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include <fstream>

lemon::Tolerance<double> g_tol(1e-4);

std::mt19937 g_rng(0);

VerbosityLevel g_verbosity = VERBOSE_ESSENTIAL;

int g_lineNumber = 0;

std::string getLineNumber()
{
  char buf[1024];
  
  snprintf(buf, 1024, "Line: %d. ", g_lineNumber);
  
  return std::string(buf);
}

std::istream& getline(std::istream& is, std::string& t)
{
  ++g_lineNumber;
  
  // source: http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
  t.clear();
  
  // The characters in the stream are read one-by-one using a std::streambuf.
  // That is faster than reading them one-by-one using the std::istream.
  // Code that uses streambuf this way must be guarded by a sentry object.
  // The sentry object performs various tasks,
  // such as thread synchronization and updating the stream state.
  
  std::istream::sentry se(is, true);
  std::streambuf* sb = is.rdbuf();
  
  for(;;) {
    int c = sb->sbumpc();
    switch (c) {
      case '\n':
        return is;
      case '\r':
        if(sb->sgetc() == '\n')
          sb->sbumpc();
        return is;
      case EOF:
        // Also handle the case when the last line has no line ending
        if(t.empty())
          is.setstate(std::ios::eofbit);
        return is;
      default:
        t += (char)c;
    }
  }
}

std::ostream& operator<<(std::ostream& out, const DoubleMatrix& M)
{
  int m = M.size();
  int n = M.empty() ? -1 : M[0].size();
  
  out << m << std::endl;
  out << n << std::endl;
  
  for (int i = 0; i < m; ++i)
  {
    for (int k = 0; k < n; ++k)
    {
      out << M[i][k] << " ";
    }
    out << std::endl;
  }
  
  return out;
}

std::istream& operator>>(std::istream& in, DoubleMatrix& M)
{
  int m = -1, n = -1;
  
  std::string line;
  getline(in, line);
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
  
  M = DoubleMatrix(m, DoubleVector(n, 0));
  for (int i = 0; i < m; ++i)
  {
    getline(in, line);
    ss.clear();
    ss.str(line);
    
    for (int j = 0; j < n; ++j)
    {
      ss >> M[i][j];
    }
  }
  
  return in;
}

std::ostream& operator<<(std::ostream& out, const DoubleTensor& M)
{
  int k = M.size();
  int m = k == 0 ? 0 : M.front().size();
  int n = m == 0 ? 0 : M.front().front().size();
  
  out << k << " #k" << std::endl;
  out << m << " #m" << std::endl;
  out << n << " #n" << std::endl;
  
  for (int i = 0; i < k; ++i)
  {
    for (int p = 0; p < m; ++p)
    {
      bool first = true;
      for (int c = 0; c < n; ++c)
      {
        if (first)
          first = false;
        else
          out << " ";
        
        out << M[i][p][c];
      }
      out << std::endl;
    }
  }
  
  return out;
}

std::istream& operator>>(std::istream& in, DoubleTensor& M)
{
  int k = -1, m = -1, n = -1;
  
  std::string line;
  getline(in, line);
  std::stringstream ss(line);
  
  ss >> k;
  if (k <= 0)
  {
    throw std::runtime_error("Error: k should be nonnegative");
  }
  
  getline(in, line);
  ss.clear();
  ss.str(line);
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
  
  M = DoubleTensor(k, DoubleMatrix(m, DoubleVector(n, 0)));
  for (int i = 0; i < k; ++i)
  {
    for (int p = 0; p < m; ++p)
    {
      getline(in, line);
      ss.clear();
      ss.str(line);
      
      for (int c = 0; c < n; ++c)
      {
        ss >> M[i][p][c];
      }
    }
  }
  
  return in;
}
