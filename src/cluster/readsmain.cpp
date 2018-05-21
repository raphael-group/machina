/*
 * readsmain.cpp
 *
 *  Created on: 7-sep-2017
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include <lemon/arg_parser.h>
#include <fstream>
#include "cluster/readmatrix.h"

IntMatrix readClustering(std::istream& in)
{
  IntMatrix clustering;
  std::string line;
  while (in.good())
  {
    getline(in, line);
    if (line.empty()) continue;
    StringVector s;
    boost::split(s, line, boost::is_any_of(";"));
    IntVector c;
    for (const std::string& ss : s)
    {
      c.push_back(boost::lexical_cast<int>(ss));
    }
    clustering.push_back(c);
  }
  return clustering;
}

int main(int argc, char** argv)
{
  double alpha = 0.000001;
  std::string clusteringFilename;
  lemon::ArgParser ap(argc, argv);
  ap.refOption("a", "Confidence interval used for clustering (default: 0.000001)", alpha)
    .refOption("C", "Clustering filename", clusteringFilename)
    .other("R", "Read matrix");
  ap.parse();
  
  if (ap.files().size() != 1)
  {
    std::cerr << "Error: <R> must be specificed" << std::endl;
    return 1;
  }
  
  std::string filenameR = ap.files()[0].c_str();
  
  ReadMatrix R;
  try
  {
    if (filenameR != "-")
    {
      std::ifstream inR(filenameR.c_str());
      if (!inR.good())
      {
        std::cerr << "Error: failed to open '" << ap.files()[0] << "' for reading" << std::endl;
        return 1;
      }
      
      inR >> R;
    }
    else
    {
      std::cin >> R;
    }
  }
  catch (std::runtime_error& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  
  if (clusteringFilename == "")
  {
    FrequencyMatrix F = R.toFrequencyMatrix(alpha, 2);
    std::cout << F;
  }
  else
  {
    std::ifstream inC(clusteringFilename);
    if (!inC.good())
    {
      std::cerr << "Error: failed to open '" << clusteringFilename << "' for reading" << std::endl;
      return 1;
    }
    
    IntMatrix C;
    try
    {
      C = readClustering(inC);
    }
    catch (std::exception& e)
    {
      std::cerr << e.what() << std::endl;
      return 1;
    }
    
    std::cout << R.poolReads(C, false).toFrequencyMatrix(alpha, 2);
  }
  
  return 0;
}
