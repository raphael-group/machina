/*
 * clustermain.cpp
 *
 *  Created on: 7-sep-2017
 *      Author: M. El-Kebir
 */

#include "cluster.h"
#include "utils.h"
#include <lemon/arg_parser.h>
#include <fstream>
#include "cluster/readmatrix.h"
#include "cluster/cluster.h"

int main(int argc, char** argv)
{
  double alpha = 0.001;
  double beta = 0.01;
  bool relabel = false;
  bool outputAncesTree = false;
  std::string clusteringFilename;
  double fwr = -1;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("a", "Confidence interval used for clustering (default: 0.001)", alpha)
    .refOption("FWR", "Family-wise error rate", fwr)
    .refOption("b", "Confidence interval used for pooled frequency matrix (default: 0.01)", beta)
    .refOption("A", "Output AncesTree input file", outputAncesTree)
    .refOption("C", "Clustering input filename", clusteringFilename)
    .refOption("r", "Relabel mutation clusters", relabel)
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
  
  if (fwr != -1)
  {
    alpha = fwr * R.getNrSamples() * R.getNrCharacters();
  }
  
  Cluster cluster(R, alpha, relabel);
  if (clusteringFilename.empty())
  {
    cluster.clusterCC(beta);
  }
  else
  {
    std::ifstream inC(clusteringFilename.c_str());
    if (!inC.good())
    {
      std::cerr << "Error: could not open '" << clusteringFilename << "' for reading" << std::endl;;
      return 1;
    }
    cluster.readClustering(inC, beta);
  }
  cluster.writeClustering(std::cerr);
  
  if (outputAncesTree)
  {
    cluster.writeAncesTreeInput(std::cout);
  }
  else
  {
    std::cout << cluster.getClusteredF();
  }
  
  return 0;
}
