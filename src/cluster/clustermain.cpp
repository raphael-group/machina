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
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("a", "Confidence interval used for clustering (default: 0.001)", alpha)
    .refOption("b", "Confidence interval used for pooled frequency matrix (default: 0.01)", beta)
    .refOption("A", "Output AncesTree input file", outputAncesTree)
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
  
  Cluster cluster(R, alpha, relabel);
  cluster.clusterCC(beta);
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
