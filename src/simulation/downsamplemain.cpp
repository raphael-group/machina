/*
 * downsample.cpp
 *
 *  Created on: 6-dec-2017
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include <lemon/arg_parser.h>
#include <fstream>
#include "cluster/readmatrix.h"

int main(int argc, char** argv)
{
  int seed = 0;
  double purity = 1.0;
  int nrSamplesPerAnatomicalSite = -1;
  int coverage = -1;
  double seqErrorRate = 0;
  double fractionSNVs = 1;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("s", "Random number generator seed (default: 0)", seed, false)
    .refOption("C", "Target coverage (default: -1 maintains orginal counts)", coverage, false)
    .refOption("P", "Purity (default: 1.0)", purity, false)
    .refOption("k", "Number of samples per anatomical site (default: -1, maintains original samples)", nrSamplesPerAnatomicalSite, false)
    .refOption("E", "Per base sequencing error rate (default: 0)", seqErrorRate, false)
    .refOption("snvFrac", "Fraction of SNVs to retain (default: 1)", fractionSNVs, false)
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
  
  g_rng = std::mt19937(seed);
  
  ReadMatrix newR = R.downSample(nrSamplesPerAnatomicalSite, coverage, purity, seqErrorRate, fractionSNVs);
  
  std::cout << newR;
  
  return 0;
}
