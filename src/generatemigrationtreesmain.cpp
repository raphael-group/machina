/*
 * generatecondensedsampletrees.cpp
 *
 *  Created on: 18-oct-2016
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "gabowmyers.h"
#include "migrationtree.h"
#include <stdlib.h>
#include <fstream>

int main(int argc, char** argv)
{
  if (argc < 4)
  {
    std::cerr << "Usage: " << argv[0] << " <ANATOMICAL_SITES>" << std::endl;
    return 1;
  }
  
  int nrSamples = argc - 1;
  if (nrSamples < 2)
  {
    std::cerr << "Error: #ANATOMICAL_SITES must be at least 2" << std::endl;
    return 1;
  }
  
  StringSet mets;
  for (int i = 0; i < nrSamples - 1; ++i)
  {
    mets.insert(argv[i + 2]);
  }
  
  MigrationTreeList listBarS;
  MigrationTree::enumerate(argv[1], mets, listBarS);
  
  std::cout << listBarS;
  
  return 0;
}
