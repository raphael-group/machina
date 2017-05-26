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
    std::cerr << "Usage: " << argv[0] << " <SAMPLES>" << std::endl;
    return 1;
  }
  
  int nrSamples = argc - 1;
  if (nrSamples < 2)
  {
    std::cerr << "Error: #SAMPLES must be at least 2" << std::endl;
    return 1;
  }
  
  StringSet mets;
  for (int i = 0; i < nrSamples - 1; ++i)
  {
    mets.insert(argv[i + 2]);
  }
  
  MigrationTreeList listBarS;
  MigrationTree::enumerate(argv[1], mets, listBarS);
  
  std::cout << nrSamples << " #samples" << std::endl;
  std::cout << listBarS.size() << " #condensed sample trees" << std::endl;
  
  int index = 1;
  char buf[1024];
  for (const MigrationTree& barS : listBarS)
  {
    StringPairList edgeList = barS.getEdgeList();

    snprintf(buf, 1024, "%d.txt", index);
    std::ofstream outFile(buf);
    for (const StringPair& st : edgeList)
    {
      outFile << st.first << " " << st.second << std::endl;
    }
    outFile.close();
    
    ++index;
  }
  
  return 0;
}
