/*
 * rfmain.cpp
 *
 *  Created on: 31-aug-2017
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "clonetree.h"
#include <fstream>

bool readCloneTree(const std::string& filenameT,
                   const std::string& filenameLeafLabelingT,
                   CloneTree& T)
{
  std::ifstream inT(filenameT.c_str());
  if (!inT.good())
  {
    std::cerr << "Could not open '" << filenameT << "' for reading" << std::endl;
    return false;
  }
  
  std::ifstream inLeafLabeling(filenameLeafLabelingT.c_str());
  if (!inLeafLabeling.good())
  {
    std::cerr << "Could not open '" << filenameLeafLabelingT << "' for reading" << std::endl;
    return false;
  }
  
  if (!T.read(inT)) return false;
  if (!T.readLeafLabeling(inLeafLabeling)) return false;
  
  return true;
}

int main(int argc, char** argv)
{
  if (argc != 5 && argc != 3)
  {
    std::cerr << "Usage: " << argv[0] << " <T1> <lT1> (<T2> <lT2>)" << std::endl;
    return 1;
  }
  
  std::string filenameT1 = argv[1];
  std::string filenameLeafLabelingT1 = argv[2];

  std::string filenameT2 = argc == 5 ? argv[3] : "";
  std::string filenameLeafLabelingT2 = argc == 5 ? argv[4] : "";
  
  CloneTree T1, T2;
  if (!readCloneTree(filenameT1, filenameLeafLabelingT1, T1))
    return 1;
  if (argc == 5 && !readCloneTree(filenameT2, filenameLeafLabelingT2, T2))
    return 1;
  
  SplitSet splitSet1 = T1.getSplits();
  std::cout << "Splits in T1:" << std::endl;
  for (const Split& split : splitSet1)
  {
    bool first = true;
    for (const StringSet& X : split)
    {
      if (first)
        first = false;
      else
        std::cout << " ;";

      for (const std::string& l : X)
      {
        std::cout << " " << l;
      }
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  
  SplitSet splitSet2 = T2.getSplits();
  std::cout << "Splits in T2:" << std::endl;
  for (const Split& split : splitSet2)
  {
    bool first = true;
    for (const StringSet& X : split)
    {
      if (first)
        first = false;
      else
        std::cout << " ;";

      for (const std::string& l : X)
      {
        std::cout << " " << l;
      }
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  
  SplitSet result;
  std::set_symmetric_difference(splitSet1.begin(), splitSet1.end(),
                                splitSet2.begin(), splitSet2.end(),
                                std::inserter(result, result.begin()));
  
  std::cout << "Symmetric difference:" << std::endl;
  for (const Split& split : result)
  {
    if (splitSet1.count(split) == 1)
      std::cout << "T1:";
    if (splitSet2.count(split) == 1)
      std::cout << "T2:";
    
    bool first = true;
    for (const StringSet& X : split)
    {
      if (first)
        first = false;
      else
        std::cout << " ;";
      for (const std::string& l : X)
      {
        std::cout << " " << l;
      }
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  
  
  std::cout << "Robinson-Foulds distance: " << result.size() << std::endl;
  
  return 0;
}
