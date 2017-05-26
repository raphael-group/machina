/*
 * visualizeclonetreemain.cpp
 *
 *  Created on: 19-oct-2016
 *      Author: M. El-Kebir
 */

#include <iostream>
#include "utils.h"
#include "nonbinaryclonetree.h"
#include <fstream>
#include <lemon/arg_parser.h>

int main(int argc, char** argv)
{
  std::string filenameColorMap;
  std::string filenameVertexLabeling;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("c", "Color map file", filenameColorMap)
    .other("T", "Clone tree")
    .other("leaf_labeling", "Leaf labeling")
    .refOption("l", "Vertex labeling", filenameVertexLabeling);
  ap.parse();
  
  if (ap.files().size() != 2)
  {
    std::cerr << "Error: <T> and <leaf_labeling> must be specified" << std::endl;
    return 1;
  }
  
  std::string filenameT = ap.files()[0];
  std::string filenameLeafLabeling = ap.files()[1];
  
  std::ifstream inT(filenameT.c_str());
  if (!inT.good())
  {
    std::cerr << "Could not open '" << filenameT << "' for reading" << std::endl;
    return 1;
  }
  
  std::ifstream inLeafLabeling(filenameLeafLabeling.c_str());
  if (!inLeafLabeling.good())
  {
    std::cerr << "Could not open '" << filenameLeafLabeling << "' for reading" << std::endl;
    return 1;
  }
  
  NonBinaryCloneTree T;
  if (!T.read(inT)) return 1;
  if (!T.readLeafLabeling(inLeafLabeling)) return 1;
  
  StringToIntMap colorMap;
  if (filenameColorMap.empty())
  {
    colorMap = T.generateColorMap();
  }
  else
  {
    std::ifstream inColorMap(filenameColorMap.c_str());
    if (!inColorMap.good())
    {
      std::cerr << "Could not open '" << filenameColorMap << "' for reading" << std::endl;
      return 1;
    }

    if (!BaseTree::readColorMap(inColorMap, colorMap))
    {
      return 1;
    }
  }
  
  if (!filenameVertexLabeling.empty())
  {
    StringNodeMap lPlus(T.tree());
    if (filenameVertexLabeling != "-")
    {
      std::ifstream inVertexLabeling(filenameVertexLabeling.c_str());
      if (!inVertexLabeling.good())
      {
        std::cerr << "Could not open '" << filenameVertexLabeling << "' for reading" << std::endl;
        return 1;
      }
      if (!T.readVertexLabeling(inVertexLabeling, T, lPlus))
      {
        return 1;
      }
    }
    else
    {
      if (!T.readVertexLabeling(std::cin, T, lPlus))
      {
        return 1;
      }
    }
    
    T.writeDOT(std::cout, lPlus, colorMap);
  }
  else
  {
    T.writeDOT(std::cout, colorMap);
  }
  
  return 0;
}
