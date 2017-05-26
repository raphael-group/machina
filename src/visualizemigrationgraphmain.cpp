/*
 * visualizemigrationgraphmain.cpp
 *
 *  Created on: 24-oct-2016
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "nonbinaryclonetree.h"
#include "migrationgraph.h"
#include <lemon/arg_parser.h>
#include <fstream>

int main(int argc, char** argv)
{
  std::string filenameColorMap;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("c", "Color map file", filenameColorMap)
    .other("T", "Clone tree")
    .other("leaf_labeling", "Leaf labeling")
    .other("vertex_labeling", "Vertex labeling");
  ap.parse();
  
  if (ap.files().size() != 3)
  {
    std::cerr << "Error: <T>, <leaf_labeling> and <vertex_labeling> must be specified" << std::endl;
    return 1;
  }
  
  std::string filenameT = ap.files()[0];
  std::string filenameLeafLabeling = ap.files()[1];
  std::string filenameVertexLabeling = ap.files()[2];
  
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

  StringNodeMap vertexLabel(T.tree());

  if (filenameVertexLabeling != "-")
  {
    std::ifstream inVertexLabeling(filenameVertexLabeling.c_str());
    if (!inVertexLabeling.good())
    {
      std::cerr << "Could not open '" << filenameVertexLabeling << "' for reading" << std::endl;
      return 1;
    }
    if (!T.readVertexLabeling(inVertexLabeling, T, vertexLabel))
    {
      return 1;
    }
  }
  else
  {
    if (!T.readVertexLabeling(std::cin, T, vertexLabel))
    {
      return 1;
    }
  }
  
  MigrationGraph S(T, vertexLabel);
  
  if (filenameColorMap.empty())
  {
    S.writeDOT(std::cout);
  }
  else
  {
    std::ifstream inColorMap(filenameColorMap.c_str());
    if (!inColorMap.good())
    {
      std::cerr << "Could not open '" << filenameColorMap << "' for reading" << std::endl;
      return 1;
    }
    
    StringToIntMap colorMap;
    if (!BaseTree::readColorMap(inColorMap, colorMap))
    {
      return 1;
    }
    
    S.writeDOT(std::cout, colorMap);
  }
  
  return 0;
}
