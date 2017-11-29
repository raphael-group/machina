/*
 * visualizeclonetreemain.cpp
 *
 *  Created on: 29-nov-2017
 *      Author: M. El-Kebir
 */

#include <iostream>
#include "utils.h"
#include "clonetree.h"
#include "migrationgraph.h"
#include <fstream>
#include <lemon/arg_parser.h>

int main(int argc, char** argv)
{
  if (argc != 3)
  {
    std::cerr << "Usage: " << argv[0] << " <T> <vertex_labeling>" << std::endl;
    return 1;
  }
  
  std::ifstream inT(argv[1]);
  if (!inT.good())
  {
    std::cerr << "Could not open '" << argv[1] << "' for reading" << std::endl;
    return 1;
  }
  
  std::ifstream inLeafLabeling(argv[2]);
  if (!inLeafLabeling.good())
  {
    std::cerr << "Could not open '" << argv[2] << "' for reading" << std::endl;
    return 1;
  }
  
  CloneTree T;
  try
  {
    if (!T.read(inT)) return 1;
    if (!T.readLeafLabeling(inLeafLabeling)) return 1;
  }
  catch (std::runtime_error& e)
  {
    std::cerr << e.what() << std::endl;
    return 1;
  }
  inT.close();
  inLeafLabeling.close();
  
  StringNodeMap lPlus(T.tree());
  std::ifstream inVertexLabeling(argv[2]);
  if (!inVertexLabeling.good())
  {
    std::cerr << "Could not open '" << argv[2] << "' for reading" << std::endl;
    return 1;
  }
  if (!T.readVertexLabeling(inVertexLabeling, T, lPlus))
  {
    return 1;
  }
  inVertexLabeling.close();
  
  MigrationGraph G(T, lPlus);
  
  int mu = G.getNrMigrations();
  int gamma = G.getNrComigrations(T, lPlus);
  int sigma = G.getNrSeedingSites();
  
  std::cout << mu << "\t"
    << gamma << "\t"
    << sigma << "\t"
    << G.getPatternString(G.getPattern(), G.isMonoclonal())
    << std::endl;
  
  return 0;
}
