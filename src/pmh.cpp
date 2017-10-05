/*
 * pmh.cpp
 *
 *  Created on: 2-feb-2017
 *      Author: M. El-Kebir
 */

#include <iostream>
#include "utils.h"
#include "clonetree.h"
#include <fstream>
#include <lemon/arg_parser.h>
#include "old_ilps/ilpsolver.h"
#include "ilppmhsolver.h"
#include "migrationgraph.h"
#include <boost/algorithm/string.hpp>
#include "migrationtree.h"

int main(int argc, char** argv)
{
  std::string filenameColorMap;
  std::string filenameVertexLabeling;
  std::string filenameSearchGraph;
  bool gurobiLog = false;
  bool outputILP = false;
  bool outputSearchGraph = false;
  std::string filenameOutMigrationGraph;
  std::string filenameOutCloneTree;
  std::string outputDirectory;
  std::string primary;
  std::string pattern = "0,1,2,3";
  int nrThreads = -1;
  int timeLimit = -1;
  IntTriple bounds = std::make_pair(-1, std::make_pair(-1, -1));
  std::string migrationTreeFile;
  bool oldMode = false;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("c", "Color map file", filenameColorMap, true)
    .other("T", "Clone tree")
    .other("leaf_labeling", "Leaf labeling")
    .refOption("g", "Output search graph", outputSearchGraph)
    .refOption("log", "Gurobi logging", gurobiLog)
    .refOption("t", "Number of threads (default: -1, #cores)", nrThreads)
    .refOption("o", "Output prefix" , outputDirectory)
    .refOption("OLD", "Use old ILP (typically much slower)", oldMode)
    .refOption("m", "Allowed migration patterns:\n"\
                    "       0 : PS\n"\
                    "       1 : PS, S\n" \
                    "       2 : PS, S, M\n" \
                    "       3 : PS, S, M, R\n" \
                    "     If no pattern is specified, all allowed patterns will be enumerated.", pattern)
    .refOption("G", "Optional file with migration graphs", migrationTreeFile)
    .refOption("e", "Export ILP", outputILP)
    .refOption("p", "Primary anatomical site", primary, true)
    .refOption("UB_mu", "Upper bound on the migration number (default: -1, disabled)", bounds.first)
    .refOption("UB_gamma", "Upper bound on the comigration number (default: -1, disabled)", bounds.second.first)
    .refOption("UB_sigma", "Upper bound on the seeding site number (default: -1, disabled)", bounds.second.second)
    .refOption("l", "Time limit in seconds (default: -1, no time limit)", timeLimit);
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
    std::cerr << "Could not open '" << filenameLeafLabeling << "' for reading"
      << std::endl;
    return 1;
  }
  
  CloneTree T;
  if (!T.read(inT)) return 1;
  if (!T.readLeafLabeling(inLeafLabeling)) return 1;
  
  // Read color map
  StringToIntMap colorMap;
  try
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
  catch (std::runtime_error& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  
  // Read migration tree (optional)
  EdgeListVector migrationTrees;
  try
  {
    if (!migrationTreeFile.empty())
    {
      std::ifstream inBarS(migrationTreeFile.c_str());
      if (!inBarS.good())
      {
        std::cerr << "Error: failed to open '" << migrationTreeFile << "' for reading" << std::endl;
        return 1;
      }
      
      inBarS >> migrationTrees;
      std::cerr << "Read " << migrationTrees.size() << " migration graphs" << std::endl;
    }
  }
  catch (std::runtime_error& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  
  if (T.getAnatomicalSites().count(primary) == 0)
  {
    std::cerr << "Error: " << primary << " is not an anatomical site" << std::endl;
    return 1;
  }
  
  // Parse migration patterns
  std::list<MigrationGraph::Pattern> patterns;
  try
  {
    StringVector s;
    boost::split(s, pattern, boost::is_any_of(","));
    for (const std::string& patternStr : s)
    {
      int patternInt = boost::lexical_cast<int>(patternStr);
      patterns.push_back(static_cast<MigrationGraph::Pattern>(patternInt));
    }
  }
  catch (std::exception& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  
  char buf[1024];
  for (MigrationGraph::Pattern pattern : patterns)
  {
    if (migrationTrees.size() > 0)
    {
      int migTreeIdx = 0;
      for (const StringPairList& migrationTree : migrationTrees)
      {
        if (!oldMode)
        {
          snprintf(buf, 1024, "%d-%s-", migTreeIdx, primary.c_str());
          IlpPmhSolver::run(T,
                            primary,
                            outputDirectory,
                            buf,
                            colorMap,
                            pattern,
                            nrThreads,
                            outputILP,
                            outputSearchGraph,
                            timeLimit,
                            bounds,
                            migrationTree);
        }
        else
        {
          IlpSolver::run(T,
                         primary,
                         outputDirectory,
                         colorMap,
                         pattern,
                         nrThreads,
                         outputILP,
                         outputSearchGraph,
                         timeLimit,
                         bounds,
                         migrationTree);
        }
        ++migTreeIdx;
      }
    }
    else
    {
      if (!oldMode)
      {
        IlpPmhSolver::run(T,
                          primary,
                          outputDirectory,
                          primary + "-",
                          colorMap,
                          pattern,
                          nrThreads,
                          outputILP,
                          outputSearchGraph,
                          timeLimit,
                          bounds,
                          StringPairList());
      }
      else
      {
        IlpSolver::run(T,
                       primary,
                       outputDirectory,
                       colorMap,
                       pattern,
                       nrThreads,
                       outputILP,
                       outputSearchGraph,
                       timeLimit,
                       bounds,
                       StringPairList());
      }
    }
  }

  return 0;
}
