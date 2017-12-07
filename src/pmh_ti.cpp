/*
 * pmh_ti.cpp
 *
 *  Created on: 13-apr-2017
 *      Author: M. El-Kebir
 */

#include <iostream>
#include "utils.h"
#include "clonetree.h"
#include "frequencymatrix.h"
#include <fstream>
#include <lemon/arg_parser.h>
#include "old_ilps/ilpsolverext.h"
#include "ilppmhtisolver.h"
#include "migrationgraph.h"
#include <boost/algorithm/string.hpp>
#include "spruce/rootedcladisticnoisyenumeration.h"
#include "spruce/statetree.h"
#include "spruce/perfectphylotree.h"
#include "enumeratemutationtrees.h"

int main(int argc, char** argv)
{
  std::string filenameColorMap;
  std::string filenameVertexLabeling;
  std::string filenameSearchGraph;
  std::string filenameBarT;
  bool gurobiLog = false;
  bool outputILP = false;
  bool outputSearchGraph = false;
  std::string filenameFrequencies;
  std::string filenameOutMigrationGraph;
  std::string filenameOutCloneTree;
  std::string outputDirectory;
  std::string primary;
  std::string pattern = "0,1,2,3";
  int nrThreads = -1;
  int timeLimitILP = -1;
  int seed = 0;
  IntTriple bounds = std::make_pair(-1, std::make_pair(-1, -1));
  std::string migrationTreeFile;
  int nrMutTrees = -1;
  bool oldMode = false;
  bool disablePolytomyResolution = false;
  int mutationTreeIdx = -1;
  bool useBounds = false;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("c", "Color map file", filenameColorMap, true)
    .refOption("barT", "Mutation trees", filenameBarT, true)
    .refOption("F", "Frequencies file", filenameFrequencies, true)
    .refOption("g", "Output search graph", outputSearchGraph)
    .refOption("log", "Gurobi logging", gurobiLog)
    .refOption("t", "Number of threads (default: -1, #cores)", nrThreads)
    .refOption("o", "Output prefix" , outputDirectory)
    .refOption("noPR", "Disable polytomy resolution", disablePolytomyResolution)
    .refOption("useBounds", "Only retain optimal solution", useBounds)
    .refOption("OLD", "Use old ILP (typically much slower)", oldMode)
    .refOption("mutTreeIdx", "Mutation tree index (default: -1)", mutationTreeIdx)
    .refOption("m", "Allowed migration patterns:\n"\
               "       0 : PS\n"\
               "       1 : PS, S\n" \
               "       2 : PS, S, M\n" \
               "       3 : PS, S, M, R\n" \
               "     If no pattern is specified, all allowed patterns will be\n"
               "     enumerated (default: '0,1,2,3')", pattern)
    .refOption("G", "Optional file with migration graphs", migrationTreeFile)
    .refOption("e", "Export ILP", outputILP)
    .refOption("p", "Primary anatomical site", primary, true)
    .refOption("UB_mu", "Upper bound on the migration number (default: -1, disabled)", bounds.first)
    .refOption("UB_gamma", "Upper bound on the comigration number (default: -1, disabled)", bounds.second.first)
    .refOption("UB_sigma", "Upper bound on the seeding site number (default: -1, disabled)", bounds.second.second)
    .refOption("l", "Time limit in seconds for the ILP (default: -1, unlimited)", timeLimitILP);
  ap.parse();
  
  FrequencyMatrix F;
  try
  {
    std::ifstream inFrequencies(filenameFrequencies.c_str());
    if (!inFrequencies.good())
    {
      std::cerr << "Could not open '" << filenameFrequencies << "' for reading" << std::endl;
      return 1;
    }
    
    inFrequencies >> F;
    inFrequencies.close();
  }
  catch (std::runtime_error& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  
  if (!F.isAnatomicalSite(primary))
  {
    std::cerr << "Error: " << primary << " is not an anatomical site" << std::endl;
    return 1;
  }
  
  // Mutation tree(s)
  EnumerateMutationTrees::TreeVector mutationTrees;
  try
  {
    std::ifstream inT(filenameBarT.c_str());
    if (!inT.good())
    {
      std::cerr << "Could not open '" << filenameBarT << "' for reading" << std::endl;
      return 1;
    }
    
    inT >> mutationTrees;
  }
  catch (std::runtime_error& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  std::cerr << "Read " << mutationTrees.size() << " mutation trees." << std::endl << std::endl;
  
  if (nrMutTrees > 0)
  {
    g_rng = std::mt19937(seed);
    EnumerateMutationTrees::pick(mutationTrees, nrMutTrees);
    std::cerr << "Picked " << mutationTrees.size() << " mutation trees with seed " << seed << std::endl << std::endl;
  }
  
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
        int mutTreeIdx = 0;
        for (const CloneTree& barT : mutationTrees)
        {
          if (mutationTreeIdx == -1 || mutTreeIdx == mutationTreeIdx)
          {
            snprintf(buf, 1024, "%d-%d-", migTreeIdx, mutTreeIdx);
            
            if (!oldMode)
            {
              IntTriple res = IlpPmhTiSolver::run(barT,
                                                  F,
                                                  primary,
                                                  outputDirectory,
                                                  buf,
                                                  colorMap,
                                                  pattern,
                                                  nrThreads,
                                                  outputILP,
                                                  outputSearchGraph,
                                                  timeLimitILP,
                                                  bounds,
                                                  migrationTree,
                                                  disablePolytomyResolution);
              
              if (res.first != -1 && useBounds)
              {
                bounds = res;
              }
            }
            else
            {
              IlpSolverExt::run(barT,
                                F,
                                primary,
                                outputDirectory,
                                buf,
                                colorMap,
                                pattern,
                                nrThreads,
                                outputILP,
                                outputSearchGraph,
                                timeLimitILP,
                                bounds,
                                migrationTree);
            }
          }
          ++mutTreeIdx;
        }
        ++migTreeIdx;
      }
    }
    else
    {
      int mutTreeIdx = 0;
      for (const CloneTree& barT : mutationTrees)
      {
        if (mutationTreeIdx == -1 || mutTreeIdx == mutationTreeIdx)
        {
          snprintf(buf, 1024, "%d-", mutTreeIdx);
          if (!oldMode)
          {
            IntTriple res = IlpPmhTiSolver::run(barT,
                                                F,
                                                primary,
                                                outputDirectory,
                                                buf,
                                                colorMap,
                                                pattern,
                                                nrThreads,
                                                outputILP,
                                                outputSearchGraph,
                                                timeLimitILP,
                                                bounds,
                                                StringPairList(),
                                                disablePolytomyResolution);
            
            if (res.first != -1 && useBounds)
            {
              bounds = res;
            }
          }
          else
          {
            IlpSolverExt::run(barT,
                              F,
                              primary,
                              outputDirectory,
                              buf,
                              colorMap,
                              pattern,
                              nrThreads,
                              outputILP,
                              outputSearchGraph,
                              timeLimitILP,
                              bounds,
                              StringPairList());
          }
        }
        ++mutTreeIdx;
      }
    }
  }

  return 0;
}
