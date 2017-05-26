/*
 * pmh_pr.cpp
 *
 *  Created on: 2-feb-2017
 *      Author: M. El-Kebir
 */

#include <iostream>
#include "utils.h"
#include "nonbinaryclonetree.h"
#include <fstream>
#include <lemon/arg_parser.h>
#include "ilpsolver.h"
#include "ilpbinarizationsolver.h"
#include "migrationgraph.h"
#include <boost/algorithm/string.hpp>

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
  int pattern = -1;
  int nrThreads = -1;
  int timeLimit = -1;
  double UB = -1;
  std::string sampleTreeFile;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("c", "Color map file", filenameColorMap)
    .other("T", "Clone tree")
    .other("leaf_labeling", "Leaf labeling")
    .refOption("g", "Output search graph", outputSearchGraph)
    .refOption("log", "Gurobi logging", gurobiLog)
    .refOption("t", "Number of threads (default: -1)", nrThreads)
    .refOption("o", "Output prefix" , outputDirectory)
    .refOption("m", "Migration pattern:\n"\
                    "       0 : Single-site parallel seeding"\
                    "       1 : Single-site seeding\n" \
                    "       2 : Multi-site seeding\n" \
                    "       3 : Reseeding", pattern)
    .refOption("s", "Fix comigrations according to provided sample tree", sampleTreeFile)
    .refOption("e", "Export ILP", outputILP)
    .refOption("p", "Primary", primary)
    .refOption("UB", "Upper bound (default: -1, disabled)", UB)
    .refOption("l", "Time limit in seconds (default: -1)", timeLimit);
  ap.parse();
  
  if (ap.files().size() != 2)
  {
    std::cerr << "Error: <T> and <leaf_labeling> must be specified" << std::endl;
    return 1;
  }
  
  StringVector primaryVector;
  if (!primary.empty())
  {
    boost::split(primaryVector, primary, boost::is_any_of(","));
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
      std::cerr << "Could not open '" << filenameColorMap << "' for reading"
        << std::endl;
      return 1;
    }
    
    if (!BaseTree::readColorMap(inColorMap, colorMap))
    {
      return 1;
    }
  }
  
  std::cerr << "Clone tree has " << T.getNrSamples() << " samples" << std::endl;
  
  StringPairList forcedComigrations;
  if (!sampleTreeFile.empty()
      && !parseMigrationGraph(sampleTreeFile, T.getSamples(), forcedComigrations))
  {
    return 1;
  }

  std::list<MigrationGraph::Pattern> patterns;
  if (pattern != -1)
  {
    patterns.push_back(static_cast<MigrationGraph::Pattern>(pattern));
  }
  else
  {
    for (int i = 0; i < MigrationGraph::_nrPatterns; ++i)
    {
      patterns.push_back(static_cast<MigrationGraph::Pattern>(i));
    }
  }
  
  StringSet sampleSet = T.getSamples();
  if (primaryVector.empty())
  {
    std::cerr << "Warning: Considering all samples as primary." << std::endl;
    primaryVector = StringVector(sampleSet.begin(),
                                 sampleSet.end());
  }

  for (MigrationGraph::Pattern pattern : patterns)
  {
    for (const std::string& primary : primaryVector)
    {
      if (sampleSet.count(primary) == 0)
      {
        std::cerr << "Warning: primary sample '" << primary
          << "' missing in leaf labeling. Skipping." << std::endl;
        continue;
      }
      IlpBinarizationSolver::run(T,
                                 primary,
                                 outputDirectory,
                                 colorMap,
                                 static_cast<MigrationGraph::Pattern>(pattern) ,
                                 nrThreads,
                                 outputILP,
                                 outputSearchGraph,
                                 timeLimit, UB, forcedComigrations);
    }
  }

  return 0;
}
