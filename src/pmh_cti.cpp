/*
 * pmh_cti.cpp
 *
 *  Created on: 13-apr-2017
 *      Author: M. El-Kebir
 */

#include <iostream>
#include "utils.h"
#include "nonbinaryclonetree.h"
#include "frequencymatrix.h"
#include <fstream>
#include <lemon/arg_parser.h>
#include "ilpsolverext.h"
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
  ap.refOption("c", "Color map file", filenameColorMap, true)
    .other("T", "Mutation tree")
    .other("frequencies", "Frequencies")
    .refOption("g", "Output search graph", outputSearchGraph)
    .refOption("log", "Gurobi logging", gurobiLog)
    .refOption("t", "Number of threads (default: -1, #cores)", nrThreads)
    .refOption("o", "Output prefix" , outputDirectory)
    .refOption("m", "Allowed migration patterns:\n"\
               "       0 : PS\n"\
               "       1 : PS, S\n" \
               "       2 : PS, S, M\n" \
               "       3 : PS, S, M, R\n" \
               "     If no pattern is specified, all allowed patterns will be enumerated.", pattern)
    .refOption("s", "Fix comigrations according to provided migration graph", sampleTreeFile)
    .refOption("e", "Export ILP", outputILP)
    .refOption("p", "Primary samples separated by commas (if omitted, every sample will be\n" \
               "     considered iteratively as the primary)", primary)
    .refOption("UB", "Upper bound (default: -1, disabled)", UB)
    .refOption("l", "Time limit in seconds (default: -1)", timeLimit);
  ap.parse();
  
  if (ap.files().size() != 2)
  {
    std::cerr << "Error: <T> and <frequencies> must be specified" << std::endl;
    return 1;
  }
  
  StringVector primaryVector;
  if (!primary.empty())
  {
    boost::split(primaryVector, primary, boost::is_any_of(","));
  }
  
  std::string filenameT = ap.files()[0];
  std::string filenameFrequencies = ap.files()[1];
  
  std::ifstream inT(filenameT.c_str());
  if (!inT.good())
  {
    std::cerr << "Could not open '" << filenameT << "' for reading" << std::endl;
    return 1;
  }

  NonBinaryCloneTree T;
  if (!T.read(inT)) return 1;
  inT.close();
  
  std::ifstream inFrequencies(filenameFrequencies.c_str());
  if (!inFrequencies.good())
  {
    std::cerr << "Could not open '" << filenameFrequencies << "' for reading" << std::endl;
    return 1;
  }
  
  FrequencyMatrix F;
  inFrequencies >> F;
  inFrequencies.close();
  
  StringToIntMap colorMap;
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
  
//  std::cerr << "Clone tree has " << F.getNrSamples() << " samples" << std::endl;
  
  StringPairList forcedComigrations;
  if (!sampleTreeFile.empty()
      && !parseMigrationGraph(sampleTreeFile, F.getSamples(), forcedComigrations))
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
  
  StringSet sampleSet = F.getSamples();
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
      IlpSolverExt::run(T,
                        F,
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
