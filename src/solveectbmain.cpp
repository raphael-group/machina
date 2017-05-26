/*
 * solveectbmain.cpp
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
#include "ilpsolver.h"
#include "ilpbinarizationsolverext.h"
#include "ilpsolverext.h"
#include "migrationgraph.h"
#include <boost/algorithm/string.hpp>
#include <lemon/time_measure.h>

std::string solverModeFullString(IlpSolverExt::Mode mode)
{
  switch (mode)
  {
    case IlpSolver::SINGLE_SOURCE_SEEDING:
      return "single-site seeding";
    case IlpSolver::MULTI_SOURCE_SEEDING:
      return "multi-site seeding";
    case IlpSolver::RESEEDING:
      return "reseeding";
    case IlpSolver::PARALLEL_SINGLE_SOURCE_SEEDING:
      return "parallel single-site seeding";
    default:
      assert(false);
      return "ERROR";
  }
}

std::string solverModeString(IlpSolverExt::Mode mode)
{
  switch (mode)
  {
    case IlpSolver::SINGLE_SOURCE_SEEDING:
      return "SSS";
    case IlpSolver::MULTI_SOURCE_SEEDING:
      return "MSS";
    case IlpSolver::RESEEDING:
      return "RS";
    case IlpSolver::PARALLEL_SINGLE_SOURCE_SEEDING:
      return "PSSS";
    default:
      assert(false);
      return "ERR";
  }
}

std::string binarizationFullString(bool binarization)
{
  if (binarization)
  {
    return "binarization";
  }
  else
  {
    return "no binarization";
  }
}

std::string binarizationString(bool binarization)
{
  if (binarization)
  {
    return "-binarized";
  }
  else
  {
    return "";
  }
}

void runSolver(const NonBinaryCloneTree& T,
               const FrequencyMatrix& F,
               const std::string& primary,
               const std::string& outputPrefix,
               const StringToIntMap& colorMap,
               IlpSolverExt::Mode mode,
               bool binarization,
               int nrThreads,
               bool outputILP,
               bool outputSearchGraph,
               int timelimit,
               double UB,
               const StringPairList& forcedComigrations)
{
  StringSet mets = F.getSamples();
  mets.erase(primary);

//  MigrationTreeList listBarS;
//  MigrationTree::enumerate(primary, mets, listBarS);
//  
//  int index = 0;
//  for (MigrationTree& barS : listBarS)
//  {
//    ++index;
  
    char buf[1024];
    std::string filenameGurobiLog;
    std::string filenameSearchGraph;
    
    if (!outputPrefix.empty())
    {
      snprintf(buf, 1024, "%s/log-%s-%s%s.txt",
               outputPrefix.c_str(),
               primary.c_str(),
               solverModeString(mode).c_str(),
               binarizationString(binarization).c_str());
      
      filenameGurobiLog = buf;
      
      snprintf(buf, 1024, "%s/searchG-%s-%s%s.dot",
               outputPrefix.c_str(),
               primary.c_str(),
               solverModeString(mode).c_str(),
               binarizationString(binarization).c_str());
      
      filenameSearchGraph = buf;
    }
    
//    StringPairList forcedComigrations = barS.getEdgeList();
    if (!binarization)
    {
      IlpSolverExt solver(T,
                          F,
                          primary,
                          mode,
                          filenameGurobiLog,
                          forcedComigrations);
      solver.init(UB);
      
      if (!outputPrefix.empty() && outputILP)
      {
        snprintf(buf, 1024, "%s/ilp-%s-%s%s.lp",
                 outputPrefix.c_str(),
                 primary.c_str(),
                 solverModeString(mode).c_str(),
                 binarizationString(binarization).c_str());
        solver.exportModel(buf);
      }
      
      if (outputSearchGraph)
      {
        std::ofstream outSearchG(filenameSearchGraph.c_str());
        solver.writeDOT(outSearchG, colorMap);
        outSearchG.close();
      }
      
      lemon::Timer timer;
      if (!solver.solve(nrThreads, timelimit))
      {
        std::cerr << "No solution found (" << outputPrefix << ")" << std::endl;
        return;
      }

      MigrationGraph G(solver.T(), solver.lPlus());
      
      std::cerr << "With primary '" << primary << "', "
        << solverModeFullString(mode) << " and " << binarizationFullString(binarization) << ": "
        << G.getNrMigrations() << " migrations, "
        << G.getNrComigrations(solver.T(), solver.lPlus()) << " comigrations, "
        << G.getNrNonUniqueParentageSamples() << " non-unique parentage sites and "
        << G.getNrSeedingSamples() << " seeding sites";
      if (G.hasReseeding())
      {
        std::cerr << " including reseeding";
      }
      std::cerr << ". [LB, UB] = [" << solver.LB() << ", " << solver.UB() << "]. " << timer.realTime() << " seconds (" << outputPrefix << ")" << std::endl;
      
      if (!outputPrefix.empty())
      {
        snprintf(buf, 1024, "%s/T-%s-%s%s.dot",
                 outputPrefix.c_str(),
                 primary.c_str(),
                 solverModeString(mode).c_str(),
                 binarizationString(binarization).c_str());
        std::ofstream outT(buf);
        solver.T().writeDOT(outT,
                            solver.lPlus(),
                            colorMap,
                            solver.getF(),
                            solver.getU());
        outT.close();
        
        snprintf(buf, 1024, "%s/G-%s-%s%s.dot",
                 outputPrefix.c_str(),
                 primary.c_str(),
                 solverModeString(mode).c_str(),
                 binarizationString(binarization).c_str());
        
        std::ofstream outG(buf);
        G.writeDOT(outG, colorMap);
        outG.close();
      }
    }
    else
    {
      IlpBinarizationSolverExt solver(T,
                                      F,
                                      primary,
                                      mode,
                                      filenameGurobiLog,
                                      forcedComigrations);
      solver.init(UB);
      
      if (!outputPrefix.empty() && outputILP)
      {
        snprintf(buf, 1024, "%s/ilp-%s-%s%s.lp",
                 outputPrefix.c_str(),
                 primary.c_str(),
                 solverModeString(mode).c_str(),
                 binarizationString(binarization).c_str());
        solver.exportModel(buf);
      }
      
      if (outputSearchGraph)
      {
        std::ofstream outSearchG(filenameSearchGraph.c_str());
        solver.writeDOT(outSearchG, colorMap);
        outSearchG.close();
      }
      
      lemon::Timer timer;
      if (!solver.solve(nrThreads, timelimit))
      {
        std::cerr << "No solution found (" << outputPrefix << ")" << std::endl;
//        continue;
        return;
      }
      
      MigrationGraph G(solver.T(), solver.lPlus());
      
      std::cerr << "With primary '" << primary << "', " << solverModeFullString(mode) << " and " << binarizationFullString(binarization) << ": "
        << G.getNrMigrations() << " migrations, "
        << G.getNrComigrations(solver.T(), solver.lPlus()) << " comigrations, "
        << G.getNrNonUniqueParentageSamples() << " non-unique parentage sites and "
        << G.getNrSeedingSamples() << " seeding sites";
      if (G.hasReseeding())
      {
        std::cerr << " including reseeding";
      }
      std::cerr << ". [LB, UB] = [" << solver.LB() << ", " << solver.UB() << "]. " << timer.realTime() << " seconds (" << outputPrefix << ")" << std::endl;
      
      if (!outputPrefix.empty())
      {
        snprintf(buf, 1024, "%s/T-%s-%s%s.dot",
                 outputPrefix.c_str(),
                 primary.c_str(),
                 solverModeString(mode).c_str(),
                 binarizationString(binarization).c_str());
        std::ofstream outT(buf);
        solver.T().writeDOT(outT,
                            solver.lPlus(),
                            colorMap,
                            solver.getF(),
                            solver.getU());
        outT.close();
        
        snprintf(buf, 1024, "%s/T-%s-%s%s-condensed.dot",
                 outputPrefix.c_str(),
                 primary.c_str(),
                 solverModeString(mode).c_str(),
                 binarizationString(binarization).c_str());
        std::ofstream outCondensedT(buf);
        solver.T().writeDOT(outCondensedT,
                            solver.lPlus(),
                            colorMap,
                            solver.getU(),
                            solver.getCharacterLabel());
        outCondensedT.close();
        
        snprintf(buf, 1024, "%s/G-%s-%s%s.dot",
                 outputPrefix.c_str(),
                 primary.c_str(),
                 solverModeString(mode).c_str(),
                 binarizationString(binarization).c_str());
        
        std::ofstream outG(buf);
        G.writeDOT(outG, colorMap);
        outG.close();
      }
    }
//  }
}

bool parseSampleTree(const std::string& sampleTreeFile,
                     const FrequencyMatrix& F,
                     StringPairList& forcedComigrations)
{
  StringSet samples = F.getSamples();
  
  std::ifstream inFile(sampleTreeFile);
  if (!inFile.good())
  {
    std::cerr << "Could not open '" << sampleTreeFile << "' for reading" << std::endl;
    return false;
  }
  else
  {
    int idx = 0;
    while (inFile.good())
    {
      ++idx;
      std::string line;
      getline(inFile, line);
      
      if (line.empty()) continue;
      
      StringVector s;
      boost::split(s, line, boost::is_any_of("\t "));
      
      if (s.size() < 2)
      {
        std::cerr << "Line " << idx << " is invalid in '" << sampleTreeFile << "'" << std::endl;
        return false;
      }
      if (samples.count(s[0]) != 1)
      {
        std::cerr << "Line " << idx << " is invalid in '" << sampleTreeFile << "'. Sample '" << s[0] << "' is incorrect." << std::endl;
        return false;
      }
      if (samples.count(s[1]) != 1)
      {
        std::cerr << "Line " << idx << " is invalid in '" << sampleTreeFile << "'. Sample '" << s[1] << "' is incorrect." << std::endl;
        return false;
      }
      forcedComigrations.push_back(StringPair(s[0], s[1]));
    }
  }
  
  return true;
}

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
  std::string outputPrefix;
  std::string primary;
  int mode = -1;
  int nrThreads = -1;
  bool binarization = false;
  int timeLimit = -1;
  double UB = -1;
  std::string sampleTreeFile;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("c", "Color map file", filenameColorMap, true)
    .other("T", "Clone tree")
    .other("frequencies", "Frequencies")
    .refOption("g", "Output search graph", outputSearchGraph)
    .refOption("log", "Gurobi logging", gurobiLog)
    .refOption("t", "Number of threads (default: -1)", nrThreads)
    .refOption("b", "Binarization", binarization)
    .refOption("o", "Output prefix" , outputPrefix)
    .refOption("m", "Migration pattern:\n"\
                    "       0 : Single-site parallel seeding\n"\
                    "       1 : Single-site seeding\n" \
                    "       2 : Multi-site seeding\n" \
                    "       3 : Reseeding", mode)
    .refOption("s", "Fix comigrations according to provided sample tree", sampleTreeFile)
    .refOption("e", "Export ILP", outputILP)
    .refOption("p", "Primary", primary)
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
      && !parseSampleTree(sampleTreeFile, F, forcedComigrations))
  {
    return 1;
  }
  if (!forcedComigrations.empty())
  {
    mode = 1;
  }

  if (mode != -1)
  {
    if (primaryVector.empty())
    {
      for (const std::string& primary : T.getSamples())
      {
        runSolver(T,
                  F,
                  primary,
                  outputPrefix,
                  colorMap,
                  static_cast<IlpSolverExt::Mode>(mode),
                  binarization,
                  nrThreads,
                  outputILP,
                  outputSearchGraph,
                  timeLimit, UB, forcedComigrations);
      }
    }
    else
    {
      for (const std::string& primary : primaryVector)
      {
        if (!F.isSample(primary))
        {
          std::cerr << "Warning: primary sample '" << primary << "' missing. Skipping." << std::endl;
        }
        else
        {
          runSolver(T,
                    F,
                    primary,
                    outputPrefix,
                    colorMap,
                    static_cast<IlpSolverExt::Mode>(mode),
                    binarization,
                    nrThreads,
                    outputILP,
                    outputSearchGraph,
                    timeLimit, UB, forcedComigrations);
        }
      }
    }
  }
  else
  {
    if (primaryVector.empty())
    {
      for (const std::string& primary : T.getSamples())
      {
        for (int mode = 0; mode < IlpSolver::_nrModes; ++mode)
        {
          runSolver(T,
                    F,
                    primary,
                    outputPrefix,
                    colorMap,
                    static_cast<IlpSolverExt::Mode>(mode),
                    false,
                    nrThreads,
                    outputILP,
                    outputSearchGraph,
                    timeLimit, UB, forcedComigrations);
          runSolver(T,
                    F,
                    primary,
                    outputPrefix,
                    colorMap,
                    static_cast<IlpSolverExt::Mode>(mode),
                    true,
                    nrThreads,
                    outputILP,
                    outputSearchGraph,
                    timeLimit, UB, forcedComigrations);
        }
      }
    }
    else
    {
      for (const std::string& primary : primaryVector)
      {
        for (int mode = 0; mode < IlpSolver::_nrModes; ++mode)
        {
          if (!F.isSample(primary))
          {
            std::cerr << "Warning: primary sample '" << primary
                      << "' missing. Skipping."
                      << std::endl;
          }
          else
          {
            runSolver(T,
                      F,
                      primary,
                      outputPrefix,
                      colorMap,
                      static_cast<IlpSolverExt::Mode>(mode),
                      false,
                      nrThreads,
                      outputILP,
                      outputSearchGraph,
                      timeLimit, UB, forcedComigrations);
            runSolver(T,
                      F,
                      primary,
                      outputPrefix,
                      colorMap,
                      static_cast<IlpSolverExt::Mode>(mode),
                      true,
                      nrThreads,
                      outputILP,
                      outputSearchGraph,
                      timeLimit, UB, forcedComigrations);
          }
        }
      }
    }
  }

  return 0;
}
