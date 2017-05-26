/*
 * solvemain.cpp
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
#include <lemon/time_measure.h>

std::string solverModeFullString(IlpSolver::Mode mode)
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

std::string solverModeString(IlpSolver::Mode mode)
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
               const std::string& primary,
               const std::string& outputPrefix,
               const StringToIntMap& colorMap,
               IlpSolver::Mode mode,
               bool binarization,
               int nrThreads,
               bool outputILP,
               bool outputSearchGraph,
               int timelimit,
               double UB,
               const StringPairList& forcedComigrations)
{
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
  
  if (!binarization)
  {
    IlpSolver solver(T,
                     primary,
                     mode,
                     filenameGurobiLog,
                     forcedComigrations);
    solver.init(UB);
    
    if (!outputPrefix.empty() && outputILP)
    {
      snprintf(buf, 1024, "%s/ilp-%s-%s%s.ilp",
               outputPrefix.c_str(),
               primary.c_str(),
               solverModeString(mode).c_str(),
               binarizationString(binarization).c_str());
      solver.exportModel(buf);
    }
    
    lemon::Timer timer;
    if (!solver.solve(nrThreads, timelimit))
    {
      std::cerr << "No solution found" << std::endl;
      return;
    }
    MigrationGraph G(T, solver.lPlus());
    
    std::cerr << "With primary '" << primary << "', " << solverModeFullString(mode) << " and " << binarizationFullString(binarization) << ": "
      << G.getNrMigrations() << " migrations, "
      << G.getNrComigrations(T, solver.lPlus()) << " comigrations, "
      << G.getNrNonUniqueParentageSamples() << " non-unique parentage sites and "
      << G.getNrSeedingSamples() << " seeding sites";
    if (G.hasReseeding())
    {
      std::cerr << " including reseeding";
    }
    std::cerr << ". [LB, UB] = [" << solver.LB() << ", " << solver.UB() << "]. " << timer.realTime() << " seconds" << std::endl;
    
    if (!outputPrefix.empty())
    {
      snprintf(buf, 1024, "%s/T-%s-%s%s.dot",
               outputPrefix.c_str(),
               primary.c_str(),
               solverModeString(mode).c_str(),
               binarizationString(binarization).c_str());
      std::ofstream outT(buf);
      T.writeDOT(outT, solver.lPlus(), colorMap);
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
    NonBinaryCloneTree TT = T;
    TT.mergeSameSampleSiblingLeaves();
    IlpBinarizationSolver solver(TT,
                                 primary,
                                 mode,
                                 filenameGurobiLog,
                                 forcedComigrations);
    solver.init(UB);
    
    if (outputSearchGraph)
    {
      std::ofstream outSearchG(filenameSearchGraph.c_str());
      solver.writeDOT(outSearchG, colorMap);
      outSearchG.close();
    }
    
    if (!outputPrefix.empty() && outputILP)
    {
      snprintf(buf, 1024, "%s/ilp-%s-%s%s.ilp",
               outputPrefix.c_str(),
               primary.c_str(),
               solverModeString(mode).c_str(),
               binarizationString(binarization).c_str());
      solver.exportModel(buf);
    }
    
    lemon::Timer timer;
    if (!solver.solve(nrThreads, timelimit))
    {
      std::cerr << "No solution found" << std::endl;
      return;
    }
    MigrationGraph G(solver.getCloneTree(), solver.lPlus());
    
    std::cerr << "With primary '" << primary << "', " << solverModeFullString(mode) << " and " << binarizationFullString(binarization) << ": "
      << G.getNrMigrations() << " migrations, "
      << G.getNrComigrations(solver.getCloneTree(), solver.lPlus()) << " comigrations, "
      << G.getNrNonUniqueParentageSamples() << " non-unique parentage sites and "
      << G.getNrSeedingSamples() << " seeding sites";
    if (G.hasReseeding())
    {
      std::cerr << " including reseeding";
    }
    std::cerr << ". [LB, UB] = [" << solver.LB() << ", " << solver.UB() << "]. " << timer.realTime() << " seconds" << std::endl;
    
    if (!outputPrefix.empty())
    {
      snprintf(buf, 1024, "%s/T-%s-%s%s.dot",
               outputPrefix.c_str(),
               primary.c_str(),
               solverModeString(mode).c_str(),
               binarizationString(binarization).c_str());
      std::ofstream outT(buf);
      solver.getCloneTree().writeDOT(outT, solver.lPlus(), colorMap);
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
}

bool parseSampleTree(const std::string& sampleTreeFile,
                     const StringSet& samples,
                     StringPairList& forcedComigrations)
{
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
  ap.refOption("c", "Color map file", filenameColorMap)
    .other("T", "Clone tree")
    .other("leaf_labeling", "Leaf labeling")
    .refOption("g", "Output search graph", outputSearchGraph)
    .refOption("log", "Gurobi logging", gurobiLog)
    .refOption("t", "Number of threads (default: -1)", nrThreads)
    .refOption("b", "Binarization", binarization)
    .refOption("o", "Output prefix" , outputPrefix)
    .refOption("m", "Migration pattern:\n"\
                    "       0 : Single-site parallel seeding"\
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
  
  std::cerr << "Clone tree has " << T.getNrSamples() << " samples" << std::endl;
  
  StringPairList forcedComigrations;
  if (!sampleTreeFile.empty()
      && !parseSampleTree(sampleTreeFile, T.getSamples(), forcedComigrations))
  {
    return 1;
  }
//  if (!forcedComigrations.empty())
//  {
//    mode = 1;
//  }

  if (mode != -1)
  {
    if (primaryVector.empty())
    {
      for (const std::string& primary : T.getSamples())
      {
        runSolver(T,
                  primary,
                  outputPrefix,
                  colorMap,
                  static_cast<IlpSolver::Mode>(mode),
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
        if (T.getSamples().count(primary) != 1)
        {
          std::cerr << "Warning: primary sample '" << primary << "' missing in leaf labeling. Skipping." << std::endl;
        }
        else
        {
          runSolver(T,
                    primary,
                    outputPrefix,
                    colorMap,
                    static_cast<IlpSolver::Mode>(mode),
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
                    primary,
                    outputPrefix,
                    colorMap,
                    static_cast<IlpSolver::Mode>(mode),
                    false,
                    nrThreads,
                    outputILP,
                    outputSearchGraph,
                    timeLimit, UB, forcedComigrations);
          runSolver(T,
                    primary,
                    outputPrefix,
                    colorMap,
                    static_cast<IlpSolver::Mode>(mode),
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
          if (T.getSamples().count(primary) != 1)
          {
            std::cerr << "Warning: primary sample '" << primary << "' missing in leaf labeling. Skipping." << std::endl;
          }
          else
          {
            runSolver(T,
                      primary,
                      outputPrefix,
                      colorMap,
                      static_cast<IlpSolver::Mode>(mode),
                      false,
                      nrThreads,
                      outputILP,
                      outputSearchGraph,
                      timeLimit, UB, forcedComigrations);
            runSolver(T,
                      primary,
                      outputPrefix,
                      colorMap,
                      static_cast<IlpSolver::Mode>(mode),
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
