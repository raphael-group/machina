/*
 * pmh_sankoff.cpp
 *
 *  Created on: 12-jan-2017
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "sankofflabeling.h"
#include "nonbinaryclonetree.h"
#include <lemon/arg_parser.h>
#include <fstream>
#include <boost/algorithm/string.hpp>

void runSankoff(const NonBinaryCloneTree& T,
                const std::string& primary,
                const std::string& outputPrefix,
                const StringToIntMap& colorMap)
{
  SankoffLabeling sankoff(T, primary);
  
  sankoff.run();
  const int nrLabelings = sankoff.getNrLabelings();
  std::cerr << "Found " << nrLabelings << " maximum parsimony labelings with primary '" << primary << "'" << std::endl;
  
  SankoffLabeling::IntTripleToIntMap result = sankoff.classify();
  for (const auto& kv : result)
  {
    std::cerr << "Found " << kv.second << " labelings with " << kv.first.first << " comigrations, " << kv.first.second.first << " seeding sites and " << MigrationGraph::getPatternString(static_cast<MigrationGraph::Pattern>(kv.first.second.second)) << std::endl;
  }
  
  char buf[1024];
  for (int solIdx = 0; solIdx < nrLabelings; ++solIdx)
  {
    MigrationGraph G = sankoff.getMigrationGraph(solIdx);
    std::cerr << "Labeling " << solIdx << ": "
      << G.getNrMigrations() << " migrations, "
      << G.getNrComigrations(T, sankoff.getLabeling(solIdx)) << " comigrations, "
      << G.getNrSeedingSamples() << " seeding sites and "
      << G.getPatternString(G.getPattern(T, sankoff.getLabeling(solIdx)));
    std::cerr << std::endl;
    
    if (!outputPrefix.empty())
    {
      snprintf(buf, 1024, "%s/T-%s-%d.dot", outputPrefix.c_str(), primary.c_str(), solIdx);
      std::ofstream outT(buf);
      T.writeDOT(outT, sankoff.getLabeling(solIdx), colorMap);
      outT.close();
      
      snprintf(buf, 1024, "%s/G-%s-%d.dot", outputPrefix.c_str(), primary.c_str(), solIdx);
      std::ofstream outG(buf);
      G.writeDOT(outG, colorMap);
      outG.close();
    }
  }
}

int main(int argc, char** argv)
{
  std::string primarySample;
  std::string outputPrefix;
  std::string filenameColorMap;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("p", "Primary sample", primarySample)
    .refOption("o", "Output prefix", outputPrefix)
    .refOption("c", "Color map file", filenameColorMap)
    .other("T", "Clone tree")
    .other("leaf_labeling", "Leaf labeling");
  ap.parse();
  
  StringVector primaryVector;
  if (!primarySample.empty())
  {
    boost::split(primaryVector, primarySample, boost::is_any_of(","));
  }

  if (ap.files().size() != 2)
  {
    std::cerr << "Error: <T> and <leaf_labeling> must be specified" << std::endl;
    return 1;
  }
  
  std::string filenameT = ap.files()[0];
  std::ifstream inT(filenameT.c_str());
  if (!inT.good())
  {
    std::cerr << "Could not open '" << filenameT << "' for reading" << std::endl;
    return 1;
  }
  
  NonBinaryCloneTree T;
  if (!T.read(inT)) return 1;
  
  std::string filenameLabeling = ap.files()[1];
  std::ifstream inLabeling(filenameLabeling.c_str());
  if (!inLabeling.good())
  {
    std::cerr << "Could not open '" << filenameLabeling << "' for reading" << std::endl;
    return 1;
  }
  
  if (!T.readLeafLabeling(inLabeling)) return 1;
  
  StringToIntMap colorMap;
  if (!filenameColorMap.empty())
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
  else
  {
    colorMap = T.generateColorMap();
  }
  
  std::cerr << "Clone tree has " << T.getNrSamples() << " samples" << std::endl;

  if (primaryVector.empty())
  {
    for (const std::string& primary : T.getSamples())
    {
      runSankoff(T, primary, outputPrefix, colorMap);
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
        runSankoff(T, primary, outputPrefix, colorMap);
      }
    }
  }
  
  return 0;
}
