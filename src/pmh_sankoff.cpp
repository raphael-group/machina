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

int main(int argc, char** argv)
{
  std::string primarySample;
  std::string outputDirectory;
  std::string filenameColorMap;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("p", "Primary sample (if omitted, every sample is considered iteratively as the primary)", primarySample)
    .refOption("o", "Output prefix", outputDirectory)
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
      SankoffLabeling::run(T, primary, outputDirectory, colorMap);
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
        SankoffLabeling::run(T, primary, outputDirectory, colorMap);
      }
    }
  }
  
  return 0;
}
