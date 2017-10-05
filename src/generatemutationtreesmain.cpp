/*
 *  generatemutationtrees.cpp
 *
 *   Created on: 24-aug-2017
 *       Author: M. El-Kebir
 */

#include "utils.h"
#include "frequencymatrix.h"
#include "rootedcladisticnoisysparseenumeration.h"
#include "enumeratecanonicalclonetrees.h"
#include "enumeratemutationtrees.h"
#include "clonetree.h"
#include <lemon/arg_parser.h>
#include <fstream>

int main(int argc, char** argv)
{
  std::string outputDirectory;
  bool canonicalCloneTrees = false;
  int nrThreads = 1;
  int limit = -1;
  int timeLimit = -1;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("o", "Output directory" , outputDirectory);
  ap.refOption("l", "Maximum number of mutation trees to enumerate (default: -1, unlimited)" , limit);
  ap.refOption("tl", "Time limit in seconds (default: -1, unlimited)" , timeLimit);
  ap.refOption("C", "Output canonical clone trees", canonicalCloneTrees);
  ap.refOption("t", "Number of threads (default: 1)", nrThreads);
  ap.other("frequencies", "Frequencies");
  ap.parse();
  
  if (ap.files().size() != 1)
  {
    std::cerr << "Error: <frequencies> must be specified" << std::endl;
    return 1;
  }
  
  std::string filenameFrequencies = ap.files()[0];
  
  FrequencyMatrix F;
  try
  {
    if (filenameFrequencies != "-")
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
    else
    {
      std::cin >> F;
    }
  }
  catch (std::runtime_error& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  
  g_verbosity = VERBOSE_NON_ESSENTIAL;
  if (canonicalCloneTrees)
  {
    EnumerateCanonicalCloneTrees::TreeVector trees;
    EnumerateCanonicalCloneTrees enumerate(F);
    enumerate.enumerate(outputDirectory, nrThreads, limit, trees);
  }
  else
  {
    EnumerateMutationTrees::TreeVector trees;
    EnumerateMutationTrees enumerate(F);
    enumerate.enumerate(outputDirectory,
                        nrThreads,
                        limit,
                        timeLimit,
                        trees);
    
    if (outputDirectory.empty())
    {
      std::cout << "# " << enumerate.getNrCharactersInTrees()
                << " out of " << F.getNrCharacters() << " mutations"
                << std::endl;
      std::cout << trees;
    }
  }
  
  return 0;
}
