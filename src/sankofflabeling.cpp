/*
 * sankofflabeling.cpp
 *
 *  Created on: 12-jan-2017
 *      Author: M. El-Kebir
 */

#include "sankofflabeling.h"
#include <fstream>

void SankoffLabeling::run()
{
  assert(_stateToSample[0] == _primary);
  Sankoff sankoff(_charT);
  sankoff.run(0);
  
  do
  {
    StringNodeMap* pLabeling = new StringNodeMap(_T.tree());
    for (NodeIt v(_charT.tree()); v != lemon::INVALID; ++v)
    {
      const Node vv = _T.getNodeByLabel(_charT.label(v));
      assert(vv != lemon::INVALID);
      
      const int s = sankoff.state(v, 0);
      assert(0 <= s && s < _stateToSample.size());
      
      const std::string& sample = _stateToSample[s];
      pLabeling->set(vv, sample);
    }
    _mpLabelings.push_back(pLabeling);
  } while (sankoff.nextSolution());
}

void SankoffLabeling::run(const NonBinaryCloneTree& T,
                          const std::string& primary,
                          const std::string& outputDirectory,
                          const StringToIntMap& colorMap)
{
  SankoffLabeling sankoff(T, primary);
  
  sankoff.run();
  const int nrLabelings = sankoff.getNrLabelings();
  std::cerr << "Found " << nrLabelings << " maximum parsimony labelings with primary '" << primary << "'" << std::endl;
  
  SankoffLabeling::IntTripleToIntMap result = sankoff.classify();
  for (const auto& kv : result)
  {
    MigrationGraph::Pattern pattern = static_cast<MigrationGraph::Pattern>(kv.first.second.second);
    std::cerr << "Found " << kv.second
      << " labelings with " << kv.first.first << " comigrations, "
      << kv.first.second.first << " seeding sites and "
      << MigrationGraph::getPatternString(pattern)
      << std::endl;
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
    
    if (!outputDirectory.empty())
    {
      snprintf(buf, 1024, "%s/T-%s-%d.dot", outputDirectory.c_str(), primary.c_str(), solIdx);
      std::ofstream outT(buf);
      T.writeDOT(outT, sankoff.getLabeling(solIdx), colorMap);
      outT.close();
      
      snprintf(buf, 1024, "%s/G-%s-%d.dot", outputDirectory.c_str(), primary.c_str(), solIdx);
      std::ofstream outG(buf);
      G.writeDOT(outG, colorMap);
      outG.close();
    }
  }
}
