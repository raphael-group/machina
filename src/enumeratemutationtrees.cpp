/*
 * enumeratemutationtrees.cpp
 *
 *  Created on: 07-sep-2017
 *      Author: M. El-Kebir
 */

#include "enumeratemutationtrees.h"
#include "rootedcladisticnoisysparseenumeration.h"
#include "spruce/rootedcladisticnoisyenumeration.h"

EnumerateMutationTrees::EnumerateMutationTrees(const FrequencyMatrix& F)
  : _F(F)
  , _nrCharactersInTrees(0)
{
}

void EnumerateMutationTrees::enumerate(const std::string& outputDirectory,
                                       const int nrThreads,
                                       const int limit,
                                       const int timeLimit,
                                       TreeVector& mutationTrees)
{
  // 1. Construct frequency tensors F_lb and F_ub
  RealTensor F_lb, F_ub;
  getFrequencyTensor(F_lb, F_ub);
  
  // 2. Construct state tree vector, each state tree is 0 -> 1
  IntVector pi(2);
  pi[0] = -1;
  pi[1] = 0;
  StateTree S(pi);
  StateTreeVector stateTrees(_F.getNrCharacters(), S);
  
  // 3. Enumerate all mutation trees
  gm::RootedCladisticNoisyAncestryGraph G(stateTrees, F_lb, F_ub);
  G.init();
  G.setLabels(F_lb);
  
  gm::RootedCladisticNoisyEnumeration enumerate(G, limit, timeLimit,
                                                nrThreads,
                                                1, //_F.getNrCharacters(),
                                                true,
                                                false,
                                                IntSet());
  
  enumerate.run();
  
  _nrCharactersInTrees = enumerate.objectiveValue();
  
  gm::SolutionSet solutionSet;
  enumerate.populateSolutionSet(solutionSet);
  
  // 4. Transform enumerated mutation trees to the right format
  mutationTrees.clear();
  
  int nrSolutions = solutionSet.solutionCount();
  for (int idx = 0; idx < nrSolutions; ++idx)
  {
    const gm::Solution& sol = solutionSet.solution(idx);
    gm::PerfectPhyloTree phyloT(sol.A(), sol.S());
    
    const Digraph& T = phyloT.T();
    NodeNodeMap old2new(T, lemon::INVALID);
    
    Digraph newT;
    Node newRoot = lemon::INVALID;
    StringNodeMap idMap(newT, "");
    StringNodeMap l(newT, "");
    
    // Copy vertices
    for (NodeIt v(T); v != lemon::INVALID; ++v)
    {
      if (v != phyloT.root())
      {
        Node vv = newT.addNode();
        old2new[v] =  vv;
        assert(phyloT.nodeToCharState(v).second == 1);
        idMap[vv] = _F.indexToCharacter(phyloT.nodeToCharState(v).first);
      }
    }
    
    for (ArcIt a(T); a != lemon::INVALID; ++a)
    {
      Node u = T.source(a);
      Node v = T.target(a);
      
      if (u != phyloT.root())
      {
        Node uu = old2new[u];
        Node vv = old2new[v];
        
        newT.addArc(uu, vv);
      }
      else
      {
        // phyloT needs to have monoclonal origin
        // (i.e. v_(*,0) has only a single child)
        assert(newRoot == lemon::INVALID);
        
        newRoot = old2new[v];
      }
    }
    
    mutationTrees.push_back(CloneTree(newT, newRoot, idMap, l));
    if (!outputDirectory.empty())
    {
      char buf[1024];
      snprintf(buf, 1024, "%s/barT%d.tree", outputDirectory.c_str(), idx);
      
      std::ofstream outBarT(buf);
      
      mutationTrees.back().write(outBarT);
      
      outBarT.close();
    }
  }
  
  std::cerr << "Found " << mutationTrees.size() << " mutation trees with "
            << _nrCharactersInTrees << " out of "
            << _F.getNrCharacters() << " mutations" << std::endl;
}

void EnumerateMutationTrees::getFrequencyTensor(RealTensor& F_lb,
                                                RealTensor& F_ub) const
{
  const int n = _F.getNrCharacters();
  const int k = _F.getNrSamples();
  
  F_lb = RealTensor(2, k, n);
  F_ub = RealTensor(2, k, n);
  
  for (int p = 0; p < k; ++p)
  {
    F_lb.setRowLabel(p, _F.indexToSample(p));
    F_ub.setRowLabel(p, _F.indexToSample(p));
    for (int i = 0; i < n; ++i)
    {
      if (p == 0)
      {
        F_lb.setColLabel(i, _F.indexToCharacter(i));
        F_ub.setColLabel(i, _F.indexToCharacter(i));
      }
      F_lb.set(1, p, i, _F.min(p, i));
      F_ub.set(1, p, i, _F.max(p, i));
      F_lb.set(0, p, i, 1 - _F.max(p, i));
      F_ub.set(0, p, i, 1 - _F.min(p, i));
    }
  }
}

std::ostream& operator<<(std::ostream& out,
                         const EnumerateMutationTrees::TreeVector& trees)
{
  out << trees.size() << " #trees" << std::endl;
  int idx = 0;
  for (const CloneTree& T : trees)
  {
    out << lemon::countArcs(T.tree()) << " #edges, tree " << idx + 1 << std::endl;
    T.write(out);
    ++idx;
  }
  return out;
}

std::istream& operator>>(std::istream& in,
                         EnumerateMutationTrees::TreeVector& trees)
{
  int nrTrees = -1;
  
  std::string line;
  getline(in, line); // skip first line
  getline(in, line);
  
  std::stringstream ss(line);
  ss >> nrTrees;
  
  if (nrTrees < 0)
  {
    throw std::runtime_error("Error: number of trees should be nonnegative");
  }
  
  for (int i = 0; i < nrTrees; ++i)
  {
    int nrEdges = -1;
    
    getline(in, line);
    ss.clear();
    ss.str(line);
    ss >> nrEdges;

    if (nrEdges < 0)
    {
      throw std::runtime_error("Error: number of edges should be nonnegative");
    }
    
    ss.clear();
    ss.str("");
    for (int j = 0; j < nrEdges; ++j)
    {
      getline(in, line);
      ss << line << std::endl;
    }
    
    CloneTree T;
    if (!T.read(ss))
      throw std::runtime_error("");
    
    trees.push_back(T);
  }
  
  return in;
}
