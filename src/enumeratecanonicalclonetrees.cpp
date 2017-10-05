/*
 * enumeratecanonicalclonetrees.cpp
 *
 *  Created on: 07-sep-2017
 *      Author: M. El-Kebir
 */

#include "enumeratecanonicalclonetrees.h"
#include "rootedcladisticnoisysparseenumeration.h"
#include "spruce/rootedcladisticnoisyenumeration.h"

EnumerateCanonicalCloneTrees::EnumerateCanonicalCloneTrees(const FrequencyMatrix& F)
  : _F(F)
{
}

void EnumerateCanonicalCloneTrees::enumerate(const std::string& outputDirectory,
                                             const int nrThreads,
                                             const int limit,
                                             TreeVector& canonicalCloneTrees) const
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
  
  gm::RootedCladisticNoisySparseEnumeration enumerate(G, limit, -1,
                                                      nrThreads,
                                                      1, //_F.getNrCharacters(),
                                                      true,
                                                      false,
                                                      IntSet());
  
  enumerate.run();
  
  gm::SolutionSet solutionSet;
  enumerate.populateSolutionSet(solutionSet);
  
  // 4. Transform enumerated mutation trees to the right format
  canonicalCloneTrees.clear();
  
  int nrSolutions = solutionSet.solutionCount();
  for (int idx = 0; idx < nrSolutions; ++idx)
  {
    const gm::Solution& sol = solutionSet.solution(idx);
    gm::PerfectPhyloTree phyloT(sol.A(), sol.S());
    
    CanonicalCloneTree* pTree = constructCanonicalCloneTree(sol);
    canonicalCloneTrees.push_back(pTree->_T);
    
    if (!outputDirectory.empty())
    {
      char buf[1024];
      snprintf(buf, 1024, "%s/T%d.dot", outputDirectory.c_str(), idx);
      
      std::ofstream outBarT(buf);
      pTree->_T.writeDOT(outBarT);
      outBarT.close();
      
      snprintf(buf, 1024, "%s/T%d.tree", outputDirectory.c_str(), idx);
      std::ofstream outT(buf);
      pTree->_T.write(outT);
      outT.close();
      
      snprintf(buf, 1024, "%s/T%d.labeling", outputDirectory.c_str(), idx);
      std::ofstream outHatL(buf);
      pTree->_T.writeLeafLabeling(outHatL);
      outHatL.close();
    }
    
    delete pTree;
  }
  
  std::cerr << "Found " << canonicalCloneTrees.size() << " canonical clone trees." << std::endl;
}

EnumerateCanonicalCloneTrees::CanonicalCloneTree* EnumerateCanonicalCloneTrees::constructCanonicalCloneTree(const gm::Solution& solution) const
{
  const int m = _F.getNrAnatomicalSites();
  const int n = _F.getNrCharacters();
  const int k = _F.getNrSamples();
  
  gm::PerfectPhyloTree T(solution.A(), solution.S());
  
  Digraph newTree;
  DoubleVectorNodeMap frequency(newTree);
  DoubleVectorNodeMap usage(newTree);
  
  const Digraph& tree = T.T();
  NodeNodeMap oldTreeToNewTree(tree, lemon::INVALID);
  Node newRoot = lemon::INVALID;
  IntNodeMap vertexToMutation(newTree);
  NodeVector mutationToVertex(n, lemon::INVALID);
  NodeVector anatomicalSiteToVertex(m, lemon::INVALID);
  std::vector<NodeSet> anatomicalSiteToVertices(m);
  StringNodeMap nodeLabel(newTree);
  
  // Copy the vertices:
  // ==================
  // (1) Assign a unique mutation from [n] to every vertex
  // (2) Assign frequencies
  // (3) Assign usages and identify leaves
  for (NodeIt v(tree); v != lemon::INVALID; ++v)
  {
    if (v != T.root())
    {
      Node vv = newTree.addNode();
      oldTreeToNewTree[v] = vv;
      IntPair ci = T.nodeToCharState(v);
      vertexToMutation[vv] = ci.first;
      nodeLabel[vv] = _F.indexToCharacter(ci.first);
      mutationToVertex[ci.first] = vv;
      frequency[vv] = DoubleVector(k, 0);
      usage[vv] = DoubleVector(k, 0);
      for (int p = 0; p < k; ++p)
      {
        frequency[vv][p] = solution.inferredF()(1, p, ci.first);
        
        const int col_idx = ci.first == 0 && ci.second == 0 ? 0 : (n + 1) * (ci.second - 1) + ci.first + 1;
        usage[vv][p] = solution.U()(p, col_idx);
        if (g_tol.nonZero(usage[vv][p]))
        {
          int s = _F.sampleIndexToAnatomicalSiteIndex(p);
          anatomicalSiteToVertices[s].insert(vv);
        }
      }
    }
  }
  
  for (ArcIt uv(tree); uv != lemon::INVALID; ++uv)
  {
    Node u = tree.source(uv);
    Node v = tree.target(uv);
    if (u != T.root())
    {
      Node uu = oldTreeToNewTree[u];
      Node vv = oldTreeToNewTree[v];
      
      newTree.addArc(uu, vv);
    }
    else
    {
      newRoot = oldTreeToNewTree[v];
    }
  }
  
  // Determine leaf labeling
  StringNodeMap label(newTree);
  for (int s = 0; s < m; ++s)
  {
    for (Node v : anatomicalSiteToVertices[s])
    {
      Node leaf = newTree.addNode();
      newTree.addArc(v, leaf);
      usage[leaf] = DoubleVector(0);
      for (int p : _F.anatomicalSiteIndexToSampleIndices(s))
      {
        usage[leaf].push_back(usage[v][p]);
      }
      frequency[leaf].clear();
      label[leaf] = _F.indexToAnatomicalSite(s);
      nodeLabel[leaf] = _F.indexToCharacter(vertexToMutation[v]) + "_" + _F.indexToAnatomicalSite(s);
    }
  }
  
  CanonicalCloneTree* pSolution = new CanonicalCloneTree(newTree, newRoot, nodeLabel, label);
  
  for (NodeIt vv(newTree); vv != lemon::INVALID; ++vv)
  {
    const std::string label_vv = nodeLabel[vv];
    Node new_v = pSolution->_T.getNodeByLabel(label_vv);
    pSolution->_frequency[new_v] = frequency[vv];
    if (_F.characterToIndex(label_vv) != -1)
    {
      pSolution->_characterNodeMap[new_v] = label_vv;
    }
    if (OutArcIt(newTree, vv) == lemon::INVALID)
    {
      pSolution->_usage[new_v] = usage[vv];
    }
  }
  
  return pSolution;
}

void EnumerateCanonicalCloneTrees::getFrequencyTensor(RealTensor& F_lb,
                                         RealTensor& F_ub) const
{
  const int n = _F.getNrCharacters();
  const int k = _F.getNrSamples();
  
  F_lb = RealTensor(2, k, n);
  F_ub = RealTensor(2, k, n);
  
  for (int p = 0; p < k; ++p)
  {
    for (int i = 0; i < n; ++i)
    {
      F_lb.set(1, p, i, _F.min(p, i));
      F_ub.set(1, p, i, _F.max(p, i));
      F_lb.set(0, p, i, 1 - _F.max(p, i));
      F_ub.set(0, p, i, 1 - _F.min(p, i));
    }
  }
}
