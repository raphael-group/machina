/*
 *  msenumeration.cpp
 *
 *   Created on: 17-aug-2017
 *       Author: M. El-Kebir
 */

#include "msenumeration.h"
#include "spruce/rootedcladisticnoisyancestrygraph.h"
#include "rootedcladisticnoisysparseenumeration.h"
#include "migrationgraph.h"

MSEnumeration::MSEnumeration(const FrequencyMatrix& F,
                             const std::string& primary,
                             const std::string& outputDirectory,
                             const StringToIntMap& colorMap)
  : _F(F)
  , _primary(primary)
  , _outputDirectory(outputDirectory)
  , _colorMap(colorMap)
  , _m(F.getNrAnatomicalSites())
  , _k(F.getNrSamples())
  , _n(F.getNrCharacters())
  , _Sigma(_F.getAnatomicalSites())
  , _F_lb()
  , _F_ub()
  , _S()
  , _barS()
  , _solutions()
{
  initRealTensors();
  initStateTrees();
}

void MSEnumeration::initRealTensors()
{
  _F_lb = RealTensor(_m, _k, _n + 1);
  _F_ub = RealTensor(_m, _k, _n + 1);
  
  for (int s = 0; s < _m; ++s)
  {
    for (int p = 0; p < _k; ++p)
    {
      double max_f_p_lb = 0;
      double max_f_p_ub = 0;
      for (int i = 0; i < _n; ++i)
      {
        if (_F.min(p, i) > max_f_p_lb)
        {
          max_f_p_lb = _F.min(p, i);
        }
        if (_F.max(p, i) > max_f_p_ub)
        {
          max_f_p_ub = _F.max(p, i);
        }
        
        if (s == 0)
        {
          _F_lb.set(s, p, i, 1 - _F.max(p, i));
          _F_ub.set(s, p, i, 1 - _F.min(p, i));
        }
        else if (s == 1)
        {
          _F_lb.set(s, p, i, _F.min(p, i));
          _F_ub.set(s, p, i, _F.max(p, i));
        }
        else
        {
          _F_lb.set(s, p, i, 0);
          _F_ub.set(s, p, i, 0);
        }
      }
      if (s == 0 && _F.anatomicalSiteIndexToSampleIndices(s).count(p) == 1)
      {
        _F_lb.set(s, p, _n, 1);
        _F_ub.set(s, p, _n, 1);
      }
      else if (_F.anatomicalSiteIndexToSampleIndices(s).count(p) == 1)
      {
        _F_lb.set(s, p, _n, max_f_p_lb);
        _F_ub.set(s, p, _n, max_f_p_ub);
      }
      else
      {
        _F_lb.set(s, p, _n, 0);
        _F_ub.set(s, p, _n, 0);
      }
    }
  }

  for (int i = 0; i < _n; ++i)
  {
    _F_lb.setColLabel(i, _F.indexToCharacter(i));
    _F_ub.setColLabel(i, _F.indexToCharacter(i));
  }
  _F_lb.setColLabel(_n, "SITE");
  _F_ub.setColLabel(_n, "SITE");
  
  for (int p = 0; p < _k; ++p)
  {
    _F_lb.setRowLabel(p, _F.indexToSample(p));
    _F_ub.setRowLabel(p, _F.indexToSample(p));
  }
}

void MSEnumeration::initStateTrees()
{
  IntVector pi(_m);
  pi[0] = -1;
  pi[1] = 0;
  for (int i = 2; i < _m; ++i)
  {
    pi[i] = -2;
  }
  StateTree S(pi);
  
  for (int i = 0; i < _n; ++i)
  {
    _S.push_back(S);
  }
}

int MSEnumeration::run(const MigrationTree& migrationTree,
                       bool force_mS,
                       const std::string& migrationTreeString,
                       int limit)
{
  char buf[1024];
  gm::SolutionSet solutions;
  
  _S.push_back(StateTree(migrationTree, _F.getIndexToAnatomicalSites()));
  
  gm::RootedCladisticNoisyAncestryGraph G(_S, _F_lb, _F_ub);
  G.initMultiState(_F.getSampleIndexToAnatomicalSiteIndex(),
                   _F.getAnatomicalSiteIndexToSampleIndices());
  G.setLabels(_F_lb);
  
  if (!_outputDirectory.empty())
  {
    std::ofstream outMigrationTree(_outputDirectory + "/" + "S" + migrationTreeString + ".dot");
    if (_colorMap.empty())
    {
      migrationTree.writeDOT(outMigrationTree);
    }
    else
    {
      migrationTree.writeDOT(outMigrationTree, _colorMap);
    }
    outMigrationTree.close();
    
    std::ofstream outAncestryGraph(_outputDirectory + "/" + "searchG" + migrationTreeString + ".dot");
    G.writeDOT(outAncestryGraph);
    outAncestryGraph.close();
  }

  IntSet whiteList;
  if (force_mS)
  {
    whiteList.insert(_n);
  }
  
  g_verbosity = VERBOSE_NON_ESSENTIAL;
  gm::RootedCladisticNoisySparseEnumeration enumerate(G, limit, -1, 1,
                                                      force_mS ? _n + 1 : _n,
                                                      true, false, whiteList);
  enumerate.run();
//  std::cout << "Objective value: " << enumerate.objectiveValue() << std::endl;
  
  if (enumerate.objectiveValue() == _n + 1)
  {
    enumerate.populateSolutionSet(solutions);
  }
  
  _S.pop_back();
  
  SolutionVector localSolutions;
  
  int nrSolutions = solutions.solutionCount();
  for (int i = 0; i < nrSolutions; ++i)
  {
    SolutionEntry* pSolution = constructCloneTree(solutions.solution(i));
    localSolutions.push_back(pSolution);
  }
  
  if (!_outputDirectory.empty())
  {
    int i = 0;
    for (const SolutionEntry* entry : localSolutions)
    {
      snprintf(buf, 1024, "%s/T%s%d.dot", _outputDirectory.c_str(), migrationTreeString.c_str(), i);
      std::ofstream outFile(buf);
      if (_colorMap.empty())
      {
        entry->writeDOT(outFile);
      }
      else
      {
        entry->writeDOT(outFile, _colorMap);
      }
      outFile.close();
      
      snprintf(buf, 1024, "%s/T%s%d.tree", _outputDirectory.c_str(), migrationTreeString.c_str(), i);
      std::ofstream outT(buf);
      entry->_T.write(outT);
      outT.close();
      
      snprintf(buf, 1024, "%s/T%s%d.labeling", _outputDirectory.c_str(), migrationTreeString.c_str(), i);
      std::ofstream outLabeling(buf);
      entry->_T.writeVertexLabeling(outT, entry->_lPlus);
      outLabeling.close();
      
//      snprintf(buf, 1024, "%s/phyloT%s%d.dot", _outputDirectory.c_str(), migrationTreeString.c_str(), i);
//      std::ofstream outFileT(buf);
//      gm::PerfectPhyloTree phyloT(solutions.solution(i).A(), solutions.solution(i).S());
//      phyloT.writeDOT(outFileT);
//      outFileT.close();

      snprintf(buf, 1024, "%s/G%s%d.dot", _outputDirectory.c_str(), migrationTreeString.c_str(), i);
      std::ofstream outFileG(buf);
      MigrationGraph G(entry->_T, entry->_lPlus);
      
      if (_colorMap.empty())
      {
        G.writeDOT(outFileG);
      }
      else
      {
        G.writeDOT(outFileG, _colorMap);
      }
      outFileG.close();
      
      ++i;
    }
  }
  
  _solutions.insert(_solutions.end(), localSolutions.begin(), localSolutions.end());
  
  return localSolutions.size();
}

void MSEnumeration::run(bool force_mS, int limit)
{
  char buf[1024];
  
  std::cerr << std::endl;
  std::cerr << "Initializing migration trees..." << std::endl;
  StringSet metastases = _F.getAnatomicalSites();
  assert(metastases.count(_primary) == 1);
  metastases.erase(_primary);
  MigrationTree::enumerate(_primary, metastases, _barS);
  std::cerr << "Found " << _barS.size() << " migration trees" << std::endl;
  
  int count = 0;
  for (const MigrationTree& S_c : _barS)
  {
    std::cerr << "Considering migration tree #" << count + 1 << std::endl;
    snprintf(buf, 1024, "%d_", count);
    int n = run(S_c, force_mS, buf, limit);
    ++count;
    std::cerr << "Found " << n << " clone trees with an mS pattern" << std::endl;
  }
}

MSEnumeration::SolutionEntry* MSEnumeration::constructCloneTree(const gm::Solution& solution)
{
  gm::PerfectPhyloTree T(solution.A(), solution.S());
  
  Digraph newTree;
  Node newRoot = lemon::INVALID;
  DoubleVectorNodeMap frequency(newTree);
  DoubleVectorNodeMap usage(newTree);
  IntNodeMap vertexToMutation(newTree);
  StringNodeMap nodeLabel(newTree);
  
  const Digraph& tree = T.T();
  NodeNodeMap oldTreeToNewTree(tree, lemon::INVALID);
  NodeVector mutationToVertex(_n, lemon::INVALID);
  std::vector<NodeSet> anatomicalSiteToVertices(_m, NodeSet());
  
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

      frequency[vv] = DoubleVector(_k, 0);
      usage[vv] = DoubleVector(_k, 0);
      
      const IntPair ci = T.nodeToCharState(v);
      const int col_idx = ci.first == 0 && ci.second == 0 ? 0 : (_n + 1) * (ci.second - 1) + ci.first + 1;
      
      vertexToMutation[vv] = ci.first;
      nodeLabel[vv] = ci.first == _n ? _F.indexToAnatomicalSite(ci.second) : _F.indexToCharacter(ci.first);
      if (ci.first < _n)
      {
        mutationToVertex[ci.first] = vv;
      }

      for (int p = 0; p < _k; ++p)
      {
        const int s = _F.sampleIndexToAnatomicalSiteIndex(p);
        
        frequency[vv][p] = solution.inferredF()(1, p, ci.first);
        usage[vv][p] = solution.U()(p, col_idx);
        if (g_tol.nonZero(usage[vv][p]))
        {
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
  for (int s = 0; s < _m; ++s)
  {
    for (Node v : anatomicalSiteToVertices[s])
    {
      if (OutArcIt(newTree, v) != lemon::INVALID)
      {
        Node leaf = newTree.addNode();
        newTree.addArc(v, leaf);
        usage[leaf] = DoubleVector();
        for (int p : _F.anatomicalSiteIndexToSampleIndices(s))
        {
          usage[leaf].push_back(usage[v][p]);
        }
        frequency[leaf].clear();
        label[leaf] = _F.indexToAnatomicalSite(s);
        
        if (vertexToMutation[v] == _n)
        {
          InArcIt a(newTree, v);
          while (a != lemon::INVALID)
          {
            int d = vertexToMutation[newTree.target(a)];
            if (d != _n)
            {
              nodeLabel[leaf] = _F.indexToCharacter(d) + "_" + _F.indexToAnatomicalSite(s);
            }
            Node u = newTree.source(a);
            a = InArcIt(newTree, u);
          }
        }
        else
        {
          nodeLabel[leaf] = _F.indexToCharacter(vertexToMutation[v]) + "_" + _F.indexToAnatomicalSite(s);
        }
      }
      else
      {
        label[v] = _F.indexToAnatomicalSite(s);
        DoubleVector tmp = usage[v];
        usage[v].clear();
        for (int p : _F.anatomicalSiteIndexToSampleIndices(s))
        {
          usage[v].push_back(tmp[p]);
        }
      }
    }
  }
  
  SolutionEntry* pSolution = new SolutionEntry(newTree, newRoot, nodeLabel, label);
  labelVerticesByAnatomicalSites(pSolution->_T,
                                 pSolution->_T.root(),
                                 _primary,
                                 pSolution->_lPlus);
  
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

void MSEnumeration::labelVerticesByAnatomicalSites(const CloneTree& T,
                                                   Node v,
                                                   const std::string& sStr,
                                                   StringNodeMap& lPlus) const
{
  if (_Sigma.count(T.label(v)) == 1)
  {
    lPlus[v] = T.label(v);
  }
  else
  {
    lPlus[v] = sStr;
  }
  
  for (OutArcIt vw(T.tree(), v); vw != lemon::INVALID; ++vw)
  {
    labelVerticesByAnatomicalSites(T, T.tree().target(vw), lPlus[v], lPlus);
  }
}
