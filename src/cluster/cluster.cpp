/*
 * cluster.cpp
 *
 *  Created on: 7-sep-2017
 *      Author: M. El-Kebir
 */

#include "cluster.h"
#include <lemon/connectivity.h>

Cluster::Cluster(const ReadMatrix& R,
                 double alpha,
                 int varThreshold,
                 bool relabel)
  : _R(R)
  , _F(R.toFrequencyMatrix(alpha, varThreshold))
  , _profile(R.getNrSamples(),
             ProfileVector(R.getNrCharacters(),
                           ABSENT))
  , _clustering()
  , _relabel(relabel)
{
  classify();
}

void Cluster::writeAncesTreeInput(std::ostream& out) const
{
  const int n = _newR.getNrCharacters();
  const int k = _newR.getNrSamples();
  
  out << "cluster";
  for (int p = 0; p < k; ++p)
  {
    out << "\t" << _newR.indexToSample(p)
        << "\t" << _newR.indexToSample(p);
  }
  out << std::endl;
  
  for (int i = 0; i < n; ++i)
  {
    out << _newR.indexToCharacter(i);
    for (int p = 0; p < k; ++p)
    {
      out << "\t" << _newR.getRef(p, i)
          << "\t" << _newR.getVar(p, i);
    }
    out << std::endl;
  }
}

void Cluster::classify()
{
  const int k = _F.getNrSamples();
  const int n = _F.getNrCharacters();

  for (int p = 0; p < k; ++p)
  {
    for (int i = 0; i < n; ++i)
    {
      double f_lb_pi = _F.min(p, i);
      double f_ub_pi = _F.max(p, i);
      
      if (f_lb_pi == 0)
      {
        _profile[p][i] = ABSENT;
      }
      else
      {
        _profile[p][i] = CLONAL;
        for (int j = 0; j < n; ++j)
        {
          if (i == j) continue;
          if (f_ub_pi < _F.min(p, j))
          {
            _profile[p][i] = SUBCLONAL;
            break;
          }
        }
      }
    }
  }
}

void Cluster::clusterClonalityStatus(double beta)
{
  const int n = _F.getNrCharacters();
  
  std::map<ProfileVector, IntVector> C;
  for (int i = 0; i < n; ++i)
  {
    ProfileVector profile_i = getProfile(i);
    C[profile_i].push_back(i);
  }
  
  _clustering.clear();
  for (const auto& kv : C)
  {
    bool allAbsent = true;
    for (MutationType type : kv.first)
    {
      if (type != ABSENT)
      {
        allAbsent = false;
        break;
      }
    }
    
    if (!allAbsent)
    {
      _clustering.push_back(kv.second);
    }
  }
  
  _newR = _R.poolReads(_clustering, _relabel);
  _newF = _newR.toFrequencyMatrix(beta, 3);
}

void Cluster::clusterCC(double beta)
{
  Graph G;
  Graph::NodeMap<int> nodeToMutation(G);

  std::vector<Graph::Node> mutationToNode(_F.getNrCharacters(), lemon::INVALID);
  
  const int n = _F.getNrCharacters();
  const int k = _F.getNrSamples();
  
  // 1. Add nodes, only if present in at least one sample
  for (int i = 0; i < n; ++i)
  {
    bool present = false;
    for (int p = 0; p < k; ++p)
    {
      if (_F.max(p, i) > 0)
      {
        present = true;
        break;
      }
    }
    
    if (present)
    {
      Graph::Node v_i = G.addNode();
      nodeToMutation[v_i] = i;
      mutationToNode[i] = v_i;
    }
  }
  
  // 2. Add edges
  for (int i = 0; i < n; ++i)
  {
    Graph::Node v_i = mutationToNode[i];
    if (v_i == lemon::INVALID) continue;

    for (int j = i + 1; j < n; ++j)
    {
      Graph::Node v_j = mutationToNode[j];
      if (v_j == lemon::INVALID) continue;

      bool ok = true;
      for (int p = 0; p < k; ++p)
      {
        double i_lb = _F.min(p, i);
        double i_ub = _F.max(p, i);

        double j_lb = _F.min(p, j);
        double j_ub = _F.max(p, j);
        
        double lb = std::max(i_lb, j_lb);
        double ub = std::min(i_ub, j_ub);
        
        if (lb > ub)
        {
          ok = false;
          break;
        }
      }
      
      if (ok)
      {
        G.addEdge(v_i, v_j);
      }
    }
  }
  
//  // 3. DEBUG write graph
//  std::cout << "graph G {" << std::endl;
//  for (Graph::NodeIt v_i(G); v_i != lemon::INVALID; ++v_i)
//  {
//    int i = nodeToMutation[v_i];
//    std::cout << "\t" << i
//              << " [label=\""
//              << _F.indexToCharacter(i)
//              << "\"]" << std::endl;
//  }
//  for (Graph::EdgeIt a_ij(G); a_ij != lemon::INVALID; ++a_ij)
//  {
//    Graph::Node v_i = G.u(a_ij);
//    Graph::Node v_j = G.v(a_ij);
//
//    int i = nodeToMutation[v_i];
//    int j = nodeToMutation[v_j];
//
//    std::cout << "\t" << i << " -- " << j << std::endl;
//  }
//  std::cout << "}" << std::endl;
  
  // 4. Find connected components
  Graph::NodeMap<int> compMap(G, -1);
  int nComps = lemon::connectedComponents(G, compMap);
  _clustering = IntMatrix(nComps);
  for (int i = 0; i < n; ++i)
  {
    Graph::Node v_i = mutationToNode[i];
    if (v_i == lemon::INVALID) continue;
    _clustering[compMap[v_i]].push_back(i);
  }
  
  _newR = _R.poolReads(_clustering, _relabel);
  _newF = _newR.toFrequencyMatrix(beta, 3);
}

void Cluster::readClustering(std::istream& in,
                             double beta)
{
  _clustering.clear();
  
  while (in.good())
  {
    std::string line;
    getline(in, line);

    if (line.empty()) continue;
    
    _clustering.push_back(IntVector());
    
    StringVector s;
    boost::split(s, line, boost::is_any_of(";"));
    
    IntVector ss;
    for (const std::string& cStr : s)
    {
      int c = _R.characterToIndex(cStr);
      _clustering.back().push_back(c);
    }
  }
  
  _newR = _R.poolReads(_clustering, _relabel);
  _newF = _newR.toFrequencyMatrix(beta, 3);
}

void Cluster::writeClustering(std::ostream& out) const
{
  int idx = 0;
  for (const IntVector& C : _clustering)
  {
    if (_relabel)
    {
      out << "cluster_" << idx + 1 << "\t";
    }
    bool first = true;
    for (int i : C)
    {
      if (first)
        first = false;
      else
        out << ";";
      out << _F.indexToCharacter(i);
    }
    out << std::endl;
    ++idx;
  }
//  for (const auat& kv : _clustering)
//  {
//    out << "{";
//    for (MutationType type : kv.first)
//    {
//      switch (type)
//      {
//        case ABSENT:
//          out << " A";
//          break;
//        case SUBCLONAL:
//          out << " S";
//          break;
//        case CLONAL:
//          out << " C";
//          break;
//      }
//    }
//    out << " } : ";
//    for (int i : kv.second)
//    {
//      out << " " << i;
//    }
//    out << std::endl;
//  }
}
