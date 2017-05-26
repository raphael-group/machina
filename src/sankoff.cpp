/*
 * sankoff.cpp
 *
 *  Created on: 10-jan-2017
 *      Author: M. El-Kebir
 */

#include "sankoff.h"

Sankoff::Sankoff(const CharacterTree& T)
  : MPBase(T)
  , _M(T.getNrCharacters(), NULL)
{
}

void Sankoff::run(int rootState)
{
  const int nrCharacters = _T.getNrCharacters();
  
  // clear and initialize M
  for (int c = 0; c < nrCharacters; ++c)
  {
    delete _M[c];
    _M[c] = new EntryMatrix(_T.tree(), EntryVector(_T.getNrStates(c)));
  }
  
  Node root =  _T.root();
  for (int c = 0; c < nrCharacters; ++c)
  {
    const int nrStates = _T.getNrStates(c);
    
    run(c, root);

    // backtrace
    NodeStatePairList frontier;
    if (rootState == -1)
    {
      // 1. find minimum cost state for the root vertex
      int min_cost = std::numeric_limits<int>::max();
      int min_s = -1;
      for (int s = 0; s < nrStates; ++s)
      {
        const Entry& entry_s = entry(c, root, s);
        if (entry_s._cost < min_cost)
        {
          min_cost = entry_s._cost;
          min_s = s;
        }
      }
      
      assert(min_s != -1);
      for (int s = 0; s < nrStates; ++s)
      {
        if (entry(c, root, s)._cost == min_cost)
        {
          frontier.push_back(std::make_pair(root, s));
        }
      }
    }
    else
    {
      frontier.push_back(std::make_pair(root, rootState));
    }
    
    IntNodeMap stateVector(_T.tree(), -1);
    for (Node u : _T.leafSet())
    {
      stateVector[u] = _T.state(u, c);
    }
    constructBackTrace(c, stateVector, frontier);
  }
  
  updateHomoplasy();
}

void Sankoff::constructBackTrace(int c,
                                 IntNodeMap& stateVector,
                                 NodeStatePairList& frontier)
{
  const Digraph& T = _T.tree();
  
  if (frontier.empty())
  {
    // report solution
    IntNodeMap* pStateVector = new IntNodeMap(T, -1);
    _stateVector[c].push_back(pStateVector);
    _homoplasyVector[c].push_back(true);
    
    for (NodeIt w(_T.tree()); w != lemon::INVALID; ++w)
    {
      assert(stateVector[w] >= 0);
      pStateVector->set(w, stateVector[w]);
    }
  }
  else
  {
    do
    {
      NodeStatePair vs = frontier.front();
      stateVector[vs.first] = vs.second;
      frontier.pop_front();
      const Entry& entry_vs = entry(c, vs.first, vs.second);
      
      // remove elements from frontier with the same vertex v
      NodeStatePairList newFrontier = frontier;
      bool removed = false;
      for (auto it = newFrontier.begin(); it != newFrontier.end();)
      {
        if (it->first == vs.first)
        {
          it = newFrontier.erase(it);
          removed = true;
        }
        else
        {
          ++it;
        }
      }
      
      // add children and their options to frontier
      for (const auto& wS : entry_vs._previousEntries)
      {
        if (_T.isLeaf(wS.first)) continue;
        for (int t : wS.second)
        {
          newFrontier.push_front(std::make_pair(wS.first, t));
        }
      }
      constructBackTrace(c, stateVector, newFrontier);
      
      if (!removed)
        break;
    } while (!frontier.empty());
  }
}

void Sankoff::constructBackTrace(int c, Node u, int s)
{
  const Digraph& T = _T.tree();
  
  (*_stateVector[c][0])[u] = s;
  const Entry& entry_cus = entry(c, u, s);
  
  if (!_T.isLeaf(u))
  {
    assert(!entry_cus._previousEntries.empty());
    
    for (Digraph::OutArcIt a(T, u); a != lemon::INVALID; ++a)
    {
      Node v = T.target(a);
      assert(entry_cus._previousEntries.count(v) > 0);
      
      int t = *entry_cus._previousEntries.find(v)->second.begin();
      constructBackTrace(c, v, t);
    }
  }
}

void Sankoff::run(int c, Node u)
{
  const Digraph& T = _T.tree();
  const int nrStates = _T.getNrStates(c);
  
  if (_T.isLeaf(u))
  {
    for (int s = 0; s < nrStates; ++s)
    {
      Entry& entry_cus = entry(c, u, s);
      if (s == _T.state(u, c))
      {
        entry_cus._cost = 0;
      }
      else
      {
        entry_cus._cost = std::numeric_limits<int>::max();
      }
    }
  }
  else
  {
    // solve children
    for (Digraph::OutArcIt a(T, u); a != lemon::INVALID; ++a)
    {
      Node v = T.target(a);
      run(c, v);
    }
    
    // compute _M[characterIndex][u][s]
    for (int s = 0; s < nrStates; ++s)
    {
      Entry& entry_cus = entry(c, u, s);
      entry_cus._cost = 0;
      
      for (Digraph::OutArcIt a(T, u); a != lemon::INVALID; ++a)
      {
        // compute
        // \sum_{v \in \delta_T(u)} \min_{t \in \Sigma} { c_{s,t} + M[v,t] }
        Node v = T.target(a);
        
        int minimum_cost = std::numeric_limits<int>::max();
        for (int t = 0; t < nrStates; ++t)
        {
          int cost_vt = computeCost(c, u, s, v, t);
          if (cost_vt < minimum_cost)
          {
            minimum_cost = cost_vt;
          }
        }
        
        // now update _M[characterIndex][u][s]
        entry_cus._cost += minimum_cost;
        
        for (int t = 0; t < nrStates; ++t)
        {
          int cost_vt = computeCost(c, u, s, v, t);
          if (cost_vt == minimum_cost)
          {
            entry_cus._previousEntries[v].insert(t);
          }
        }
      }
    }
  }
}

void Sankoff::write(int c, std::ostream& out) const
{
  assert(0 <= c && c < _T.getNrCharacters());
  
  const Digraph& T = _T.tree();
  const int nrStates = _T.getNrStates(c);
  
  for (NodeIt u(T); u != lemon::INVALID; ++u)
  {
    out << _T.label(u);
    for (int s = 0; s < nrStates; ++s)
    {
      out << "\t" << entry(c, u, s)._cost;
    }
    out << std::endl;
  }
  
  for (NodeIt u(T); u != lemon::INVALID; ++u)
  {
    out << _T.label(u);
    for (int s = 0; s < nrStates; ++s)
    {
      out << "\t";
      bool first = true;
      for (const auto& vt : entry(c, u, s)._previousEntries)
      {
        if (first)
          first = false;
        else
          out << ",";
        
        out << "(" << _T.label(vt.first);
        
        const IntSet& states = vt.second;
        for (int t : states)
        {
          out << "," << t;
        }
        out << ")";
      }
    }
    out << std::endl;
  }
}
