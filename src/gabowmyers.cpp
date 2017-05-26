/*
 * gabowmyers.cpp
 *
 *  Created on: 18-oct-2016
 *      Author: M. El-Kebir
 */

#include "gabowmyers.h"
#include <lemon/bfs.h>

GabowMyers::GabowMyers(const Digraph& G, Node root)
  : _G(G)
  , _root(root)
  , _result()
{
}

void GabowMyers::run()
{
  BoolNodeMap filterNodesT(_G, false);
  BoolArcMap filterArcsT(_G, false);
  SubDigraph T(_G, filterNodesT, filterArcsT);
  
  BoolNodeMap filterNodesG(_G, true);
  BoolArcMap filterArcsG(_G, true);
  SubDigraph subG(_G, filterNodesG, filterArcsG);
  
  ArcList F;
  
  init(subG, T, F);
  grow(subG, T, F);
  
  std::cerr << std::endl;
}

void GabowMyers::init(SubDigraph& subG,
                      SubDigraph& T,
                      ArcList& F)
{
  T.enable(_root);
  F.clear();
  for (OutArcIt a(_G, _root); a != lemon::INVALID; ++a)
  {
    F.push_back(a);
  }
}

void GabowMyers::grow(SubDigraph& G,
                      SubDigraph& T,
                      ArcList& F)
{
  if (F.empty())
  {
    finalize(T);
  }
  else
  {
    ArcList FF;
    
    do
    {
      assert(!F.empty());
      
      Arc a_cidj = F.back();
      F.pop_back();
      
      Node v_ci = G.source(a_cidj);
      Node v_dj = G.target(a_cidj);
      
      assert(T.status(v_ci));
      assert(!T.status(v_dj));
      assert(!T.status(a_cidj));
      
      // add a_cidj to T
      addArc(T, a_cidj);
      
      ArcList newF = F;
      
      // remove each arc wv where w in T from F
      for (ArcListNonConstIt it = newF.begin(); it != newF.end();)
      {
        if (G.target(*it) == v_dj)
        {
          assert(T.status(G.source(*it)));
          it = newF.erase(it);
        }
//        else if (G.source(*it) == v_ci)
//        {
//          it = newF.erase(it);
//        }
        else
        {
          ++it;
        }
      }
      
      // push each arc a_djel where v_el not in V(T) onto F
      for (SubOutArcIt a_djel(G, v_dj); a_djel != lemon::INVALID; ++a_djel)
      {
        Node v_el = G.target(a_djel);
        
        // violation of tree constraint (no cycles)
        if (T.status(v_el))
          continue;
        
        newF.push_back(a_djel);
      }
      
      grow(G, T, newF);
      
      G.disable(a_cidj);
      
      removeArc(T, a_cidj);
      
      FF.push_back(a_cidj);
    } while (!F.empty());
    
    for (ArcListRevIt it = FF.rbegin(); it != FF.rend(); ++it)
    {
      Arc a = *it;
      assert(!G.status(a));
      
      F.push_back(*it);
      G.enable(a);
    }
  }
}

void GabowMyers::addArc(SubDigraph& T,
                        Arc a_cidj) const
{
  const Node v_dj = T.target(a_cidj);
  
  // add a_cidj to T
  T.enable(v_dj);
  T.enable(a_cidj);
  
  assert(isArborescence(T));
}

void GabowMyers::removeArc(SubDigraph& T,
                           Arc a_cidj) const
{
  assert(T.status(a_cidj));
  
  const Node v_dj = T.target(a_cidj);
  
  // remove a_cidj from T
  T.disable(a_cidj);
  T.disable(v_dj);
  
  assert(isArborescence(T));
}

bool GabowMyers::isArborescence(const SubDigraph& T) const
{
  assert(T.status(_root));
  lemon::Bfs<SubDigraph> bfs(T);
  bfs.run(_root);
  
  for (SubNodeIt v_ci(T); v_ci != lemon::INVALID; ++v_ci)
  {
    if (!bfs.reached(v_ci))
    {
      return false;
    }
  }
  
  return true;
}

void GabowMyers::finalize(const SubDigraph& T)
{
  assert(isArborescence(T));
//  writeDOT(T, std::cout);
  
  _result.push_back(ArcList());
  ArcList& res = _result.back();
  
  for (SubArcIt a_cidj(T); a_cidj != lemon::INVALID; ++a_cidj)
  {
    res.push_back(a_cidj);
  }
  std::cerr << "\rNumber of trees of size " << res.size() << ": " << _result.size() << std::flush;
}

void GabowMyers::writeDOT(const StringNodeMap& label,
                          std::ostream& out) const
{
  out << "digraph G {" << std::endl;
  
  for (NodeIt v(_G); v != lemon::INVALID; ++v)
  {
    out << "\t" << _G.id(v) << " [label=\"" << label[v] << "\"]" << std::endl;
  }
  
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    Node u = _G.source(a);
    Node v = _G.target(a);
    
    out << "\t" << _G.id(u) << " -> " << _G.id(v) << std::endl;
  }
  
  out << "}" << std::endl;
}

void GabowMyers::writeDOT(const SubDigraph& T, std::ostream& out) const
{
  out << "digraph T {" << std::endl;
  
  for (SubNodeIt v(T); v != lemon::INVALID; ++v)
  {
    out << "\t" << _G.id(v) << std::endl;
  }
  
  for (SubArcIt a(T); a != lemon::INVALID; ++a)
  {
    Node u = _G.source(a);
    Node v = _G.target(a);
    
    out << "\t" << _G.id(u) << " -> " << _G.id(v) << std::endl;
  }
  
  out << "}" << std::endl;
}

void GabowMyers::result(SubDigraph& T,
                        int i) const
{
  assert(0 <= i && i < _result.size());
  
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    T.disable(a);
  }
  
  for (NodeIt v(_G); v != lemon::INVALID; ++v)
  {
    T.enable(v);
  }
  
  ArcListList::const_iterator it = _result.begin();
  for (int j = 0; j < i; ++j)
  {
    ++it;
  }
  
  for (Arc a : *it)
  {
    T.enable(a);
  }
}
