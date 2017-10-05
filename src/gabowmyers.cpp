/*
 * gabowmyers.cpp
 *
 *  Created on: 18-oct-2016
 *      Author: M. El-Kebir
 */

#include "gabowmyers.h"
#include <lemon/bfs.h>

GabowMyers::GabowMyers(const Digraph& G, Node root, int limit)
  : _G(G)
  , _limit(limit)
  , _root(root)
  , _result()
  , _nrVertices(lemon::countNodes(_G))
{
}

void GabowMyers::run()
{
  _result.clear();
  
  BoolNodeMap nodesG(_G, true);
  BoolNodeMap nodesT(_G, false);
  BoolNodeMap nodesL(_G, false);
  
  BoolArcMap arcsG(_G, true);
  BoolArcMap arcsT(_G, false);
  BoolArcMap arcsL(_G, false);
  
  SubDigraph G(_G, nodesG, arcsG);
  SubDigraph T(_G, nodesT, arcsT);
  SubDigraph L(_G, nodesL, arcsL);
  
  SubBfs bfsL(L);
  bfsL.init();
  
  ArcList F;
  for (SubOutArcIt a(G, _root); a != lemon::INVALID; ++a)
  {
    F.push_back(a);
  }
  
  T.enable(_root);
  grow(G, T, L, bfsL, F);
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

bool GabowMyers::grow(SubDigraph& G,
                      SubDigraph& T,
                      SubDigraph& L,
                      SubBfs& bfsL,
                      ArcList& F)
{
  if (lemon::countNodes(T) == _nrVertices)
  {
    // clear L
    for (NodeIt v(_G); v != lemon::INVALID; ++v)
    {
      L.disable(v);
    }
    for (ArcIt a(_G); a != lemon::INVALID; ++a)
    {
      L.disable(a);
    }
    
    // report T and copy to L
    ArcList arcsL;
    L.enable(_root);
    for (SubArcIt a(T); a != lemon::INVALID; ++a)
    {
      arcsL.push_back(a);
      L.enable((Arc)a);
      L.enable((Node)T.target(a));
    }
    _result.push_back(arcsL);
    
    if (_limit == -1 || _result.size() < _limit)
      return true;
    else
      return false;
//    std::cerr << "\rNumber of trees of size " << arcsL.size() << ": " << _result.size() << std::flush;
  }
  else
  {
    ArcList FF;
    
    bool done;
    do
    {
      assert(!F.empty());
      
      Arc uv = F.back();
      F.pop_back();
      Node u = _G.source(uv);
      Node v = _G.target(uv);
      
      assert(T.status(u));
      assert(!T.status(v));
      assert(!T.status(uv));
      
      // add uv to T
      T.enable(uv);
      T.enable(v);

      ArcList newF = F;
      
      // push each arc vw where w not in V(T) onto F
      for (SubOutArcIt vw(G, v); vw != lemon::INVALID; ++vw)
      {
        Node w = G.target(vw);
        if (!T.status(w))
        {
          newF.push_back(vw);
        }
      }
      
      // remove each arc wv where w in T from F
      for (ArcListNonConstIt it = newF.begin(); it != newF.end();)
      {
        if (_G.target(*it) == v && T.status(_G.source(*it)))
        {
          it = newF.erase(it);
        }
        else
        {
          ++it;
        }
      }
      
      if (!grow(G, T, L, bfsL, newF))
        return false;
      
      G.disable(uv);
      T.disable(uv);
      T.disable(v);
      
      FF.push_back(uv);
      
      done = true;
      if (lemon::countNodes(L) != 0)
      {
        bfsL.run(v);
        
        for (SubInArcIt wv(G, v); wv != lemon::INVALID; ++wv)
        {
          Node w = G.source(wv);
          if (!bfsL.reached(w))
            done = false;
        }
      }
    } while (!done);
    
    for (ArcListRevIt it = FF.rbegin(); it != FF.rend(); ++it)
    {
      Arc a = *it;
      assert(!G.status(a));
      
      F.push_back(*it);
      G.enable(a);
    }
  }
  
  return true;
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
