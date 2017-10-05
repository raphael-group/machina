/*
 * solutiongraph.cpp
 *
 *  Created on: 21-aug-2017
 *      Author: M. El-Kebir
 */

#include "solutiongraph.h"
#include "gabowmyers.h"
#include <lemon/connectivity.h>
#include <fstream>

SolutionGraph::SolutionGraph(const CloneTree& T,
                             const IntPairNodeMap& leafLabelingT,
                             const IntPairToNodeMap& lcaT,
                             const IntPairSetNodeMap& vertexToStateSet,
                             const IntPairToNodeSetMap& X,
                             int primaryIndex,
                             const StringVector& indexToAnatomicalSite)
  : _T(T)
  , _leafLabelingT(leafLabelingT)
  , _lcaT(lcaT)
  , _vertexToStateSet(vertexToStateSet)
  , _X(X)
  , _primaryIndex(primaryIndex)
  , _indexToAnatomicalSite(indexToAnatomicalSite)
  , _G()
  , _rootG(lemon::INVALID)
  , _vertexLabelingG(_G)
  , _scToG()
  , _stateTrees()
  , _refinements()
{
  initG();
  
  VerbosityLevel oldVerbosity = g_verbosity;
  g_verbosity = VERBOSE_NONE;
  std::ofstream outG("/tmp/G.dot");
  writeDOT(outG);
  outG.close();
  enumerate();
  g_verbosity = oldVerbosity;
}

SolutionGraph::~SolutionGraph()
{
  for (SubTree* pS : _stateTrees)
  {
    delete pS;
  }
  
  for (Refinement* pRefinement : _refinements)
  {
    delete pRefinement;
  }
}

void SolutionGraph::writeDOT(std::ostream& out) const
{
  out << "digraph G {" << std::endl;
  for (NodeIt v_sc(_G); v_sc != lemon::INVALID; ++v_sc)
  {
    IntPair sc = _vertexLabelingG[v_sc];
    out << "\t" << _G.id(v_sc)
    << " [label=\"" << _indexToAnatomicalSite[sc.first]
    << " , " << sc.second << "\"]" << std::endl;
  }
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    out << "\t" << _G.id(_G.source(a)) << " -> " << _G.id(_G.target(a)) << std::endl;
  }
  out << "}" << std::endl;
}

void SolutionGraph::writeDOT(std::ostream& out,
                             const SubDigraph& S) const
{
  out << "digraph S {" << std::endl;
  for (SubNodeIt v_sc(S); v_sc != lemon::INVALID; ++v_sc)
  {
    IntPair sc = _vertexLabelingG[v_sc];
    out << "\t" << _G.id(v_sc)
    << " [label=\"" << _indexToAnatomicalSite[sc.first]
    << " , " << sc.second << "\"]" << std::endl;
  }
  for (SubArcIt a(S); a != lemon::INVALID; ++a)
  {
    out << "\t" << _G.id(_G.source(a)) << " -> " << _G.id(_G.target(a)) << std::endl;
  }
  out << "}" << std::endl;
}

void SolutionGraph::initG()
{
  std::set<IntPair> leafLabelSet;
  for (Node v_i : _T.leafSet())
  {
    leafLabelSet.insert(_leafLabelingT[v_i]);
  }
  
  // add vertices
  for (IntPair sc : leafLabelSet)
  {
    Node v_sc = _G.addNode();
    _vertexLabelingG[v_sc] = sc;
    _scToG[sc] = v_sc;
    
    if (sc.first == _primaryIndex && sc.second == 0)
    {
      _rootG = v_sc;
    }
  }
  
  // add edges
  for (NodeIt v_sc(_G); v_sc != lemon::INVALID; ++v_sc)
  {
    IntPair sc = _vertexLabelingG[v_sc];
    for (NodeIt v_td(_G); v_td != lemon::INVALID; ++v_td)
    {
      if (v_td == _rootG) continue;
      if (v_sc == v_td) continue;
      IntPair td = _vertexLabelingG[v_td];
      
      Node lca_sc = _lcaT[sc];
      Node lca_td = _lcaT[td];
      assert(lca_sc != lemon::INVALID);
      assert(lca_td != lemon::INVALID);
      if (_lcaT[sc] == _lcaT[td])
      {
        _G.addArc(v_sc, v_td);
      }
      else
      {
        // traverse from lca_td to lca_sc
        Node parent = _T.parent(lca_td);
        while (parent != lemon::INVALID)
        {
          if (parent == lca_sc)
          {
            _G.addArc(v_sc, v_td);
            break;
          }
          else if (!_vertexToStateSet[parent].empty())
          {
            assert(_vertexToStateSet[parent].count(sc) == 0);
            break;
          }
          parent = _T.parent(parent);
        }
      }
    }
  }
}

void SolutionGraph::enumerate()
{
  const int nrVertices = lemon::countNodes(_G);
  
  GabowMyers gm(_G, _rootG, 10);
  gm.run();
  
  int n = gm.getNrTrees();
  for (int i = 0; i < n; ++i)
  {
    SubTree* pSubTree = new SubTree(_G);
    
    gm.result(pSubTree->subT(), i);
    if (lemon::countArcs(pSubTree->subT()) == nrVertices - 1)
    {
      _stateTrees.push_back(pSubTree);
    }
    else
    {
      writeDOT(std::cerr);
      writeDOT(std::cerr, pSubTree->subT());
      assert(false);
      delete pSubTree;
    }
  }
  
  for (SubTree* pS : _stateTrees)
  {
    refine(pS->subT());
  }
}

void SolutionGraph::refine(const SubDigraph& S)
{
  // 0. Start checking whether overlap is indeed at most 1.
  std::set<IntPair> leafLabelSet;
  for (Node v_i : _T.leafSet())
  {
    leafLabelSet.insert(_leafLabelingT[v_i]);
  }
  
#ifdef DEBUG
  IntPairSetNodeMap invLabels(_T.tree());
  for (IntPair sc : leafLabelSet)
  {
    const NodeSet& X_sc = _X.find(sc)->second;
    for (IntPair td : leafLabelSet)
    {
      if (sc == td) continue;
      
      const NodeSet& X_td = _X.find(td)->second;
      NodeSet intersection;
      std::set_intersection(X_sc.begin(), X_sc.end(),
                            X_td.begin(), X_td.end(),
                            std::inserter(intersection, intersection.begin()));
      
      if (intersection.size() > 1)
      {
        std::cerr << "(s,c) = (" << _indexToAnatomicalSite[sc.first] << "," << sc.second << ") ; ";
        std::cerr << "X_sc =";
        for (Node v : X_sc)
        {
          std::cerr << " " << _T.label(v);
        }
        std::cerr << std::endl;
        std::cerr << "(t,d) = (" << _indexToAnatomicalSite[td.first] << "," << td.second << ") ; ";
        std::cerr << "X_td =";
        for (Node v : X_td)
        {
          std::cerr << " " << _T.label(v);
        }
        std::cerr << std::endl;
      }
      assert(intersection.size() <= 1);
    }
  }
#endif //DEBUG
  
  // 1. refine
  Digraph Tprime;
  Node root_Tprime = Tprime.addNode();
  StringNodeMap label(Tprime);
  StringNodeMap lPlus(Tprime);
  label[root_Tprime] = _T.label(_T.root());
  lPlus[root_Tprime] = _indexToAnatomicalSite[_vertexLabelingG[_rootG].first];

  refine(S, _T.root(), Tprime, root_Tprime, label, lPlus);
  
  // 2. construct clone tree
  Refinement* pRefinement = new Refinement(Tprime, root_Tprime, label, lPlus);
  _refinements.push_back(pRefinement);
}

bool SolutionGraph::isConnected(const SubDigraph& S,
                                const IntPairSet& Sigma_u,
                                IntPair& root_sc) const
{
  // Check whether Sigma_u induces a connected subgraph of S
  BoolNodeMap filter(_G, false);
  for (IntPair sc : Sigma_u)
  {
    assert(_scToG.count(sc) == 1);
    filter[_scToG.find(sc)->second] = true;
  }

  root_sc = std::make_pair(-1, -1);
  for (IntPair sc : Sigma_u)
  {
    assert(_scToG.count(sc) == 1);
    Node v_sc = _scToG.find(sc)->second;
    Arc a = InArcIt(_G, v_sc);
    if (a == lemon::INVALID)
    {
      root_sc = sc;
    }
    else if (!filter[_G.source(a)])
    {
      if (root_sc == std::make_pair(-1, -1))
      {
        root_sc = sc;
      }
      else
      {
        root_sc = std::make_pair(-1, -1);
        return false;
      }
    }
  }
  assert (root_sc != std::make_pair(-1, -1));
  return true;
}

bool SolutionGraph::areSiblings(const SubDigraph& S,
                                const IntPairSet& Sigma_u,
                                IntPair& parent_sc) const
{
  parent_sc = std::make_pair(-1, -1);
  for (IntPair sc : Sigma_u)
  {
    assert(_scToG.count(sc) == 1);
    Node v_sc = _scToG.find(sc)->second;
    
    if (v_sc == _rootG)
    {
      return false;
    }
    
    Node v_dt = _G.source(InArcIt(_G, v_sc));
    if (parent_sc == std::make_pair(-1, -1))
    {
      parent_sc = _vertexLabelingG[v_dt];
    }
    else if (parent_sc != _vertexLabelingG[v_dt])
    {
      parent_sc = std::make_pair(-1, -1);
      return false;
    }
  }
  
  return true;
}

void SolutionGraph::refine(const SubDigraph& S,
                           Node u,
                           Digraph& Tprime,
                           Node uu,
                           StringNodeMap& label,
                           StringNodeMap& lPlus)
{
  const IntPairSet& Sigma_u = _vertexToStateSet[u];
  if (Sigma_u.size() <= 1)
  {
    // leave vertex u unperturbed
    if (_T.isLeaf(u))
    {
      lPlus[uu] = _indexToAnatomicalSite[_leafLabelingT[u].first];
    }
    else
    {
      if (u != _T.root())
      {
        if (Sigma_u.empty())
        {
          Node pi_uu = Tprime.source(InArcIt(Tprime, uu));
          lPlus[uu] = lPlus[pi_uu];
        }
        else
        {
          IntPair sc = *(Sigma_u.begin());
          lPlus[uu] = _indexToAnatomicalSite[sc.first];
        }
      }
      for (OutArcIt a(_T.tree(), u); a != lemon::INVALID; ++a)
      {
        Node v = _T.tree().target(a);
        Node vv = Tprime.addNode();
        Tprime.addArc(uu, vv);
        
        label[vv] = _T.label(v);
        refine(S, v, Tprime, vv, label, lPlus);
      }
    }
  }
  else
  {
    // Check whether Sigma_u induces a connected subgraph of S
    IntPair pi_sc;
    IntPairToNodeMap toBackBone;
    if (isConnected(S, Sigma_u, pi_sc))
    {
      // add vertices
      for (IntPair sc : Sigma_u)
      {
        if (sc != pi_sc)
        {
          Node v_sc = Tprime.addNode();
          toBackBone[sc] = v_sc;
          label[v_sc] = _T.label(u) + "_" + _indexToAnatomicalSite[sc.first];
          lPlus[v_sc] = _indexToAnatomicalSite[sc.first];
        }
        else
        {
          toBackBone[sc] = uu;
          lPlus[uu] = _indexToAnatomicalSite[sc.first];
        }
      }
      
      // add edges
      for (IntPair sc : Sigma_u)
      {
        if (sc != pi_sc)
        {
          Node v_sc = toBackBone[sc];
          
          IntPair td = _vertexLabelingG[_G.source(SubInArcIt(S, _scToG[sc]))];
          Node v_td = toBackBone[td];
          Tprime.addArc(v_td, v_sc);
        }
      }
    }
    else if (areSiblings(S, Sigma_u, pi_sc))
    {
      lPlus[uu] = _indexToAnatomicalSite[pi_sc.first];

      // add vertices
      for (IntPair sc : Sigma_u)
      {
        Node v_sc = Tprime.addNode();
        toBackBone[sc] = v_sc;
        label[v_sc] = _T.label(u) + "_" + _indexToAnatomicalSite[sc.first];
        lPlus[v_sc] = _indexToAnatomicalSite[sc.first];
        Tprime.addArc(uu, v_sc);
      }
    }
    else
    {
      writeDOT(std::cout);
      writeDOT(std::cout, S);
      for (IntPair sc : Sigma_u)
      {
        std::cout << _indexToAnatomicalSite[sc.first] << " , " << sc.second << std::endl;
      }
      assert(false);
    }
    
    // add original children
    for (OutArcIt a(_T.tree(), u); a != lemon::INVALID; ++a)
    {
      Node v = _T.tree().target(a);
      Node vv = Tprime.addNode();
      
      label[vv] = _T.label(v);
      IntPairSet intersection;
      std::set_intersection(Sigma_u.begin(), Sigma_u.end(),
                            _vertexToStateSet[v].begin(), _vertexToStateSet[v].end(),
                            std::inserter(intersection, intersection.begin()));
      
      if (intersection.size() == 1)
      {
        IntPair sc = *intersection.begin();
        Tprime.addArc(toBackBone[sc], vv);
      }
      else
      {
        assert(intersection.size() == 0);
        Tprime.addArc(uu, vv);
      }
      
      refine(S, v, Tprime, vv, label, lPlus);
    }
  }
}
