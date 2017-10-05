/*
 * rootedcladisticenumeration.cpp
 *
 *  Created on: 29-sep-2015
 *      Author: M. El-Kebir
 */

#include "rootedcladisticenumeration.h"
#include <lemon/bfs.h>

namespace gm {

RootedCladisticEnumeration::RootedCladisticEnumeration(const RootedCladisticAncestryGraph& G,
                                                       int limit,
                                                       int timeLimit,
                                                       int threads,
                                                       int lowerbound,
                                                       bool monoclonal,
                                                       bool fixTrunk,
                                                       const IntSet& whiteList)
  : _G(G)
  , _result()
  , _objectiveValue(0)
  , _limit(limit)
  , _timeLimit(timeLimit)
  , _threads(threads)
  , _lowerbound(lowerbound)
  , _counter(0)
  , _mutex()
  , _sem(threads)
  , _threadGroup()
  , _timer()
  , _monoclonal(monoclonal)
  , _fixTrunk(fixTrunk)
  , _whiteList(whiteList)
{
}
  
RootedCladisticEnumeration::~RootedCladisticEnumeration()
{
}

std::string RootedCladisticEnumeration::newick(int solIdx) const
{
  assert(0 <= solIdx && solIdx < _result.size());
  
  ArcListList::const_iterator it = _result.begin();
  int tmp = solIdx;
  for (; tmp > 0; tmp--) ++it;
  
  const ArcList& L = *it;
  const Digraph& G = _G.G();
  
  BoolNodeMap filterNodesT(G, false);
  BoolArcMap filterArcsT(G, false);
  SubDigraph T(G, filterNodesT, filterArcsT);
  
  for (ArcListIt it = L.begin(); it != L.end(); ++it)
  {
    Arc a_cidj = *it;
    Node v_ci = G.source(a_cidj);
    Node v_dj = G.target(a_cidj);
    
    T.enable(a_cidj);
    T.enable(v_ci);
    T.enable(v_dj);
  }
  
  std::string str = "[";
  newick(T, _G.root(), str);
  str += "]";
  
  return str;
}
  
void RootedCladisticEnumeration::newick(const SubDigraph& T,
                                        Node v_ci,
                                        std::string& str) const
{
  str += _G.label(v_ci);
  
  bool children = SubOutArcIt(T, v_ci) != lemon::INVALID;
  if (children)
  {
    str += ";[";
    
    bool first = true;
    for (SubOutArcIt a_cidj(T, v_ci); a_cidj != lemon::INVALID; ++a_cidj)
    {
      if (first)
        first = false;
      else
        str += ";";
      
      Node v_dj = T.target(a_cidj);
      newick(T, v_dj, str);
    }
    str += "]";
  }
}
  
bool RootedCladisticEnumeration::prune(const SubDigraph& T,
                                       const ArcList& F) const
{
  int size_L = 0;
  {
    boost::interprocess::scoped_lock<boost::mutex> lock(_mutex);
    size_L = _lowerbound;
  }
  
  const Digraph& G = _G.G();
  BoolNodeMap filterNodesTT(G, true);
  BoolArcMap filterArcsTT(G, true);
  
  SubDigraph TT(G, filterNodesTT, filterArcsTT);
  
  // Remove arcs incoming to vertices in T and F
  BoolArcMap whiteList(G, false);  // these are arcs that occur in F, and need to be retained in TT
  BoolNodeMap blackList(G, false); // these are nodes that are the target nodes of an arc in F
  for (ArcListIt it = F.begin(); it != F.end(); ++it)
  {
    Arc a_cidj = *it;
    whiteList[a_cidj] = true;
    
    Node v_dj = G.target(a_cidj);
    blackList[v_dj] = true;
  }
  
  // determine set C of characters present in the tree (disregarding root)
  IntSet C = _whiteList; // initialize to whitelist, these characters must be present in the tree
  for (SubNodeIt v_ci(T); v_ci != lemon::INVALID; ++v_ci)
  {
    if (v_ci != _G.root())
    {
      for (const IntPair& ci : _G.nodeToCharState(v_ci))
      {
        C.insert(ci.first);
      }
    }
  }

  //
  for (ArcIt a_cidj(G); a_cidj != lemon::INVALID; ++a_cidj)
  {
    // retain arcs in T and F
    if (T.status(a_cidj) || whiteList[a_cidj])
      continue;
    
    Node v_ci = G.source(a_cidj);
    Node v_dj = G.target(a_cidj);
    
    // pre: a_cidj not in F
    // pre: a_cidj not in T
    if (T.status(v_ci))
    {
      // disable arcs incoming to vertices in T
      TT.disable(a_cidj);
    }
    if (blackList[v_dj])
    {
      // disable arcs incoming to target vertices in F
      TT.disable(a_cidj);
    }
  }
  
  // now do a BFS in TT
  assert(TT.status(_G.root()));
  lemon::Bfs<SubDigraph> bfs(TT);
  bfs.run(_G.root());
  
  for (NodeIt v_ci(_G.G()); v_ci != lemon::INVALID; ++v_ci)
  {
    if (!bfs.reached(v_ci))
      TT.disable(v_ci);
  }
  
  // check if there is a character that is state incomplete
  for (IntSetIt it = C.begin(); it != C.end(); ++it)
  {
    int c = *it;
    if (!isStateComplete(TT, c))
    {
      return true;
    }
  }
  
  int size_TT = makeStateComplete(TT);
  if (size_TT < size_L)
  {
//    std::cout << "TT: " << size_TT << ", L: " << size_L << ", F: " << F.size() << std::endl;
    return true;
  }
  
  return false;
}
  
void RootedCladisticEnumeration::init(Arc a_00dj,
                                      SubDigraph& subG,
                                      SubDigraph& T,
                                      ArcList& F)
{
  const Digraph& G = _G.G();
  Node root = _G.root();

  assert(G.source(a_00dj) == root);
  Node v_dj = G.target(a_00dj);
  
  T.enable(root);
  T.enable(v_dj);
  T.enable(a_00dj);
  
  F.clear();
  
  // don't add arcs that other threads will consider
  OutArcIt a(G, root);
  while (a != a_00dj)
  {
//    _mutex.lock();
//    std::cerr << "Disabling " << _G.label(G.target(a)) << " [" << _G.label(v_dj) << "]" << std::endl;
//    _mutex.unlock();
    subG.disable(a);
    ++a;
  }

  if (!_monoclonal)
  {
    for (; a != lemon::INVALID; ++a)
    {
      if (a != a_00dj && isValid(T, a))
        F.push_back(a);
    }
  }
  else
  {
    // disable all other arcs if monoclonal
    for (; a != lemon::INVALID; ++a)
    {
      subG.disable(a);
    }
  }
  
  for (OutArcIt a(G, v_dj); a != lemon::INVALID; ++a)
  {
    if (isValid(T, a))
      F.push_back(a);
  }
  
  F.sort(Compare(_G));
  assert(isValid(T));
}
  
void RootedCladisticEnumeration::init(SubDigraph& subG,
                                      SubDigraph& T,
                                      ArcList& F)
{
  const Digraph& G = _G.G();
  Node root = _G.root();

  T.enable(root);
  F.clear();
  for (OutArcIt a(G, root); a != lemon::INVALID; ++a)
  {
    F.push_back(a);
  }
  
  F.sort(Compare(_G));
}
  
void RootedCladisticEnumeration::runArc(Arc a_00dj)
{
  _sem.wait();
  const Digraph& G = _G.G();
  
  BoolNodeMap filterNodesT(G, false);
  BoolArcMap filterArcsT(G, false);
  SubDigraph T(G, filterNodesT, filterArcsT);
  
  BoolNodeMap filterNodesG(G, true);
  BoolArcMap filterArcsG(G, true);
  SubDigraph subG(G, filterNodesG, filterArcsG);
  
  ArcList F;
  
  init(a_00dj, subG, T, F);
  grow(subG, T, F);
  
//  for (int idx = 0; idx < _result.size(); ++idx)
//  {
//    std::cerr << newick(idx) << " " << idx << std::endl;
//  }
  _sem.post();
}
  
void RootedCladisticEnumeration::run()
{
  const Digraph& G = _G.G();
  Node root = _G.root();
  
  _timer.start();
  _counter = 0;
  if (_threads == 1)
  {
    BoolNodeMap filterNodesT(G, false);
    BoolArcMap filterArcsT(G, false);
    SubDigraph T(G, filterNodesT, filterArcsT);
    
    BoolNodeMap filterNodesG(G, true);
    BoolArcMap filterArcsG(G, true);
    SubDigraph subG(G, filterNodesG, filterArcsG);
    
    ArcList F;
    
    init(subG, T, F);
    grow(subG, T, F);
  }
  else
  {
    for (OutArcIt a_00dj(G, root); a_00dj != lemon::INVALID; ++a_00dj)
    {
//      std::cerr << _G.label(G.target(a_00dj)) << std::endl;
      _threadGroup.create_thread(boost::bind(&RootedCladisticEnumeration::runArc, this, a_00dj));
//      RootedCladisticEnumeration::runArc(a_00dj);
    }
    
    _threadGroup.join_all();
  }

  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cerr << "\r" << std::flush;
  }
}
  
bool RootedCladisticEnumeration::finalize(SubDigraph& T)
{
  assert(isValid(T));
  
//  writeDOT(std::cout, T, ArcList());
  
  boost::interprocess::scoped_lock<boost::mutex> lock(_mutex);
  int currentSize = std::max(_lowerbound, _objectiveValue);
  
  // make copy of T
  BoolNodeMap filterNodesT(_G.G(), false);
  BoolArcMap filterArcsT(_G.G(), false);
  SubDigraph TT(_G.G(), filterNodesT, filterArcsT);
  for (NodeIt v_ci(_G.G()); v_ci != lemon::INVALID; ++v_ci)
  {
    TT.status(v_ci, T.status(v_ci));
  }
  for (ArcIt a_cidj(_G.G()); a_cidj != lemon::INVALID; ++ a_cidj)
  {
    TT.status(a_cidj, T.status(a_cidj));
  }
  
  int newSizeT = makeStateComplete(TT);
  if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
  {
    std::cerr << "\r                                                            \r";
    std::cerr << newSizeT << "/" << _result.size() << "/" << _counter << "/" << _limit << " (" << _lowerbound << ")" << std::flush;
  }
  
  const int k = _G.F().k();
  for (int c : _whiteList)
  {
    for (int i = 1; i < k; ++i)
    {
      // abort if the tree misses a character in the whitelist
      if (_G.S(c).isPresent(i) && !TT.status(_G.charStateToNode(c, i)))
      {
        return false;
      }
    }
  }
  
  if (newSizeT < currentSize || newSizeT == 0)
  {
    return false;
  }
  else if (newSizeT > currentSize && !_result.empty())
  {
    _objectiveValue = _lowerbound = newSizeT;
    if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
    {
      std::cerr << std::endl;
    }
    _result.clear();
  }
  else
  {
    _objectiveValue = newSizeT;
  }
  
  ++_counter;
  
  _result.push_back(ArcList());
  ArcList& res = _result.back();
  
  for (SubArcIt a_cidj(TT); a_cidj != lemon::INVALID; ++a_cidj)
  {
    res.push_back(a_cidj);
  }
  
  return _result.size() >= _limit;
}
  
bool RootedCladisticEnumeration::grow(SubDigraph& G,
                                      SubDigraph& T,
                                      ArcList& F)
{
  if (limitReached())
  {
    return true;
  }
  else if (F.empty())
  {
    return finalize(T);
  }
  else if (prune(T, F))
  {
    return false;
  }
  else
  {
    ArcList FF;
    
    do
    {
      assert(!F.empty());
      
      Arc a_cidj = F.back();
      F.pop_back();
      
//      Arc a_cidj = F.front();
//      F.pop_front();
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
        else if (G.source(*it) == v_ci && !isValid(T, *it))
        {
          it = newF.erase(it);
        }
        else
        {
          assert(isValid(T, *it));
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
        
        // isFirstAncestor violation (consistency)
        bool isConsistent = true;
        for (const IntPair& el : _G.nodeToCharState(v_el))
        {
          int pi_l = _G.S(el.first).parent(el.second);
          assert(pi_l != -1 && pi_l != -2);
          if (!(pi_l != -1 && pi_l != -2))
          {
            abort();
          }
          Node v_e_pi_l = _G.charStateToNode(el.first, pi_l);
          assert(v_e_pi_l != lemon::INVALID);
          
          isConsistent = isConsistent && isFirstAncestor(T, el.first, v_e_pi_l, v_dj);
        }
        
        if (isValid(T, a_djel))
        {
          newF.push_back(a_djel);
        }
      }
      
      if (grow(G, T, newF))
        return true;
      
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
    
    return false;
  }
}
  
bool RootedCladisticEnumeration::isFirstAncestor(const SubDigraph& T,
                                                 int c, // character to check
                                                 Node v_ci, // is this node the first ancestor
                                                 Node v_dj) const // starting from this node
{
  assert(T.status(v_dj));
  
  if (!T.status(v_ci))
  {
    return false;
  }
  if (v_ci == v_dj)
  {
    return true;
  }
  
  Arc a(SubInArcIt(T, v_dj));
  while (a != lemon::INVALID)
  {
    v_dj = T.source(a);
    a = SubInArcIt(T, v_dj);
    
    const auto& X_dj = _G.nodeToCharState(v_dj);
    bool ok = true;
    for (const IntPair& dj : X_dj)
    {
      if (dj.first == c && v_ci != v_dj)
      {
        ok = false; // not sure, maybe there is another dj in X_dj for which the condition holds
      }
      else if (dj.first == c)
      {
        return true; // found it
      }
    }
    if (!ok)
    {
      return false;
    }
  }
  
  assert(v_ci != v_dj);
  return v_ci == v_dj;
}
  
bool RootedCladisticEnumeration::isAncestor(const SubDigraph& T,
                                            Node v_dj,
                                            Node v_ci) const
{
  assert(T.status(v_ci));
  if (!T.status(v_dj))
  {
    return false;
  }
  
  Node root = _G.root();
  while (v_ci != root && v_ci != v_dj)
  {
    Arc a(SubInArcIt(T, v_ci));
    v_ci = T.source(a);
  }
  
  return v_ci == v_dj;
}
  
void RootedCladisticEnumeration::addArc(SubDigraph& T,
                                        Arc a_cidj) const
{
  const Node v_dj = T.target(a_cidj);
  
  // add a_cidj to T
  T.enable(v_dj);
  T.enable(a_cidj);
  
  assert(isArborescence(T));
  assert(isValid(T));
}
  
void RootedCladisticEnumeration::removeArc(SubDigraph& T,
                                           Arc a_cidj) const
{
  assert(T.status(a_cidj));
  
  const Node v_dj = T.target(a_cidj);
  
  // remove a_cidj from T
  T.disable(a_cidj);
  T.disable(v_dj);
  
  assert(isArborescence(T));
  assert(isValid(T));
}
  
bool RootedCladisticEnumeration::isArborescence(const SubDigraph& T) const
{
  assert(T.status(_G.root()));
  lemon::Bfs<SubDigraph> bfs(T);
  bfs.run(_G.root());
  
  for (SubNodeIt v_ci(T); v_ci != lemon::INVALID; ++v_ci)
  {
    if (!bfs.reached(v_ci))
    {
      return false;
    }
  }
  
  return true;
}
  
bool RootedCladisticEnumeration::isValid(const SubDigraph& T) const
{
  const RealTensor& F = _G.F();
  const int m = F.m();

  for (int p = 0; p < m; ++p)
  {
    for (SubNodeIt v_ci(T); v_ci != lemon::INVALID; ++v_ci)
    {
      for (const IntPair& ci : _G.nodeToCharState(v_ci))
      {
        double f_p_ci = F.getCumFreq(p, ci.first, _G.S(ci.first).D(ci.second));
        
        double sum_of_children = 0;
        
        for (SubOutArcIt a_cidj(T, v_ci); a_cidj != lemon::INVALID; ++a_cidj)
        {
          Node v_dj = T.target(a_cidj);
          
          // find maximum
          IntPair max_dj(-1, -1);
          double max_f_cum_p_dj = -1;
          
          for (const IntPair& dj : _G.nodeToCharState(v_dj))
          {
            double f_cum_p_dj = F.getCumFreq(p, dj.first, _G.S(dj.first).D(dj.second));
            if (f_cum_p_dj > max_f_cum_p_dj)
            {
              max_f_cum_p_dj = f_cum_p_dj;
              max_dj = dj;
            }
          }
          
          sum_of_children += max_f_cum_p_dj;
        }
       
        if (g_tol.less(f_p_ci, sum_of_children))
        {
          return false;
        }
      }
    }
  }
  
  return true;
}

bool RootedCladisticEnumeration::isValid(const SubDigraph& T,
                                         Arc a_ciel) const
{
  const Node v_ci = T.source(a_ciel);
  const RealTensor& F = _G.F();
  const int m = F.m();
  Node v_el = T.target(a_ciel);

  for (const IntPair& el : _G.nodeToCharState(v_el))
  {
    for (int p = 0; p < m; ++p)
    {
      for (const IntPair& ci : _G.nodeToCharState(v_ci))
      {
        double f_p_ci = F.getCumFreq(p, ci.first, _G.S(ci.first).D(ci.second));
        
        double sum_of_children = F.getCumFreq(p, el.first, _G.S(el.first).D(el.second));
        for (SubOutArcIt a_cidj(T, v_ci); a_cidj != lemon::INVALID; ++a_cidj)
        {
          Node v_dj = T.target(a_cidj);
          
          // find maximum
          IntPair max_dj(-1, -1);
          double max_f_cum_p_dj = -1;
          
          for (const IntPair& dj : _G.nodeToCharState(v_dj))
          {
            double f_cum_p_dj = F.getCumFreq(p, dj.first, _G.S(dj.first).D(dj.second));
            if (f_cum_p_dj > max_f_cum_p_dj)
            {
              max_f_cum_p_dj = f_cum_p_dj;
              max_dj = dj;
            }
          }
          
          sum_of_children += max_f_cum_p_dj;
        }
        
        if (g_tol.less(f_p_ci, sum_of_children))
        {
          return false;
        }
      }
    }
  }
  
  return true;
}
  
Solution RootedCladisticEnumeration::solution(int solIdx) const
{
  assert(0 <= solIdx && solIdx < _result.size());
  const int n = _G.F().n();
  const Digraph& G = _G.G();
  
  ArcListList::const_iterator it = _result.begin();
  int tmp = solIdx;
  for (; tmp > 0; tmp--) ++it;
  
  BoolNodeMap nodes(G, false);
  BoolArcMap arcs(G, false);
  SubDigraph T(G, nodes, arcs);
  T.enable(_G.root());
  
  const ArcList& arcList = *it;
  for (ArcListIt it2 = arcList.begin(); it2 != arcList.end(); ++it2)
  {
    Arc a = *it2;
    T.enable(G.source(a));
    T.enable(G.target(a));
    T.enable(a);
  }
  
  assert(T.status(_G.root()));
  
  BoolVector stateComplete(n, false);
  for (int c = 0; c < n; ++c)
  {
    stateComplete[c] = isStateComplete(T, c);
  }
  
  Solution sol(_G.F(), _G.S(), _G.F().m(), _G.F().n(), _G.F().k());
  
//  writeDOT(std::cerr, T, ArcList());
  
  initA(T, _G.root(), stateComplete, sol._A);
  initF(solIdx, sol._inferredF);
  initU(T, sol._inferredF, _G.root(), stateComplete, sol._U);
  
  return sol;
}
  
void RootedCladisticEnumeration::writeDOT(std::ostream& out,
                                          const SubDigraph& T,
                                          const ArcList& H) const
{
  const Digraph& G = _G.G();
  const int m = _G.F().m();
  
  out << "digraph T {" << std::endl;
  out.precision(3);
  
  for (SubNodeIt v_ci(T); v_ci != lemon::INVALID; ++v_ci)
  {
    const auto& X_ci = _G.nodeToCharState(v_ci);
    
    out << "\t" << G.id(v_ci) << " [label=\"" << _G.label(v_ci) << "\\n";
    
    for (const IntPair& ci : X_ci)
    {
      const IntSet& D_ci = _G.S(ci.first).D(ci.second);
      out << "{";
      for (IntSetIt it = D_ci.begin(); it != D_ci.end(); ++it)
      {
        out << " " << *it;
      }
      out << " } ";
    }
    out << "\\n";
    
    for (int p = 0; p < m; ++p)
    {
      for (const IntPair& ci : X_ci)
      {
        out << _G.F().getCumFreq(p, ci.first, _G.S(ci.first).D(ci.second)) << " ";
      }
      out << "\\n";
    }
    out << "\"]" << std::endl;
  }
  
  for (SubArcIt a(T); a != lemon::INVALID; ++a)
  {
    out << "\t" << G.id(G.source(a)) << " -> " << G.id(G.target(a)) << std::endl;
  }
  
  for (ArcListIt it = H.begin(); it != H.end(); ++it)
  {
    Arc a = *it;
    Node v_ci = G.source(a);
    Node v_dj = G.target(a);
    
    const auto& X_dj = _G.nodeToCharState(v_dj);
    
    out << "\t" << G.id(v_dj) << " [label=\"" << _G.label(v_dj) << "\\n";
    
    for (const IntPair& dj : X_dj)
    {
      const IntSet& D_dj = _G.S(dj.first).D(dj.second);
      out << "{";
      for (IntSetIt it = D_dj.begin(); it != D_dj.end(); ++it)
      {
        out << " " << *it;
      }
      out << " } ";
    }
    out << "\\n";
    
    for (int p = 0; p < m; ++p)
    {
      for (const IntPair& dj : X_dj)
      {
        out << _G.F().getCumFreq(p, dj.first, _G.S(dj.first).D(dj.second)) << " ";
      }
      out << "\\n";
    }
    out << "\"]" << std::endl;

    out << "\t" << G.id(v_ci) << " -> " << G.id(v_dj) << " [color=red]" << std::endl;
  }
  
  out << "}" << std::endl;
}
  
void RootedCladisticEnumeration::initA(const SubDigraph& T,
                                       Node v_ci,
                                       const BoolVector& stateComplete,
                                       PerfectPhyloMatrix& A) const
{
  assert(T.status(v_ci));
  
  const int n = _G.F().n();

  if (v_ci == _G.root())
  {
    // initialize (*,0)
    for (int d = 0; d < n; ++d)
    {
      A.set(0, 0, d, 0);
    }
  }
  else
  {
    Node v_pi_ci = T.source(SubInArcIt(T, v_ci));
    IntPair pi_ci(-1, -1);
    if (v_pi_ci == _G.root())
    {
      pi_ci = IntPair(0, 0);
    }
    else
    {
      const auto& X_pi_ci = _G.nodeToCharStateList(v_pi_ci);
      for (auto it = X_pi_ci.rbegin(); it != X_pi_ci.rend(); ++it)
      {
        if (stateComplete[it->first])
        {
          pi_ci = *it;
          break;
        }
      }
      // T is state complete => there is no empty node
      assert(!(pi_ci.first == -1 && pi_ci.second == -1));
    }
    
    const auto& X_ci = _G.nodeToCharStateList(v_ci);
    for (const IntPair& ci : X_ci)
    {
      if (!stateComplete[ci.first])
        continue;
      
      for (int d = 0; d < n; ++d)
      {
        if (d != ci.first)
          A.set(ci, d, A(pi_ci, d));
        else
          A.set(ci, ci);
      }
      pi_ci = ci;
    }
  }
  
  for (SubOutArcIt a(T, v_ci); a != lemon::INVALID; ++a)
  {
    Node v_dj = T.target(a);
    initA(T, v_dj, stateComplete, A);
  }
}
  
void RootedCladisticEnumeration::initU(const SubDigraph& T,
                                       const RealTensor& F,
                                       Node v_ci,
                                       const BoolVector& stateComplete,
                                       RealMatrix& U) const
{
  const int m = F.m();
  const int n = F.n();
  
  if (v_ci == _G.root())
  {
    for (int p = 0; p < m; ++p)
    {
      double u_pci = 1;
      
      for (SubOutArcIt a(T, v_ci); a != lemon::INVALID; ++a)
      {
        Node v_dj = T.target(a);
        
        const auto& X_dj = _G.nodeToCharStateList(v_dj);
        for (auto it = X_dj.begin(); it != X_dj.end(); ++it)
        {
          if (stateComplete[it->first])
          {
            const IntPair& dj = *it;
            u_pci -= F.getCumFreq(p, dj.first, _G.S(dj.first).D(dj.second));
            break;
          }
        }
      }
      
      U.set(p, 0, u_pci);
    }
  }
  else
  {
    const auto& X_ci = _G.nodeToCharStateList(v_ci);
    IntPair last(-1, -1);
    for (auto it = X_ci.rbegin(); it != X_ci.rend(); ++it)
    {
      if (stateComplete[it->first])
      {
        last = *it;
        break;
      }
    }
    assert(!(last.first == -1 && last.second == -1));
    
    for (auto it_ci = X_ci.begin(); it_ci != X_ci.end(); ++it_ci)
    {
      const IntPair& ci = *it_ci;
      if (ci == last)
      {
        // deal with children
        for (int p = 0; p < m; ++p)
        {
          double u_pci = F.getCumFreq(p, ci.first, _G.S(ci.first).D(ci.second));
          for (SubOutArcIt a(T, v_ci); a != lemon::INVALID; ++a)
          {
            Node v_dj = T.target(a);
            
            const auto& X_dj = _G.nodeToCharStateList(v_dj);
            for (auto it = X_dj.begin(); it != X_dj.end(); ++it)
            {
              if (stateComplete[it->first])
              {
                const IntPair& dj = *it;
                u_pci -= F.getCumFreq(p, dj.first, _G.S(dj.first).D(dj.second));
                break;
              }
            }
          }
          
          U.set(p, n * (ci.second - 1) + ci.first + 1, u_pci);
        }
      }
      else if (stateComplete[ci.first])
      {
        auto it_dj = it_ci;
        IntPair dj = *(++it_dj);
        while (!stateComplete[dj.first])
          dj = *(++it_dj);
        
        for (int p = 0; p < m; ++p)
        {
          double u_pci = F.getCumFreq(p, ci.first, _G.S(ci.first).D(ci.second));
          u_pci -= F.getCumFreq(p, dj.first, _G.S(dj.first).D(dj.second));
          
          U.set(p, n * (ci.second - 1) + ci.first + 1, u_pci);
        }
      }
    }
  }
  
  for (SubOutArcIt a(T, v_ci); a != lemon::INVALID; ++a)
  {
    Node v_dj = T.target(a);
    initU(T, F, v_dj, stateComplete, U);
  }
}
  
void RootedCladisticEnumeration::populateSolutionSet(SolutionSet& sols) const
{
  boost::interprocess::scoped_lock<boost::mutex> lock(_mutex);
  
//  SolutionSet localSols;
  int solCount = solutionCount();
  for (int idx = 0; idx < solCount; ++idx)
  {
//    localSols.add(solution(idx));
    sols.add(solution(idx));
//    std::cerr << newick(idx) << " " << idx << std::endl;
  }
//  localSols.unique();
  
//  sols.add(localSols);
}
  
} // namespace gm
