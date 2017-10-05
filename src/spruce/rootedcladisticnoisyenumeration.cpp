/*
 * rootedcladisticnoisyenumeration.cpp
 *
 *  Created on: 11-oct-2015
 *      Author: M. El-Kebir
 */

#include "rootedcladisticnoisyenumeration.h"

namespace gm {

RootedCladisticNoisyEnumeration::RootedCladisticNoisyEnumeration(const RootedCladisticNoisyAncestryGraph& G,
                                                                 int limit,
                                                                 int timeLimit,
                                                                 int threads,
                                                                 int lowerbound,
                                                                 bool monoclonal,
                                                                 bool fixTrunk,
                                                                 const IntSet& whiteList)
  : RootedCladisticEnumeration(G, limit, timeLimit, threads,
                               lowerbound, monoclonal, fixTrunk, whiteList)
  , _noisyG(G)
{
}
  
void RootedCladisticNoisyEnumeration::init(Arc a_cidj,
                                           SubDigraph& subG,
                                           SubDigraph& T,
                                           ArcList& H,
                                           RealTensor& Fhat)
{
  const Digraph& G = _G.G();
  Node root = _G.root();
  Node v_ci = G.source(a_cidj);
  Node v_dj = G.target(a_cidj);

  // enable a_00ci and a_cidj
  Arc a = a_cidj;
  while (a != lemon::INVALID)
  {
    Node v = G.target(a);
    T.enable(a);
    T.enable(v);
    
    a = InArcIt(G, G.source(a));
  }
  T.enable(root);
  
  H.clear();
  Fhat = _noisyG.F_lb();
  
  // something weird
  if (!updateFhat(T, v_dj, Fhat))
    return;
  
  // don't add arcs that other threads will consider
  OutArcIt aa(G, v_ci);
  while (aa != a_cidj)
  {
    subG.disable(aa);
    ++aa;
  }
  
  if (!_monoclonal)
  {
    for (; aa != lemon::INVALID; ++aa)
    {
      if (aa != a_cidj && isValid(T, aa))
        H.push_back(aa);
    }
  }
//  else
//  {
//    // disable all other arcs if monoclonal
//    for (; a != lemon::INVALID; ++a)
//    {
//      subG.disable(a);
//    }
//  }
  
  // now update the frontier with arcs outgoing from v_dj
  for (OutArcIt a_djel(G, v_dj); a_djel != lemon::INVALID; ++a_djel)
  {
    Node v_el = G.target(a_djel);
    
    // isFirstAncestor violation (consistency)
    bool isConsistent = true;
    for (const IntPair& el : _G.nodeToCharState(v_el))
    {
      int pi_l = _G.S(el.first).parent(el.second);
      assert(0 <= pi_l && pi_l < Fhat.k());
      Node v_e_pi_l = _G.charStateToNode(el.first, pi_l);
      
      isConsistent = isConsistent && isFirstAncestor(T, el.first, v_e_pi_l, v_dj);
    }
    
    if (!isConsistent)
      continue;
    
    if (checkFhat(T, Fhat, a_djel))
    {
      assert(isValid(T, a_djel));
      H.push_back(a_djel);
    }
  }
  
  H.sort(Compare(_G));
  assert(isValid(T));
}
  
bool RootedCladisticNoisyEnumeration::checkFhat(SubDigraph& T,
                                                RealTensor& Fhat,
                                                Arc a_cidj) const
{
  const Node v_ci = T.source(a_cidj);
  const Node v_dj = T.target(a_cidj);
  
  assert(T.status(v_ci));
  assert(!T.status(v_dj));
  assert(!T.status(a_cidj));
  
  // add a_cidj to T
  T.enable(v_dj);
  T.enable(a_cidj);

#ifdef DEBUG
  RealTensor FhatCpy = Fhat;
#endif
  
  Fhat.setTrackChanges(true);
  bool res = updateFhat(T, v_dj, Fhat);
  Fhat.setTrackChanges(false);
  Fhat.rollBack();
  
#ifdef DEBUG
  assert(FhatCpy == Fhat);
#endif

  T.disable(a_cidj);
  T.disable(v_dj);

  return res;
}
  
bool RootedCladisticNoisyEnumeration::updateFhat(const SubDigraph& T,
                                                 Node v_ci,
                                                 RealTensor& Fhat) const
{
  assert(T.status(v_ci));
  
  const int m = Fhat.m();
  const int n = Fhat.n();
  const int k = Fhat.k();
  
  if (v_ci == _G.root())
  {
    // set \hat{f}_{p,(c,i)}
    for (int p = 0; p < m; ++p)
    {
      double sum_of_children = 0;
      for (SubOutArcIt a_cidj(T, v_ci); a_cidj != lemon::INVALID; ++a_cidj)
      {
        Node v_dj = T.target(a_cidj);
        const IntPair& dj = *(_G.nodeToCharState(v_dj).begin());
        
        double f_cum_hat_p_dj = Fhat.getCumFreq(p, dj.first, _G.S(dj.first).D(dj.second));
        sum_of_children += f_cum_hat_p_dj;
      }
      
      for (int c = 0; c < n; ++c)
      {
        double sum_of_descendant_states = 0; // proper descendants
        for (int j = 1; j < k; ++j)
        {
          if (_G.charStateToNode(c, j) != lemon::INVALID)
          {
            sum_of_descendant_states += Fhat(j, p, c);
          }
        }
        
        double f_hat_p_ci = 1 - sum_of_descendant_states;
        if (g_tol.less(f_hat_p_ci, 0) || g_tol.less(1, f_hat_p_ci) || g_tol.less(1, sum_of_children))
        {
          return false;
        }
        Fhat.set(0, p, c, 1 - sum_of_descendant_states);
      }
    }
    return true;
  }
  else
  {
    const RealTensor& F_lb = _noisyG.F_lb();
    const RealTensor& F_ub = _noisyG.F_ub();
    
    const auto& X_ci = _G.nodeToCharStateList(v_ci);
    for (auto it_ci = X_ci.rbegin(); it_ci != X_ci.rend(); ++it_ci)
    {
      const IntPair& ci = *it_ci;
      for (int p = 0; p < m; ++p)
      {
        double sum_of_children = 0;
        if (it_ci == X_ci.rbegin())
        {
          for (SubOutArcIt a_cidj(T, v_ci); a_cidj != lemon::INVALID; ++a_cidj)
          {
            Node v_dj = T.target(a_cidj);
            const IntPair& dj = *(_G.nodeToCharState(v_dj).begin());
            
            double f_cum_hat_p_dj = Fhat.getCumFreq(p, dj.first, _G.S(dj.first).D(dj.second));
            sum_of_children += f_cum_hat_p_dj;
          }
        }
        else
        {
          const IntPair& pi_ci = *std::prev(it_ci);
          sum_of_children = Fhat.getCumFreq(p, pi_ci.first, _G.S(pi_ci.first).D(pi_ci.second));
        }
        
        double sum_of_descendant_states = 0; // proper descendants
        const IntSet& D_ci = _G.S(ci.first).D(ci.second);
        for (int j = 0; j < k; ++j)
        {
          if (j == ci.second) continue;
          if (D_ci.find(j) == D_ci.end()) continue;
          assert(_G.charStateToNode(ci.first, j) != lemon::INVALID);
          sum_of_descendant_states += Fhat(j, p, ci.first);
        }
        
        double l_p_ci = F_lb(ci.second, p, ci.first);
        double u_p_ci = F_ub(ci.second, p, ci.first);
        double f_hat_p_ci = std::max(l_p_ci, sum_of_children - sum_of_descendant_states);
        Fhat.set(ci.second, p, ci.first, f_hat_p_ci);
        if (g_tol.less(f_hat_p_ci, l_p_ci) || g_tol.less(u_p_ci, f_hat_p_ci))
        {
          return false;
        }
      }
    }
    
    // recurse
    return updateFhat(T, T.source(SubInArcIt(T, v_ci)), Fhat);
  }
}
  
void RootedCladisticNoisyEnumeration::init(SubDigraph& subG,
                                           SubDigraph& T,
                                           ArcList& H,
                                           RealTensor& Fhat)
{
  RootedCladisticEnumeration::init(subG, T, H);
  
  Fhat = _noisyG.F_lb();
  for (int c = 0; c < Fhat.n(); ++c)
  {
    for (int p = 0; p < Fhat.m(); ++p)
    {
      double fhat_p_c0 = 0;
      for (int i = 1; i < Fhat.k(); ++i)
      {
        fhat_p_c0 += Fhat(i, p, c);
      }
      fhat_p_c0 = 1 - fhat_p_c0;
      Fhat.set(0, p, c, fhat_p_c0);
    }
  }
  
#ifdef DEBUG
  RealTensor Fhat2 = _noisyG.F_lb();
  isValid(T, _G.root(), Fhat2);
  assert(Fhat == Fhat2);
#endif
}
  
void RootedCladisticNoisyEnumeration::run()
{
  const Digraph& G = _G.G();
  Node root = _G.root();
  
  _timer.start();
  _counter = 0;
//  if (_threads == 1)
//  {
//    BoolNodeMap filterNodesT(G, false);
//    BoolArcMap filterArcsT(G, false);
//    SubDigraph T(G, filterNodesT, filterArcsT);
//    
//    BoolNodeMap filterNodesG(G, true);
//    BoolArcMap filterArcsG(G, true);
//    SubDigraph subG(G, filterNodesG, filterArcsG);
//    
//    ArcList H;
//    RealTensor Fhat;
//    
//    init(subG, T, H, Fhat);
//    grow(subG, T, H, Fhat);
//  }
//  else
  if (!_monoclonal)
  {
    BoolNodeMap filterNodesT(G, false);
    BoolArcMap filterArcsT(G, false);
    SubDigraph T(G, filterNodesT, filterArcsT);
    filterNodesT[root] = true;
    
    BoolNodeMap filterNodesG(G, true);
    BoolArcMap filterArcsG(G, true);
    SubDigraph subG(G, filterNodesG, filterArcsG);
    
    ArcList H;
    RealTensor Fhat;
    
    init(subG, T, H, Fhat);
    grow(subG, T, H, Fhat);
  }
  if (_monoclonal && _fixTrunk)
  {
    // find monoclonal root
    OutArcIt a_00ci(G, root);
    Node v_ci = G.target(a_00ci);
    
    for (OutArcIt a_cidj(G, v_ci); a_cidj != lemon::INVALID; ++a_cidj)
    {
      _threadGroup.create_thread(boost::bind(&RootedCladisticNoisyEnumeration::runArc, this, a_cidj));
    }
    
    _threadGroup.join_all();
  }
  else if (_monoclonal)
  {
    for (OutArcIt a_00dj(G, root); a_00dj != lemon::INVALID; ++a_00dj)
    {
      _threadGroup.create_thread(boost::bind(&RootedCladisticNoisyEnumeration::runArc, this, a_00dj));
    }
    
    _threadGroup.join_all();
  }
  
  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cerr << "\r" << std::flush;
  }
}
  
void RootedCladisticNoisyEnumeration::runArc(Arc a_00dj)
{
  _sem.wait();
  const Digraph& G = _G.G();
  
  BoolNodeMap filterNodesT(G, false);
  BoolArcMap filterArcsT(G, false);
  SubDigraph T(G, filterNodesT, filterArcsT);
  
  BoolNodeMap filterNodesG(G, true);
  BoolArcMap filterArcsG(G, true);
  SubDigraph subG(G, filterNodesG, filterArcsG);
  
  ArcList H;
  RealTensor Fhat;
  
  init(a_00dj, subG, T, H, Fhat);
  grow(subG, T, H, Fhat);
  
  _sem.post();
}
  
//void RootedCladisticNoisyEnumeration::run()
//{
//  const Digraph& G = _G.G();
//  
//  BoolNodeMap filterNodesT(G, false);
//  BoolArcMap filterArcsT(G, false);
//  SubDigraph T(G, filterNodesT, filterArcsT);
//  
//  BoolNodeMap filterNodesG(G, true);
//  BoolArcMap filterArcsG(G, true);
//  SubDigraph subG(G, filterNodesG, filterArcsG);
//  
//  ArcList F;
//  
//  init(subG, T, F);
//  grow(subG, T, F);
//  
//}
  
void RootedCladisticNoisyEnumeration::writeDOT(std::ostream& out,
                                               const SubDigraph& T,
                                               const RealTensor& F_hat) const
{
  const Digraph& G = _G.G();
  const int m = _G.F().m();
  
  out << "digraph T {" << std::endl;
  out.precision(3);
  
  for (SubNodeIt v_ci(T); v_ci != lemon::INVALID; ++v_ci)
  {
    const auto& X_ci = _G.nodeToCharState(v_ci);
    out << "\t" << _G.G().id(v_ci) << " [label=\"" << _G.label(v_ci) << "\\n";
    
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
      bool first = true;
      for (const IntPair& ci : X_ci)
      {
        if (first)
          first = false;
        else
          out << " | ";
        out << _noisyG.F_lb().getCumFreq(p, ci.first, _G.S(ci.first).D(ci.second))
            << " " << F_hat.getCumFreq(p, ci.first, _G.S(ci.first).D(ci.second))
            << " " << _noisyG.F_ub().getCumFreq(p, ci.first, _G.S(ci.first).D(ci.second));
      }
      out << "\\n";
    }
    out << "\"]" << std::endl;
  }
  
  for (SubArcIt a(T); a != lemon::INVALID; ++a)
  {
    out << "\t" << G.id(G.source(a)) << " -> " << G.id(G.target(a)) << std::endl;
  }
  out << "}" << std::endl;
}
  
bool RootedCladisticNoisyEnumeration::grow(SubDigraph& G,
                                           SubDigraph& T,
                                           ArcList& H,
                                           RealTensor& Fhat)
{
  if (limitReached())
  {
    return true;
  }
  else if (H.empty())
  {
    int m = lemon::countArcs(T);
    if (m >= std::max(_lowerbound, _objectiveValue))
      return finalize(T);
  }
  else if (prune(T, H))
  {
    return false;
  }
  else
  {
    ArcList HH;
    
    do
    {
      assert(!H.empty());
      
      Arc a_cidj = H.back();
      H.pop_back();
      
      const Node v_ci = G.source(a_cidj);
      const Node v_dj = G.target(a_cidj);
      
      assert(T.status(v_ci));
      assert(!T.status(v_dj));
      assert(!T.status(a_cidj));
      
      // add a_cidj to T
      assert(isValid(T, a_cidj));
      addArc(T, Fhat, a_cidj);
      
      ArcList newH = H;
      
      // remove each arc wv where w in T from F
      for (ArcListNonConstIt it = newH.begin(); it != newH.end();)
      {
        Arc a_elfs = *it;
        Node v_el = G.source(a_elfs);
        Node v_fs = G.target(a_elfs);
        if (T.status(v_fs))
        {
          it = newH.erase(it);
        }
        else if (v_fs == v_dj)
        {
          assert(T.status(v_el));
          it = newH.erase(it);
        }
        else if (!checkFhat(T, Fhat, a_elfs))
        {
          it = newH.erase(it);
        }
        else
        {
//          if (!isValid(T, a_elfs))
//          {
//            RootedCladisticEnumeration::writeDOT(std::cout, T, HH);
//          }
          assert(checkFhat(T, Fhat, *it));
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
          assert(0 <= pi_l && pi_l < Fhat.k());
          Node v_e_pi_l = _G.charStateToNode(el.first, pi_l);

          isConsistent = isConsistent && isFirstAncestor(T, el.first, v_e_pi_l, v_dj);
        }
        
        if (!isConsistent)
          continue;
        
        if (checkFhat(T, Fhat, a_djel))
        {
          assert(isValid(T, a_djel));
          newH.push_back(a_djel);
        }
      }
      
      // shuffle frontier, depth-first is overrated
      if (_monoclonal)
      {
        std::vector<Arc> tmpVec;
        tmpVec.insert(tmpVec.begin(), newH.begin(), newH.end());
        std::shuffle(tmpVec.begin(), tmpVec.end(), g_rng);
        newH.clear();
        for (Arc a : tmpVec)
        {
          newH.push_back(a);
        }
      }
      
//      RootedCladisticEnumeration::writeDOT(std::cout, T, newH);
      
      if (grow(G, T, newH, Fhat))
        return true;
      
      G.disable(a_cidj);
      
      removeArc(T, Fhat, a_cidj);
      
      HH.push_back(a_cidj);
    } while (!H.empty());
    
    for (ArcListRevIt it = HH.rbegin(); it != HH.rend(); ++it)
    {
      Arc a = *it;
      assert(!G.status(a));
      
      H.push_back(*it);
      G.enable(a);
    }
    
    if (limitReached())
      return true;
    else
      return false;
  }
  
  if (limitReached())
    return true;
  else
    return false;
}
  
bool RootedCladisticNoisyEnumeration::isValid(const SubDigraph& T) const
{
  RealTensor F_hat = _noisyG.F_lb();
  bool res = isValid(T, _G.root(), F_hat);
  
  return res;
}
  
bool RootedCladisticNoisyEnumeration::isValid(const SubDigraph& T,
                                              Arc a_cidj) const
{
  assert(!T.status(a_cidj));
  assert(!T.status(T.target(a_cidj)));
  
  T.enable(a_cidj);
  T.enable(T.target(a_cidj));
  
  bool res = isValid(T);
  T.disable(a_cidj);
  T.disable(T.target(a_cidj));
  return res;
}
  
bool RootedCladisticNoisyEnumeration::isValid(const SubDigraph& T,
                                              Node v_ci,
                                              RealTensor& F_hat) const
{
  assert(T.status(v_ci));
  
  const RealTensor& F_lb = _noisyG.F_lb();
  const RealTensor& F_ub = _noisyG.F_ub();
  const int m = F_lb.m();
  const int n = F_lb.n();
  const int k = F_lb.k();
  
  for (SubOutArcIt a_cidj(T, v_ci); a_cidj != lemon::INVALID; ++a_cidj)
  {
    Node v_dj = T.target(a_cidj);
    if (!isValid(T, v_dj, F_hat))
    {
      return false;
    }
  }
  
  if (v_ci == _G.root())
  {
    // set \hat{f}_{p,(c,i)}
    for (int p = 0; p < m; ++p)
    {
      double sum_of_children = 0;
      for (SubOutArcIt a_cidj(T, v_ci); a_cidj != lemon::INVALID; ++a_cidj)
      {
        Node v_dj = T.target(a_cidj);
        
        const IntPair& dj = *(_G.nodeToCharState(v_dj).begin());
        
        double f_cum_hat_p_dj = F_hat.getCumFreq(p, dj.first, _G.S(dj.first).D(dj.second));
        sum_of_children += f_cum_hat_p_dj;
      }
      
      for (int c = 0; c < n; ++c)
      {
        double sum_of_descendant_states = 0;
        for (int j = 1; j < k; ++j)
        {
          if (_G.charStateToNode(c, j) != lemon::INVALID)
          {
            sum_of_descendant_states += F_hat(j, p, c);
          }
        }

        double f_hat_p_ci = 1 - sum_of_descendant_states;
        if (g_tol.less(f_hat_p_ci, 0) || g_tol.less(1, f_hat_p_ci) || g_tol.less(1, sum_of_children))
        {
          return false;
        }
        else
        {
          // make sure that frequencies sum to 1. also, we may violate u_pci...
          F_hat.set(0, p, c, 1 - sum_of_descendant_states);
        }
      }
    }
    return true;
  }
  else
  {
    const auto& X_ci = _G.nodeToCharStateList(v_ci);
    for (auto it_ci = X_ci.rbegin(); it_ci != X_ci.rend(); ++it_ci)
    {
      const IntPair& ci = *it_ci;
      
      // set \hat{f}_{p,(c,i)}
      for (int p = 0; p < m; ++p)
      {
        double sum_of_children = 0;
        if (it_ci == X_ci.rbegin())
        {
          for (SubOutArcIt a_cidj(T, v_ci); a_cidj != lemon::INVALID; ++a_cidj)
          {
            Node v_dj = T.target(a_cidj);
            const IntPair& dj = *(_G.nodeToCharState(v_dj).begin());
            
            double f_cum_hat_p_dj = F_hat.getCumFreq(p, dj.first, _G.S(dj.first).D(dj.second));
            sum_of_children += f_cum_hat_p_dj;
          }
        }
        else
        {
          const IntPair pi_ci = *std::prev(it_ci);
          sum_of_children = F_hat.getCumFreq(p, pi_ci.first, _G.S(pi_ci.first).D(pi_ci.second));
        }
      
        double sum_of_descendant_states = 0; // proper descendants
        const IntSet& D_ci = _G.S(ci.first).D(ci.second);
        for (int j = 0; j < k; ++j)
        {
          if (j == ci.second) continue;
          if (D_ci.find(j) == D_ci.end()) continue;
          if (_G.charStateToNode(ci.first, j) != lemon::INVALID)
          {
            sum_of_descendant_states += F_hat(j, p, ci.first);
          }
        }
        
        double l_p_ci = F_lb(ci.second, p, ci.first);
        double u_p_ci = F_ub(ci.second, p, ci.first);
        assert(l_p_ci <= u_p_ci);
        double f_hat_p_ci = std::max(l_p_ci, sum_of_children - sum_of_descendant_states);
        F_hat.set(ci.second, p, ci.first, f_hat_p_ci);
        if (g_tol.less(f_hat_p_ci, l_p_ci) || g_tol.less(u_p_ci, f_hat_p_ci))
        {
          return false;
        }
      }
    }
    
    return true;
  }
}
  
void RootedCladisticNoisyEnumeration::addArc(SubDigraph& T,
                                             RealTensor& Fhat,
                                             Arc a_cidj) const
{
  const Node v_ci = T.source(a_cidj);
  const Node v_dj = T.target(a_cidj);
  
//  RealTensor FF = _noisyG.F_lb();
//  isValid(T, _G.root(), FF);
//  if (Fhat != FF)
//  {
//    writeDOT(std::cerr, T, Fhat);
//    writeDOT(std::cerr, T, FF);
//  }
  
  // add a_cidj to T
  T.enable(v_dj);
  T.enable(a_cidj);

#ifdef DEBUG
  RealTensor Fhat_old = Fhat;
#endif
  
  bool res = updateFhat(T, v_dj, Fhat);
  
#ifdef DEBUG
  RealTensor Fnew = _noisyG.F_lb();
  isValid(T, _G.root(), Fnew);
  if (Fnew != Fhat)
  {
    writeDOT(std::cerr, T, Fhat);
    writeDOT(std::cerr, T, Fnew);
    updateFhat(T, v_ci, Fhat_old);
  }
  assert(Fnew == Fhat);
#endif
  
//  writeDOT(std::cerr, T, Fhat);
//  if (!res)
//  {
//    assert(isValid(T));
//    writeDOT(std::cerr, T, Fhat_old);
//    updateFhat(T, v_ci, Fhat_old);
//    writeDOT(std::cerr, T, Fhat);
//    
//  }
  assert(res);
  
  assert(isArborescence(T));
  assert(isValid(T));
}
  
void RootedCladisticNoisyEnumeration::removeArc(SubDigraph& T,
                                                RealTensor& Fhat,
                                                Arc a_cidj) const
{
  assert(T.status(a_cidj));
  
  const Node v_ci = T.source(a_cidj);
  const Node v_dj = T.target(a_cidj);
  
//  RealTensor FF = _noisyG.F_lb();
//  isValid(T, _G.root(), FF);
//  if (Fhat != FF)
//  {
//    writeDOT(std::cerr, T, Fhat);
//    writeDOT(std::cerr, T, FF);
//  }
  
  // remove a_cidj from T
  T.disable(a_cidj);
  T.disable(v_dj);
  
//  RealTensor Fhat_old = Fhat;
  
  bool res = updateFhat(T, v_ci, Fhat);
  
//  RealTensor Fnew = _noisyG.F_lb();
//  isValid(T, _G.root(), Fnew);
//  if (Fnew != Fhat)
//  {
//    writeDOT(std::cerr, T, Fhat);
//    writeDOT(std::cerr, T, Fnew);
//    updateFhat(T, v_ci, Fhat_old);
//  }
  
//  bool res = updateFhat(T, v_ci, Fhat);
  assert(res);
  
  assert(isArborescence(T));
  assert(isValid(T));
}

  
} // namespace gm
