/*
 * rootedcladisticenumeration.h
 *
 *  Created on: 29-sep-2015
 *      Author: M. El-Kebir
 */

#ifndef ROOTEDCLADISTICENUMERATION_H
#define ROOTEDCLADISTICENUMERATION_H

#include <lemon/adaptors.h>
#include <lemon/bfs.h>
#include <lemon/time_measure.h>
#include <boost/asio/signal_set.hpp>
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/interprocess/sync/interprocess_semaphore.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>
#include "utils.h"
#include "rootedcladisticancestrygraph.h"
#include "solution.h"
#include "solutionset.h"

namespace gm {

class RootedCladisticEnumeration
{
public:
  DIGRAPH_TYPEDEFS(Digraph);
  
  RootedCladisticEnumeration(const RootedCladisticAncestryGraph& G,
                             int limit,
                             int timeLimit,
                             int threads,
                             int lowerbound,
                             bool monoclonal,
                             bool fixTrunk,
                             const IntSet& whitelist);
  
  virtual ~RootedCladisticEnumeration();
  
  virtual void run();
  
  int solutionCount() const
  {
    return _result.size();
  }
  
  void stop()
  {
    _mutex.lock();
    _threadGroup.interrupt_all();
  }
  
  Solution solution(int solIdx) const;
  
  void populateSolutionSet(SolutionSet& sols) const;
  
  int objectiveValue() const
  {
    return _objectiveValue;
  }
  
  std::string newick(int solIdx) const;

protected:
  typedef std::list<Arc> ArcList;
  typedef ArcList::const_iterator ArcListIt;
  typedef ArcList::iterator ArcListNonConstIt;
  typedef ArcList::const_reverse_iterator ArcListRevIt;
  typedef std::vector<ArcList> ArcListVector;
  typedef std::list<ArcList> ArcListList;
  typedef lemon::SubDigraph<const Digraph> SubDigraph;
  typedef SubDigraph::ArcIt SubArcIt;
  typedef SubDigraph::NodeIt SubNodeIt;
  typedef SubDigraph::OutArcIt SubOutArcIt;
  typedef SubDigraph::InArcIt SubInArcIt;
  typedef Digraph::NodeMap<IntSet> IntSetNodeMap;
  
  virtual void runArc(Arc a_00dj);
  
  void init(SubDigraph& subG, SubDigraph& T, ArcList& F);
  void init(Arc a_00dj, SubDigraph& subG, SubDigraph& T, ArcList& F);
  
  bool prune(const SubDigraph& T,
             const ArcList& F) const;

  bool isFirstAncestor(const SubDigraph& T,
                       int c,
                       Node v_ci,
                       Node v_dj) const;
  
  bool isAncestor(const SubDigraph& T,
                  Node v_dj,
                  Node v_ci) const;
  
  bool isStateComplete(const SubDigraph& T,
                       int c) const
  {
    const int n = _G.F().n();
    const int k = _G.F().k();
    
    assert(0 <= c && c < n);
    for (int i = 1; i < k; ++i)
    {
//      _G.S(c).writeDOT(std::cout);
//      _G.writeDOT(std::cout);
      if (_G.S(c).isPresent(i) &&
          !T.status(_G.charStateToNode(c, i)))
      {
        return false;
      }
    }
    
    return true;
  }
  
  void writeDOT(std::ostream& out,
                const SubDigraph& T,
                const ArcList& H) const;
  
  void newick(const SubDigraph& T,
              Node v_ci,
              std::string& str) const;

  bool isStateComplete(const SubDigraph& T) const
  {
    const int n = _G.F().n();
    
    for (int c = 0; c < n; ++c)
    {
      if (!isStateComplete(T, c))
        return false;
    }
    
    return true;
  }
  
  int makeStateComplete(SubDigraph& T) const
  {
    typedef std::set<Node> NodeSet;
    typedef NodeSet::const_iterator NodeSetIt;
    
    const int n = _G.F().n();
    BoolVector stateComplete(n, false);
    for (int c = 0; c < n; ++c)
    {
      stateComplete[c] = isStateComplete(T, c);
    }
    
    NodeSet nodesToRemove;
    for (SubNodeIt v_ci(T); v_ci != lemon::INVALID; ++v_ci)
    {
      // don't remove the root node!
      if (v_ci == _G.root())
        continue;
      
      // is there a character-state pair in v_ci whose character c is state complete?
      bool somethingStateComplete_v_ci = false;
      for (const IntPair& ci : _G.nodeToCharState(v_ci))
      {
        somethingStateComplete_v_ci = somethingStateComplete_v_ci || stateComplete[ci.first];
      }

      if (!somethingStateComplete_v_ci)
        nodesToRemove.insert(v_ci);
    }
    
    for (NodeSetIt it = nodesToRemove.begin();
         it != nodesToRemove.end(); ++it)
    {
      T.disable(*it);
    }
    
    // T is state complete now, final step is to count # char-state pairs in T
    BoolNodeMap reached(_G.G(), false);
    lemon::bfs(T).reachedMap(reached).run(_G.root());
    
    bool nodeRemoved = false;
    for (NodeIt v_ci(_G.G()); v_ci != lemon::INVALID; ++v_ci)
    {
//      if (reached[v_ci])
//      {
//        for (const IntPair& ci : _G.nodeToCharState(v_ci))
//        {
//          if (stateComplete[ci.first])
//          {
//            ++res;
//          }
//        }
//      }
//      else
      if (!reached[v_ci] && T.status(v_ci))
      {
        T.disable(v_ci);
        nodeRemoved = true;
      }
    }
    
    int res = 0;
    for (int c = 0; c < n; ++c)
    {
      if (isStateComplete(T, c))
        ++res;
    }
    
    // if we removed disconnected vertices, then we need to recheck again, as we may have broken state-completeness for a different character
    if (nodeRemoved)
    {
      return makeStateComplete(T);
    }
    else
    {
      return res;
    }
  }
  
  bool isArborescence(const SubDigraph& T) const;
  
  void addArc(SubDigraph& T,
              Arc a_cidj) const;
  void removeArc(SubDigraph& T,
                 Arc a_cidj) const;
  
  bool finalize(SubDigraph& T);
  
  bool limitReached() const
  {
    boost::interprocess::scoped_lock<boost::mutex> lock(_mutex);
    return (_limit != -1 && _counter >= _limit) || (_timeLimit != -1 && _timer.realTime() > _timeLimit);
  }
  
  struct Compare
  {
  public:
    Compare(const RootedCladisticAncestryGraph& G)
      : _G(G)
    {
    }
    
    bool operator()(Arc a1, Arc a2)
    {
      Node v1 = _G.G().target(a1);
      const IntPair& ci_1 = *(_G.nodeToCharState(v1).begin());
      const IntSet& D1 = _G.S()[ci_1.first].D(ci_1.second);
      
      Node v2 = _G.G().target(a2);
      const IntPair& ci_2 = *(_G.nodeToCharState(v2).begin());
      const IntSet& D2 = _G.S()[ci_2.first].D(ci_2.second);
      
      const int m = _G.F().m();
      double max1 = 0, max2 = 0;
      for (int p = 0; p < m; ++p)
      {
        max1 = std::max(max1, _G.F().getCumFreq(p, ci_1.first, D1));
        max2 = std::max(max2, _G.F().getCumFreq(p, ci_2.first, D2));
      }
      
      return max1 > max2;
    }
    
  private:
    const RootedCladisticAncestryGraph& _G;
  };
  
  struct Compare2
  {
  public:
    Compare2(const RootedCladisticAncestryGraph& G)
      : _G(G)
    {
    }
    
    bool operator()(Arc a1, Arc a2)
    {
      Node v1 = _G.G().target(a1);
      const IntPair& ci_1 = *(_G.nodeToCharState(v1).begin());
      
      Node v2 = _G.G().target(a2);
      const IntPair& ci_2 = *(_G.nodeToCharState(v2).begin());
      
      return ci_1.second > ci_2.second;
    }
    
  private:
    const RootedCladisticAncestryGraph& _G;
  };
  
private:
  bool grow(SubDigraph& G,
            SubDigraph& T,
            ArcList& F);
  
  virtual bool isValid(const SubDigraph& T) const;
  virtual bool isValid(const SubDigraph& T, Arc a_ciel) const;
  
  void initA(const SubDigraph& T,
             Node v_ci,
             const BoolVector& stateComplete,
             PerfectPhyloMatrix& A) const;
  void initU(const SubDigraph& T,
             const RealTensor& F,
             Node v_ci,
             const BoolVector& stateComplete,
             RealMatrix& U) const;
  
  virtual void initF(int solIdx, RealTensor& F) const
  {
    F = _G.F();
  }
  
protected:
  const RootedCladisticAncestryGraph& _G;
  
  ArcListList _result;
  int _objectiveValue;
  
  int _limit;
  int _timeLimit;
  int _threads;
  int _lowerbound;
  int _counter;
  
  mutable boost::mutex _mutex;
  boost::interprocess::interprocess_semaphore _sem;
  boost::thread_group _threadGroup;
  
  lemon::Timer _timer;
  bool _monoclonal;
  bool _fixTrunk;
  
  const IntSet& _whiteList;
};
  
} // namespace gm

#endif // ROOTEDCLADISTICENUMERATION_H
