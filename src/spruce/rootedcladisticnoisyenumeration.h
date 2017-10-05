/*
 * rootedcladisticnoisyenumeration.h
 *
 *  Created on: 11-oct-2015
 *      Author: M. El-Kebir
 */

#ifndef ROOTEDCLADISTICNOISYENUMERATION_H
#define ROOTEDCLADISTICNOISYENUMERATION_H

#include "rootedcladisticenumeration.h"
#include "rootedcladisticnoisyancestrygraph.h"

namespace gm {
  
class RootedCladisticNoisyEnumeration : public RootedCladisticEnumeration
{
public:
  RootedCladisticNoisyEnumeration(const RootedCladisticNoisyAncestryGraph& G,
                                  int limit,
                                  int timeLimit,
                                  int threads,
                                  int lowerbound,
                                  bool monoclonal,
                                  bool fixTrunk,
                                  const IntSet& whiteList);
  
  const RootedCladisticNoisyAncestryGraph& noisyG() const
  {
    return _noisyG;
  }
  
  virtual void run();
  
protected:
  void addArc(SubDigraph& T,
              RealTensor& Fhat,
              Arc a_cidj) const;
  void removeArc(SubDigraph& T,
                 RealTensor& Fhat,
                 Arc a_cidj) const;
  
private:
  bool grow(SubDigraph& G,
            SubDigraph& T,
            ArcList& H,
            RealTensor& Fhat);
   
  void writeDOT(std::ostream& out,
                const SubDigraph& T,
                const RealTensor& F_hat) const;
  
  void init(SubDigraph& subG,
            SubDigraph& T,
            ArcList& H,
            RealTensor& Fhat);
  void init(Arc a_00dj,
            SubDigraph& subG,
            SubDigraph& T,
            ArcList& H,
            RealTensor& Fhat);
  
  bool updateFhat(const SubDigraph& T,
                  Node v_ci,
                  RealTensor& Fhat) const;
  
  bool checkFhat(SubDigraph& T,
                 RealTensor& Fhat,
                 Arc a_cidj) const;
  
  void initA(const SubDigraph& T, Node v_ci, PerfectPhyloMatrix& A) const;
  void initU(const SubDigraph& T, Node v_ci, RealMatrix& U) const;
    
  bool isValid(const SubDigraph& T) const;
  
  bool isValid(const SubDigraph& T,
               Arc a_cidj) const;
  
  virtual void runArc(Arc a_00dj);
  
  virtual void initF(int solIdx, RealTensor& F) const
  {
    const Digraph& G = _G.G();
    
    BoolNodeMap filterNodesT(G, false);
    BoolArcMap filterArcsT(G, false);
    
    SubDigraph T(G, filterNodesT, filterArcsT);
    
    ArcListList::const_iterator it = _result.begin();
    int tmp = solIdx;
    for (; tmp > 0; --tmp, ++it);
    const ArcList& arcs = *it;
    for (ArcListIt it = arcs.begin(); it != arcs.end(); ++it)
    {
      Arc a_cidj = *it;
      Node v_ci = G.source(a_cidj);
      Node v_dj = G.target(a_cidj);
      T.enable(v_ci);
      T.enable(v_dj);
      T.enable(a_cidj);
    }
    
    isValid(T, _G.root(), F);

    for (int p = 0; p < F.m(); ++p)
    {
      F.setRowLabel(p, _G.F().getRowLabel(p));
    }
    for (int c = 0; c < F.n(); ++c)
    {
      F.setColLabel(c, _G.F().getColLabel(c));
    }
  }
  
protected:
  const RootedCladisticNoisyAncestryGraph& _noisyG;
  
  bool isValid(const SubDigraph& T,
               Node v_ci,
               RealTensor& F_hat) const;
};
  
} // namespace gm

#endif // ROOTEDCLADISTICNOISYENUMERATION_H
