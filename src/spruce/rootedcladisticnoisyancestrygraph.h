/*
 *  rootedcladisticnoisyancestrygraph.h
 *
 *   Created on: 11-oct-2015
 *       Author: M. El-Kebir
 */

#ifndef ROOTEDCLADISTICNOISYANCESTRYGRAPH_H
#define ROOTEDCLADISTICNOISYANCESTRYGRAPH_H

#include "utils.h"
#include "rootedcladisticancestrygraph.h"

namespace gm {
  
class RootedCladisticNoisyAncestryGraph : public RootedCladisticAncestryGraph
{
public:
  RootedCladisticNoisyAncestryGraph(const StateTreeVector& S,
                                    const RealTensor& F_lb,
                                    const RealTensor& F_ub);
  
  virtual void init();
  
  void initMultiState(const IntVector& sampleToAnatomicalSite,
                      const IntSetVector& anatomicalSiteToSamples);
  
  void writeDOT(std::ostream& out) const;
  
  const RealTensor& F_lb() const
  {
    return _F_lb;
  }
  
  const RealTensor& F_ub() const
  {
    return _F_ub;
  }
  
private:
  const RealTensor& _F_lb;
  const RealTensor& _F_ub;
};
  
} // namespace gm

#endif // ROOTEDCLADISTICNOISYANCESTRYGRAPH_H
