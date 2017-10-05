/*
 * cell.cpp
 *
 *  Created on: 25-aug-2017
 *      Author: M. El-Kebir
 */

#include "cell.h"

Cell::Cell(const IntVector& passengerMutations,
           int mutation,
           int anatomicalSite)
  : _passengerMutations(passengerMutations)
  , _mutation(mutation)
  , _anatomicalSite(anatomicalSite)
{
}

Cell::Cell()
  : _passengerMutations()
  , _mutation(-1)
  , _anatomicalSite(-1)
{
}

Cell::Outcome Cell::performGeneration(double logisticFactor,
                                      int nrDriverMutations) const
{
  /*
   * one generation: possible replication (+mutation) and possible death
   * cell division are binary: cell i divides into i and j
   * j.parent = i
   * only cell j acquires a single new mutation
   * this function returns the daughter cells that are present in the next generation
   */
  
  std::uniform_real_distribution<> unif(0, 1);

  const double s_0 = 0.1 * logisticFactor;
  
  double birthRate = 0.5 * (1 + s_0);
  for (int i = 0; i < nrDriverMutations; ++i)
  {
    birthRate *= (1 + s_0);
  }
  
  if (birthRate < 0)
    birthRate = 0;
  else if (birthRate > 1)
    birthRate = 1;
  
  double r = unif(g_rng);
  if (r < birthRate)
  {
    return REPLICATION;
  }
  else
  {
    return DEATH;
  }
}
