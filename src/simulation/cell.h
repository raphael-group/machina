/*
 * cell.h
 *
 *  Created on: 25-aug-2017
 *      Author: M. El-Kebir
 */

#ifndef CELL_H
#define CELL_H

#include "utils.h"

/// This class models a cell, the central unit in the agent-based model
class Cell
{
public:
  /// Constructor
  Cell(const IntVector& passengerMutations,
       int mutation,
       int anatomicalSite);
  
  /// Default constructor
  Cell();
  
  /// Cell cycle outcome
  enum Outcome
  {
    /// Replication
    REPLICATION,
    /// Death
    DEATH
  };
  
  /// Perform a cell cycle generation
  ///
  /// @param _anatomicalSiteFactors Record of the clonal
  /// expansions in this anatomical site
  Outcome performGeneration(double logisticFactor,
                            int nrDriverMutations) const;
  
  /// Migrate to specified anatomical site
  ///
  /// @param newAnatomicalSite
  void migrate(int newAnatomicalSite)
  {
    _anatomicalSite = newAnatomicalSite;
  }
  
  /// Return anatomical site
  int getAnatomicalSite() const
  {
    return _anatomicalSite;
  }
  
  /// Return passenger mutations
  const IntVector& getPassengerMutations() const
  {
    return _passengerMutations;
  }
  
  /// Return the mutation introduced by this cell
  int getMutation() const
  {
    return _mutation;
  }
  
protected:
  /// Set of passenger mutations
  IntVector _passengerMutations;
  /// Introduced mutation
  int _mutation;
  /// Anatomical site
  int _anatomicalSite;
};

#endif // CELL_H
