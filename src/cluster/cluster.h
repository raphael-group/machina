/*
 * cluster.h
 *
 *  Created on: 7-sep-2017
 *      Author: M. El-Kebir
 */

#ifndef CLUSTER_H
#define CLUSTER_H

#include "utils.h"
#include "readmatrix.h"
#include "frequencymatrix.h"

/// This class clusters SNVs given a matrix of read counts
class Cluster
{
public:
  /// Constructor
  ///
  /// @param R Read matrix
  /// @param alpha Confidence interval used for clustering
  /// @param relabel Relabel mutation clusters
  Cluster(const ReadMatrix& R,
          double alpha,
          bool relabel);
  
  /// Clonality status of a mutation
  enum MutationType
  {
    ABSENT,
    SUBCLONAL,
    CLONAL,
  };
  
  typedef std::vector<MutationType> ProfileVector;
  
  /// Return clonality status profile of the specified mutation
  ///
  /// @param i Mutation
  ProfileVector getProfile(int i) const
  {
    assert(0 <= i && i < _R.getNrCharacters());
    const int k = _F.getNrSamples();
    
    ProfileVector profile_i(k, ABSENT);
    for (int p = 0; p < k; ++p)
    {
      profile_i[p] = _profile[p][i];
    }
    
    return profile_i;
  }
  
  /// Write clustering
  ///
  /// @param out Output stream
  void writeClustering(std::ostream& out) const;
  
  /// Write AncesTree input file
  ///
  /// @param out Output stream
  void writeAncesTreeInput(std::ostream& out) const;
  
  /// Return pooled read matrix
  const ReadMatrix& getClusteredR() const
  {
    return _newR;
  }
  
  /// Return clustered frequency matrix
  const FrequencyMatrix& getClusteredF() const
  {
    return _newF;
  }
  
  /// Cluster mutations based on clonality status
  ///
  /// @param beta Confidence interval used for pooled frequency matrix
  void clusterClonalityStatus(double beta);
  
  /// Cluster mutations by finding connected components
  ///
  /// @param beta Confidence interval used for pooled frequency matrix
  void clusterCC(double beta);
  
  /// Return clustering
  const IntMatrix& getClustering() const
  {
    return _clustering;
  }
  
private:  
  void classify();
  
  typedef std::vector<ProfileVector> ProfileMatrix;
  
private:
  /// Matrix of alternate and reference read counts for each SNV in each sample
  const ReadMatrix& _R;
  /// Matrix of cell frequency confidence intervals for each SNV in each sample
  FrequencyMatrix _F;
  /// Matrix of mutation status of each SNV in each sample
  ProfileMatrix _profile;
  /// Cluster assignment of each SNV
  IntMatrix _clustering;
  /// Flag indicating whether clusters should receive new labels
  bool _relabel;
  /// Results matrix of alternate and reference read counts for each mutation cluster in each sample
  ReadMatrix _newR;
  /// Results matrix of mutation cluster frequencies in each sample
  FrequencyMatrix _newF;
};

#endif // CLUSTER_H
