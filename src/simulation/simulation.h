/*
 * simulation.h
 *
 *  Created on: 28-aug-2017
 *      Author: M. El-Kebir
 */

#ifndef SIMULATION_H
#define SIMULATION_H

#include "simulation/cell.h"
#include "migrationgraph.h"
#include "utils.h"
#include <memory>
#include <unordered_set>

/// This class implements and extend an agent-based model of tumor progression
/// to include cell migration and the infinite sites assumption
/// This model is described in:
/// [Reiter J.G., Bozic I., Chatterjee K., Nowak M.A. (2013) LNCS].
class Simulation
{
public:
  /// Migration pattern
  enum Pattern
  {
    /// Monoclonal single-source seeding (mS)
    PATTERN_mS,
    /// Single-source seeding (S)
    PATTERN_S,
    /// Mutli-source seeding (M)
    PATTERN_M,
    /// Reseeding (R)
    PATTERN_R
  };
  
  /// Constructor
  ///
  /// @param K Carrying capacity
  /// @param migrationRate Migration rate
  /// @param driverProb Driver mutation probability
  /// @param mutFreqThreshold Mutation frequency threshold
  /// @param maxNrAnatomicalSites Maximum number of detectable anatomical sites
  /// @param nrSamplesPerAnatomicalSite Number of samples per anatomical site
  /// @param nrSamplesPrimary Number of samples of primary tumor
  /// @param targetCoverage Target coverage
  /// @param pattern Migration pattern
  /// @param seqErrorRate Per base sequencing error
  /// @param purity Sample purity
  Simulation(double K,
             double migrationRate,
             double driverProb,
             double mutFreqThreshold,
             int maxNrAnatomicalSites,
             int nrSamplesPerAnatomicalSite,
             int nrSamplesPrimary,
             int targetCoverage,
             Pattern pattern,
             double seqErrorRate,
             double purity);
  
  /// Destructor
  ~Simulation();
  
  /// Simulate a cell tree. Return true if simulation was successful and false otherwise
  ///
  /// @param verbose Output simulation information to stdout
  bool simulate(bool verbose);
  
  /// Return simulated migration graph
  MigrationGraph getMigrationGraph() const;
  
  /// Return observed migration graph
  MigrationGraph getObservedMigrationGraph() const;
  
  /// Return simulated clone tree
  const CloneTree& getCloneTree() const
  {
    assert(_pCloneT);
    return *_pCloneT;
  }
  
  /// Return observed clone tree
  const CloneTree& getObservedCloneTree() const
  {
    assert(_pCloneT);
    return *_pCloneTT;
  }
  
  /// Return observed sample proportions
  const DoubleVectorNodeMap& getObservedSampleProportions() const
  {
    assert(_pSampleProportionsTT);
    return *_pSampleProportionsTT;
  }
  
  /// Return vertex labeling
  const StringNodeMap& getObservedVertexLabeling() const
  {
    return *_pVertexLabelingTT;
  }
  
  /// Return mixture proportions of leaves in simulated clone tree
  const DoubleNodeMap& getAnatomicalSiteProportions() const
  {
    assert(_pAnatomicalSiteProportions);
    return *_pAnatomicalSiteProportions;
  }
  
  /// Return number of mutations in simulated clone tree
  int getNrMutations() const
  {
    assert(_freq.size());
    assert(_freq.front().size());
    return _freq.front().front().size();
  }
  
  /// Write read counts
  ///
  /// @param out Output stream
  void writeReadCounts(std::ostream& out) const;
  
  /// Write observed mutation clusters
  ///
  /// @param out Output stream
  void writeObservedClustering(std::ostream& out) const;
  
  /// Write driver mutations
  ///
  /// @param out Output stream
  void writeDrivers(std::ostream& out) const;
  
private:
  /// Map from set of driver mutations to carrying capacity
  typedef std::map<IntSet, double> AnatomicalSiteMap;
  /// Anatomical site and drivers-specific carrying capacities
  typedef std::vector<AnatomicalSiteMap> AnatomicalSiteFactorVector;
  /// Vector of cells
  typedef std::vector<Cell> CellVector;
  /// Matrix of cells
  typedef std::vector<CellVector> CellMatrix;
  /// Map of set of driver mutations to cell vector (with those drivers)
  typedef std::vector<std::map<IntSet, CellVector> > ClonalComposition;
  
  /// Return a new mutation
  int getNewMutation()
  {
    return _nrMutations++;
  }
  
  /// Clone tree
  void constructCloneTree();
  
  /// Simulate read counts
  bool simulateReadCounts();
  
  /// Sampled clone tree
  bool constructSampledCloneTree();
  
  /// Return new growth rate adjustment
  double getNewGrowthRateDelta() const
  {
    std::uniform_real_distribution<> unif(0, 1);
    
    if (unif(g_rng) < 0.01)
    {
      std::normal_distribution<> normal(0, 0.2);
      return normal(g_rng);
    }
    else
    {
      std::normal_distribution<> normal(0, 0.01);
      return normal(g_rng);
    }
  }

  /// Initialize simulation
  void init();
  
  /// Update anatomical site factors
  void updateAnatomicalSiteFactors();
  
  /// Perform potential migrations
  void migrate();
  
  /// Return target anatomical site
  ///
  /// @param s Source anatomical site
  int getTargetAnatomicalSite(int s);
  
  /// Compute mutation frequencies
  ///
  /// @param freq Output frequencies
  /// @param s Anatomical site
  void getMutationFrequencies(DoubleVector& freq, int s) const;
  
  /// Sort and collapse repeated columns of a binary matrix
  ///
  /// @param matrix Output matrix
  /// @param originalColumns Mapping of new columns back to original columns
  static void sortAndCondense(BoolMatrix& matrix,
                              IntSetVector& originalColumns);
  
  /// Construct perfect phylogeny tree
  ///
  /// @param matrix Input matrix (sorted and condensed)
  /// @param T Output tree
  /// @param root Output root node
  /// @param leafToTaxa Assigns to each leaf a row
  /// @param arcToCharacter Assigns to each arc a column
  /// @param nodeToMutations Assigns to each node a set of mutations
  static void constructTree(const BoolMatrix& matrix,
                            Digraph& T,
                            Node& root,
                            IntNodeMap& leafToTaxa,
                            IntArcMap& arcToCharacter,
                            IntSetNodeMap& nodeToMutations);
  
  /// Construct maps
  static void constructMaps(const Digraph& newT,
                            const Node& newRoot,
                            const IntNodeMap& leafToTaxa,
                            const IntVector& anatomicalSiteMap,
                            const StringVector& anatomicalSiteLabel,
                            const IntArcMap& arcToCharacter,
                            const IntSetNodeMap& nodeToMutations,
                            const IntSetVector& toOriginalColumns,
                            IntNodeMap& anatomicalSite,
                            StringNodeMap& anatomicalSiteStr,
                            StringNodeMap& label,
                            StringToNodeMap& invMap,
                            IntSetNodeMap& stateVector);
  
  /// Internal logging
  void writeDOT(const Digraph& T,
                Node root,
                const IntSetArcMap& mutations,
                const DoubleNodeMap& cloneFrequency,
                const IntNodeMap& anatomicalSite,
                std::ostream& out) const;
  
  typedef std::map<int, int> IntIntMap;
  
  /// Associates to each arc of migration graph the corresponding cell
  /// that migrated and its driver mutations
  typedef Digraph::ArcMap<std::pair<Cell, IntSet> > CellArcMap;
  
  /// Return the number of cells in the provided anatomical site
  /// weighted by the number of driver mutations they possess
  ///
  /// @param s Anatomical site
  double getNrWeightedExtantCells(int s) const
  {
    double res = 0;
    for (const auto& kv : _extantCellsByDrivers[s])
    {
      res += kv.first.size() * kv.second.size();
    }
    return res;
  }
  
  /// Return the number of cells in the provided anatomical site
  ///
  /// @param s Anatomical site
  int getNrExtantCells(int s) const
  {
    int res = 0;
    for (const auto& kv : _extantCellsByDrivers[s])
    {
      res += kv.second.size();
    }
    return res;
  }
  
  /// Determines whether there are migration edges in _pCloneTT
  /// between the same anatomical sites that are incident to the same source vertex
  bool fakeParallelEdges() const;
  
  /// Return the maximum number of drivers possessed by any cell in the specified anatomical site
  ///
  /// @param s Anatomical site
  int getNrNestedDrivers(int s) const
  {
    int res = 0;
    
    for (const auto& kv : _extantCellsByDrivers[s])
    {
      if (res < kv.first.size())
      {
        res = kv.first.size();
      }
    }
    
    return res;
  }
  
  /// Draw cells that will migrate
  ///
  /// @param s Anatomical site
  /// @param migrationCount Number of cells to draw
  IntVector draw(int s, int migrationCount)
  {
    std::uniform_real_distribution<> unif(0, 1);
    
    DoubleVector proportions(_extantCellsByDrivers[s].size(), 0);
    IntVector maxResult(_extantCellsByDrivers[s].size(), 0);
    int idx = 0;
    double sum = 0;
    for (const auto& kv : _extantCellsByDrivers[s])
    {
      maxResult[idx] = kv.second.size();
      proportions[idx] = kv.second.size();
      sum += proportions[idx];
      ++idx;
    }
    for (idx = 0; idx < proportions.size(); ++idx)
    {
      proportions[idx] /= sum;
    }
    
    // make proportions cumulative
    for (idx = 1; idx < proportions.size(); ++idx)
    {
      proportions[idx] += proportions[idx - 1];
    }
    
    IntVector result(_extantCellsByDrivers[s].size(), 0);
    while (migrationCount > 0)
    {
      double r = unif(g_rng);
      for (int j = 0; j < proportions.size(); ++j)
      {
        if (r < proportions[j] || j == proportions.size() - 1)
        {
          if (result[j] < maxResult[j])
          {
            // we have to be careful not to draw more than is allowed
            ++result[j];
            --migrationCount;
          }
          break;
        }
      }
    }
    
    return result;
  }
  
  /// Return the set of observable mutations (in target anatomical site)
  /// of given migration edge
  ///
  /// @param a Arc in migration graph
  IntSet getMutationsOfMigrationEdge(Arc a) const
  {
    IntSet res;
    
    Node v_t = _G.target(a);
    const int t = _anatomicalSiteMap[v_t];
    
    const std::pair<Cell, IntSet>& cellDrivers = _arcToCell[a];
    const Cell& cell = cellDrivers.first;
    IntSet mutations = cellDrivers.second;
    
    int i = cell.getMutation();
    if (_cellTreeMutationToObservableMutation.count(i) == 1)
    {
      int observed_i = _cellTreeMutationToObservableMutation.find(i)->second;
      if (_sampledMutationsByAnatomicalSite[t].count(observed_i) == 1)
      {
        mutations.insert(cell.getPassengerMutations().begin(),
                         cell.getPassengerMutations().end());
        
        IntSet mappedMutations;
        for (int j : mutations)
        {
          if (_cellTreeMutationToObservableMutation.count(i))
          {
            int mapped_j = _cellTreeMutationToObservableMutation.find(j)->second;
            if (_sampledMutationsByAnatomicalSite[t].count(mapped_j))
            {
              mappedMutations.insert(mapped_j);
            }
          }
        }
        
        std::set_intersection(mappedMutations.begin(),
                              mappedMutations.end(),
                              _sampledMutationsByAnatomicalSite[t].begin(),
                              _sampledMutationsByAnatomicalSite[t].end(),
                              std::inserter(res, res.begin()));
        
        
      }
    }

    return res;
  }
  
  /// Finish the vertex labeling
  ///
  /// @param T Clone tree
  /// @param v_i Node
  void finishLabeling(const CloneTree& T,  Node v_i);
  
  /// Determine the set of descendant leaves corresponding to a migration edge of G
  ///
  /// @param arcsGtoVerticesNewT Maps arc of G to leaves of T in same anatomical
  /// site that descend from the seeding cell
  /// @param parentArcMap Maps arc of G to parent arc
  /// @param arcsGtoCumVerticesNewT Maps arc of G to all leaves of T that descend from the seeding cell
  /// @param v_s Node of G
  void computeCumVerticesNewT(const Digraph::ArcMap<NodeSet>& arcsGtoVerticesNewT,
                              const Digraph::ArcMap<Arc>& parentArcMap,
                              Digraph::ArcMap<NodeSet>& arcsGtoCumVerticesNewT,
                              Node v_s) const;
  
  /// Compute parent arc based on mutations
  ///
  /// @param arcsGtoVerticesNewT Maps arc of G to leaves of T in same anatomical
  /// site that descend from the seeding cell
  /// @param Maps arc of G to parent arc
  bool computeParent(const Digraph::ArcMap<NodeSet>& arcsGtoVerticesNewT,
                     Digraph::ArcMap<Arc>& parentArcMap) const;
  
private:
  /// Migration graph;
  Digraph _G;
  /// Migration graph root node
  Node _rootG;
  /// Node to anatomical site map
  IntNodeMap _anatomicalSiteMap;
  /// Anatomical site to node
  NodeVector _indexToVertexG;
  
  /// Seeding cells
  CellMatrix _seedingCells;
  /// Migration arc to seeding cell
  CellArcMap _arcToCell;

  /// Generation
  int _generation;
  
  /// Extant cells per anatomical site further split by drivers
  ClonalComposition _extantCellsByDrivers;
  /// Number of extant cells
  int _nrExtantCells;
  /// Number of mutations
  int _nrMutations;
  
  /// Number of anatomical sites
  int _nrAnatomicalSites;
  /// Number of active anatomical site
  int _nrActiveAnatomicalSites;
  /// Is active anatomical site
  BoolVector _isActiveAnatomicalSite;
  /// Anatomical site factors
  AnatomicalSiteFactorVector _anatomicalSiteFactors;
  /// Anatomical site label
  StringVector _anatomicalSiteLabel;
  
  /// Carrying capacity for each site
  const double _K;
  /// Migration rate
  const double _migrationRate;
  /// Driver mutation probability
  const double _driverProb;
  /// Mutation frequency threshold
  const double _mutFreqThreshold;
  /// Maximum number of anatomical sites
  const int _maxNrAnatomicalSites;
  /// Number of samples per anatomical site
  const int _nrSamplesPerAnatomicalSite;
  /// Number of samples for the primary
  const int _nrSamplesPrimary;
  /// Target coverage
  const int _targetCoverage;
  /// Migration pattern
  const Pattern _pattern;
  /// Sequencing error rate
  const double _seqErrorRate;
  /// Purity
  const double _purity;
  
  /// _populationRecord[s][i] is the
  /// number of cells at anatomical site s in generation i
  IntMatrix _populationRecord;
  /// Set of driver mutations
  IntSet _driverMutations;
  
  /// Observable mutation to cell-tree mutation
  IntVector _observableMutationToCellTreeMutation;
  /// Cell-tree mutation to observable mutation
  IntIntMap _cellTreeMutationToObservableMutation;
  /// Clone tree
  CloneTree* _pCloneT;
  /// Proportions of clone tree leaves
  DoubleNodeMap* _pAnatomicalSiteProportions;
  /// Proportions of clone tree leaves
  DoubleVectorNodeMap* _pSampleProportions;
  /// Mutations present in leaves of _pCloneT
  IntSetNodeMap* _pMutations;
  /// Anatomical site of leaves of _pCloneT
  IntNodeMap* _pAnatomicalSite;
  /// True frequencies
  DoubleTensor _freq;
  /// Variant reads
  IntTensor _var;
  /// Reference reads
  IntTensor _ref;
  /// Sampled leaves
  NodeSet _sampledLeaves;
  
  /// Sampled clone tree
  CloneTree* _pCloneTT;
  /// Usages of clone tree leaves of _pCloneTT
  DoubleVectorNodeMap* _pSampleProportionsTT;
  /// Mutations present in leaves of _pCloneTT
  IntSetNodeMap* _pMutationsTT;
  /// Anatomical site of leaves of _pCloneTT
  IntNodeMap* _pAnatomicalSiteTT;
  /// Sampled mutations
  IntSet _sampledMutations;
  /// Sampled mutations by anatomical sites
  IntSetVector _sampledMutationsByAnatomicalSite;
  /// Sampled driver mutations
  IntSet _sampledDriverMutations;
  /// Vertex labeling
  StringNodeMap* _pVertexLabelingTT;
  /// Migration graph arcs to extant clones in _pCloneTT
  Digraph::ArcMap<NodeSet> _arcGtoVerticesTT;
  /// Migration graph arcs to sampled mutations (ancestral clone)
  Digraph::ArcMap<IntSet> _arcGtoSampledMutations;
};

#endif // SIMULATION_H
