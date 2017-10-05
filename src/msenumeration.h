/*
 *  msenumeration.h
 *
 *   Created on: 17-aug-2017
 *       Author: M. El-Kebir
 */

#ifndef MSENUMERATION_H
#define MSENUMERATION_H

#include "frequencymatrix.h"
#include "migrationtree.h"
#include "clonetree.h"
#include "spruce/realtensor.h"
#include "spruce/solutionset.h"

class MSEnumeration
{
public:
  MSEnumeration(const FrequencyMatrix& F,
                const std::string& primary,
                const std::string& outputDirectory,
                const StringToIntMap& colorMap);
  
  ~MSEnumeration()
  {
    for (SolutionEntry* pEntry : _solutions)
    {
      delete pEntry;
    }
  }
  
  struct SolutionEntry
  {
  public:
    SolutionEntry(const Digraph& tree,
                  Node root,
                  const StringNodeMap& idMap,
                  const StringNodeMap& leafLabelMap)
      : _T(tree, root, idMap, leafLabelMap)
      , _lPlus(_T.tree())
      , _usage(_T.tree())
      , _frequency(_T.tree())
      , _characterNodeMap(_T.tree())
    {
    }

    void writeDOT(std::ostream& out, const StringToIntMap& colorMap) const
    {
      _T.writeDOT(out, _lPlus, colorMap, _usage, _characterNodeMap);
    }
    
    void writeDOT(std::ostream& out) const
    {
      StringToIntMap colorMap = _T.generateColorMap();
      _T.writeDOT(out, _lPlus, colorMap, _usage, _characterNodeMap);
    }
    
    CloneTree _T;
    StringNodeMap _lPlus;
    DoubleVectorNodeMap _usage;
    DoubleVectorNodeMap _frequency;
    StringNodeMap _characterNodeMap;
  };
  
  int run(const MigrationTree& migrationTree,
          bool force_mS,
          const std::string& migrationTreeString,
          int limit);

  void run(bool force_mS,
           int limit);
  
  int getNrSolutions() const
  {
    return _solutions.size();
  }
  
  const SolutionEntry& getSolution(int index) const
  {
    assert(0 <= index && index < _solutions.size());
    return *_solutions[index];
  }
  
private:
  /// Initialize frequency tensors (F^-, F+^) with an additional
  /// multi-state character
  void initRealTensors();
  
  /// Initializes state tree vector for first n characters
  void initStateTrees();
  
  SolutionEntry* constructCloneTree(const gm::Solution& solution);
  
  void labelVerticesByAnatomicalSites(const CloneTree& T,
                                      Node v,
                                      const std::string& sStr,
                                      StringNodeMap& lPlus) const;
  
  typedef std::vector<SolutionEntry*> SolutionVector;
  
private:
  /// Frequency matrix (frequencies of two-state characters)
  const FrequencyMatrix& _F;
  /// Primary
  const std::string& _primary;
  /// Output directory
  const std::string& _outputDirectory;
  /// Color map
  const StringToIntMap& _colorMap;
  /// Number of anatomical sites
  const int _m;
  /// Number of samples
  const int _k;
  /// Number of characters
  const int _n;
  /// Anatomical sites
  StringSet _Sigma;
  /// Frequency tensor with lowerbounds
  RealTensor _F_lb;
  /// Frequency tensor with upperbounds
  RealTensor _F_ub;
  /// State tree vector
  StateTreeVector _S;
  /// Migration tree list
  MigrationTreeList _barS;
  /// Solution clone trees
  SolutionVector _solutions;
};

#endif // MSENUMERATION_H
