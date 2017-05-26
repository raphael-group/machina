/*
 * migrationgraph.h
 *
 *  Created on: 24-oct-2016
 *      Author: M. El-Kebir
 */

#ifndef SAMPLEGRAPH_H
#define SAMPLEGRAPH_H

#include "utils.h"
#include "clonetree.h"
#include "nonbinaryclonetree.h"
#include <lemon/connectivity.h>

/// This class models a migration graph
class MigrationGraph
{
public:
  enum MigrationPattern
  {
    PMSSS,
    PSSS,
    SSS,
    MSS,
    RS
  };
  
  /// Default constructor
  MigrationGraph();
  
  /// Copy constructor
  MigrationGraph(const MigrationGraph& other);
  
  /// Construct migration graph given a binary clone tree and vertex labeling
  ///
  /// @param T Binary clone tree
  /// @param lPlus Vertex labeling
  MigrationGraph(const CloneTree& T,
                 const StringNodeMap& lPlus);
  
  /// Construct migration graph given a non-binary clone tree
  /// and vertex labeling
  ///
  /// @param T Non-binary clone tree
  /// @param lPlus Vertex labeling
  MigrationGraph(const NonBinaryCloneTree& T,
                 const StringNodeMap& lPlus);
  
  /// Return anatomical site label of given node
  const std::string& l(Node x) const
  {
    return _nodeToId[x];
  }
  
  /// Return the labels of all anatomical sites
  StringSet getSamples() const
  {
    StringSet res;
    for (const auto& kv : _idToNode)
    {
      res.insert(kv.first);
    }
    return res;
  }
  
  /// Return the number of anatomical sites
  int getNrSamples() const
  {
    return _idToNode.size();
  }
  
  /// Return whether the migration graph has a directed cycle
  bool hasReseeding() const
  {
    return !lemon::dag(_G);
  }
  
  /// Return the number of migrations in the migration graph
  int getNrMigrations() const
  {
    return lemon::countArcs(_G);
  }
  
  /// Return the number of comigrations in the migration graph
  ///
  /// @param T Clone tree
  /// @param lPlus Vertex labeling
  template <class CLONETREE>
  int getNrComigrations(const CLONETREE& T,
                        const StringNodeMap& lPlus) const;
  
  /// Return the migration pattern of the migration graph
  ///
  /// @param T Clone tree
  /// @param lPlus Vertex labeling
  template <class CLONETREE>
  MigrationPattern getMigrationPattern(const CLONETREE& T,
                                       const StringNodeMap& lPlus) const;
  
  /// Return a string corresponding to the migration pattern
  /// of the migration graph
  ///
  /// @param pattern Migration pattern
  static std::string getMigrationPatternString(MigrationPattern pattern);
  
  /// Return the number of seeding sites
  int getNrSeedingSamples() const
  {
    int res = 0;
    for (NodeIt x(_G); x != lemon::INVALID; ++x)
    {
      if (OutArcIt(_G, x) != lemon::INVALID)
      {
        ++res;
      }
    }
    return res;
  }
  
  /// Return the number of anatomical sites with non-unique parentage
  int getNrNonUniqueParentageSamples() const
  {
    int res = 0;

    if (InArcIt(_G, _root) != lemon::INVALID)
    {
      ++res;
    }
    
    for (NodeIt y(_G); y != lemon::INVALID; ++y)
    {
      if (y == _root)
        continue;
      
      NodeSet sourceSet;
      for (InArcIt a(_G, y); a != lemon::INVALID; ++a)
      {
        Node x = _G.source(a);
        sourceSet.insert(x);
      }
      
      if (sourceSet.size() > 1)
      {
        ++res;
      }
    }
    
    return res;
  }
  
  /// Deserialize
  ///
  /// @param in Input stream
  bool read(std::istream& in);
  
  /// Serialize
  ///
  /// @param out Output stream
  void write(std::ostream& out) const;
  
  /// Print migration graph in DOT format
  ///
  /// @param out Output stream
  void writeDOT(std::ostream& out) const;
  
  /// Print migration graph in DOT format using the given color map
  ///
  /// @param out Output stream
  /// @param colorMap Color map
  void writeDOT(std::ostream& out,
                const StringToIntMap& colorMap) const;
  
  /// Return whether the migration graph is connected (it should be!)
  bool isConnected() const;

private:
  /// Graph
  Digraph _G;
  /// Root
  Node _root;
  /// Vertex label
  StringNodeMap _nodeToId;
  /// Get vertex by label
  StringToNodeMap _idToNode;
};

template <class CLONETREE>
inline int MigrationGraph::getNrComigrations(const CLONETREE& T,
                                          const StringNodeMap& lPlus) const
{
  const Digraph& TT = T.tree();
  
  lemon::DynArcLookUp<Digraph> arcLookUp(_G);
  
  ArcSet migrationEdges = T.getMigrationEdges(lPlus);
  
  int nrComigrations = 0;
  for (NodeIt x(_G); x != lemon::INVALID; ++x)
  {
    const std::string& l_x = l(x);
    for (NodeIt y(_G); y != lemon::INVALID; ++y)
    {
      const std::string& l_y = l(y);
      if (arcLookUp(x, y) != lemon::INVALID)
      {
        // 1. collect all migration edges between (x,y)
        //        ++nrComigrations;
        
        ArcSet migrationEdges_xy;
        for (Arc uv : migrationEdges)
        {
          Node u = TT.source(uv);
          Node v = TT.target(uv);
          
          if (lPlus[u] == l_x && lPlus[v] == l_y)
          {
            migrationEdges_xy.insert(uv);
          }
        }
        
        // 2. now for each migration edge (u,v) count the number of ancestral migration edges it has
        int nrAncestralMigrationEdges = 0;
        for (Arc uv : migrationEdges_xy)
        {
          int count = 0;
          
          Arc prevArc = uv;
          while (prevArc != lemon::INVALID)
          {
            if (migrationEdges_xy.count(prevArc) == 1)
            {
              ++count;
            }
            prevArc = InArcIt(TT, TT.source(prevArc));
          }
          if (count > nrAncestralMigrationEdges)
          {
            nrAncestralMigrationEdges = count;
          }
        }
        nrComigrations += nrAncestralMigrationEdges;
      }
    }
  }
  
  assert(nrComigrations >= (_idToNode.size() - 1));
  return nrComigrations;
}

template <class CLONETREE>
inline MigrationGraph::MigrationPattern MigrationGraph::getMigrationPattern(const CLONETREE& T,
                                                                      const StringNodeMap& lPlus) const
{
  if (getNrSeedingSamples() == 1)
  {
    if (getNrMigrations() == getNrSamples() - 1)
    {
      return PMSSS;
    }
    else
    {
      return PSSS;
    }
  }
  else if (getNrComigrations(T, lPlus) == getNrSamples() - 1)
  {
    return SSS;
  }
  else if (lemon::dag(_G))
  {
    return MSS;
  }
  else
  {
    return RS;
  }
}

inline std::string MigrationGraph::getMigrationPatternString(MigrationPattern pattern)
{
  switch (pattern)
  {
    case PMSSS:
      return "PMSSS";
    case PSSS:
      return "PSSS";
    case SSS:
      return "SSS";
    case MSS:
      return "MSS";
    case RS:
      return "RS";
    default:
      assert(false);
      return "";
  }
}

#endif // SAMPLEGRAPH_H
