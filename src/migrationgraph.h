/*
 * migrationgraph.h
 *
 *  Created on: 24-oct-2016
 *      Author: M. El-Kebir
 */

#ifndef MIGRATIONGRAP_H
#define MIGRATIONGRAP_H

#include "utils.h"
#include "clonetree.h"
#include <lemon/connectivity.h>

/// This class models a migration graph
class MigrationGraph
{
public:
  /// Migration pattern
  enum Pattern
  {
    /// Parallel single-source seeding (PS)
    PS,
    /// Single-source seeding (S)
    S,
    /// Mutli-source seeding (M)
    M,
    /// Reseeding (R)
    R
  };
  
  /// Number of migration patterns
  static const int _nrPatterns;
  
  /// Default constructor
  MigrationGraph();
  
  /// Constructor
  ///
  /// @param G Directed graph
  /// @param root Root node
  /// @param id Node identifier
  MigrationGraph(const Digraph& G,
                 Node root,
                 const StringNodeMap& id);
  
  /// Copy constructor
  MigrationGraph(const MigrationGraph& other);
  
  /// Construct migration graph given a non-binary clone tree
  /// and vertex labeling
  ///
  /// @param T Clone tree
  /// @param lPlus Vertex labeling
  MigrationGraph(const CloneTree& T,
                 const StringNodeMap& lPlus);
  
  /// Return anatomical site label of given node
  const std::string& l(Node x) const
  {
    return _nodeToId[x];
  }
  
  /// Return the labels of all anatomical sites
  StringSet getAnatomicalSites() const
  {
    StringSet res;
    for (const auto& kv : _idToNode)
    {
      res.insert(kv.first);
    }
    return res;
  }
  
  /// Return the number of anatomical sites
  int getNrAnatomicalSites() const
  {
    return _idToNode.size();
  }
  
  /// Return whether the migration graph has a directed cycle
  bool hasReseeding() const
  {
    return !lemon::dag(_G);
  }
  
  /// Return whether the migration graph has a PS or S topology
  bool isSingleSourceSeeded() const
  {
    for (NodeIt y(_G); y != lemon::INVALID; ++y)
    {
      Node parent = lemon::INVALID;
      for (InArcIt a(_G, y); a != lemon::INVALID; ++a)
      {
        Node x = _G.source(a);
        if (parent == lemon::INVALID)
          parent = x;
        else if (parent != x)
          return false;
      }
    }
    
    return true;
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
  Pattern getPattern() const;
  
  /// Return a short string corresponding to the migration pattern
  /// of the migration graph
  ///
  /// @param pattern Migration pattern
  /// @param monoclonal Monoclonal
  static std::string getPatternString(Pattern pattern,
                                      bool monoclonal);
  
  /// Return a short string corresponding to the migration pattern
  /// of the migration graph
  ///
  /// @param pattern Migration pattern
  static std::string getPatternString(Pattern pattern);
  
  /// Return a long string corresponding to the migration pattern
  /// of the migration graph
  ///
  /// @param pattern Migration pattern
  /// @param monoclonal Monoclonal
  static std::string getPatternLongString(Pattern pattern,
                                          bool monoclonal);
  
  /// Return a short string corresponding to the allowed migration patterns
  /// of the migration graph
  ///
  /// @param pattern Migration pattern
  static std::string getAllowedPatternsString(Pattern pattern);
  
  /// Return the number of seeding sites
  int getNrSeedingSites() const
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
  
  /// Return whether migration graph is monoclonal
  bool isMonoclonal() const
  {
    for (NodeIt x(_G); x != lemon::INVALID; ++x)
    {
      BoolNodeMap visited(_G, false);
      for (OutArcIt a(_G, x); a != lemon::INVALID; ++a)
      {
        Node y = _G.target(a);
        if (visited[y])
          return false;
        else
          visited[y] = true;
      }
    }
    
    return true;
  }
  
  /// Return whether migration graph is polyclonal
  bool isPolyclonal() const
  {
    return !isMonoclonal();
  }
  
  /// Return the number of anatomical sites with non-unique parentage
  int getNrNonUniqueParentageSites() const
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

inline MigrationGraph::Pattern MigrationGraph::getPattern() const
{
  if (getNrSeedingSites() == 1)
  {
    return PS;
  }
  else if (isSingleSourceSeeded() && lemon::dag(_G))
  {
    return S;
  }
  else if (lemon::dag(_G))
  {
    return M;
  }
  else
  {
    return R;
  }
}

inline std::string MigrationGraph::getPatternString(Pattern pattern)
{
  switch (pattern)
  {
    case PS:
      return "PS";
    case S:
      return "S";
    case M:
      return "M";
    case R:
      return "R";
    default:
      assert(false);
      return "ERROR";
  }
}

inline std::string MigrationGraph::getPatternString(Pattern pattern,
                                                    bool monoclonal)
{
  std::string res = monoclonal ? "m" : "p";
  switch (pattern)
  {
    case PS:
      return res + "PS";
    case S:
      return res + "S";
    case M:
      return res + "M";
    case R:
      return res + "R";
    default:
      assert(false);
      return "ERROR";
  }
}

inline std::string MigrationGraph::getAllowedPatternsString(Pattern pattern)
{
  switch (pattern)
  {
    case PS:
      return "PS";
    case S:
      return "PS, S";
    case M:
      return "PS, S, M";
    case R:
      return "PS, S, M, R";
    default:
      assert(false);
      return "ERROR";
  }
}

inline std::string MigrationGraph::getPatternLongString(Pattern pattern,
                                                        bool monoclonal)
{
  std::string res = monoclonal ? "monoclonal " : "polyclonal ";
  switch (pattern)
  {
    case PS:
      return res + "parallel single-source seeding";
    case S:
      return res + "single-source seeding";
    case M:
      return res + "multi-source seeding";
    case R:
      return res + "reseeding";
    default:
      assert(false);
      return "ERROR";
  }
}

#endif // MIGRATIONGRAPH_H
