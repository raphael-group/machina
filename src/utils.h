/*
 * utils.h
 *
 *  Created on: 19-oct-2016
 *      Author: M. El-Kebir
 */

#ifndef UTILS_H
#define UTILS_H

#include <lemon/list_graph.h>
#include <lemon/adaptors.h>
#include <lemon/tolerance.h>
#include <cassert>
#include <iostream>
#include <set>
#include <algorithm>
#include <list>
#include <map>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <random>

typedef lemon::ListDigraph Digraph;
typedef lemon::ListGraph Graph;
DIGRAPH_TYPEDEFS(Digraph);
typedef Digraph::NodeMap<Digraph::Node> NodeNodeMap;
typedef Digraph::NodeMap<std::string> StringNodeMap;
typedef std::map<std::string, Node> StringToNodeMap;
typedef std::map<std::string, int> StringToIntMap;
typedef std::set<Node> NodeSet;
typedef std::set<Arc> ArcSet;
typedef std::vector<Arc> ArcVector;
typedef Digraph::NodeMap<NodeSet> NodeNodeSetMap;
typedef std::pair<int, Node> IntNodePair;
typedef std::list<IntNodePair> IntNodePairList;
typedef std::vector<Node> NodeVector;
typedef std::vector<NodeVector> NodeMatrix;
typedef std::list<Node> NodeList;
typedef NodeList::const_iterator NodeListIt;
typedef std::vector<NodeList> NodeListVector;
typedef NodeListVector::const_iterator NodeListVectorIt;
typedef std::list<NodeListIt> NodeListItList;
typedef std::vector<NodeListIt> NodeListItVector;
typedef std::list<NodeList> NodeListList;
typedef std::vector<std::string> StringVector;
typedef std::pair<Node, Node> NodePair;
typedef std::list<NodePair> NodePairList;
typedef Digraph::NodeMap<NodePair> NodePairMap;
typedef std::set<std::string> StringSet;
typedef std::multiset<StringSet> Split;
typedef std::set<Split> SplitSet;
typedef std::set<StringSet> StringSetSet;
typedef std::vector<int> IntVector;
typedef std::vector<IntVector> IntMatrix;
typedef std::vector<IntMatrix> IntTensor;
typedef std::vector<bool> BoolVector;
typedef std::vector<BoolVector> BoolMatrix;
typedef std::vector<double> DoubleVector;
typedef std::vector<DoubleVector> DoubleMatrix;
typedef std::vector<DoubleMatrix> DoubleTensor;
typedef std::set<int> IntSet;
typedef std::list<int> IntList;
typedef IntSet::const_iterator IntSetIt;
typedef std::list<std::string> StringList;
typedef std::set<std::string> StringSet;
typedef std::pair<std::string, std::string> StringPair;
typedef std::list<StringPair> StringPairList;
typedef Digraph::NodeMap<StringList> StringListNodeMap;
typedef Digraph::NodeMap<BoolVector> BoolVectorNodeMap;
typedef Digraph::NodeMap<IntVector> IntVectorNodeMap;
typedef Digraph::NodeMap<DoubleVector> DoubleVectorNodeMap;
typedef Digraph::NodeMap<IntSet> IntSetNodeMap;
typedef Digraph::ArcMap<IntSet> IntSetArcMap;
typedef std::pair<int, int> IntPair;
typedef std::set<IntPair> IntPairSet;
typedef std::vector<IntSet> IntSetVector;
typedef Digraph::NodeMap<IntPair> IntPairNodeMap;
typedef Digraph::NodeMap<IntPairSet> IntPairSetNodeMap;
typedef std::map<IntPair, Node> IntPairToNodeMap;
typedef std::map<IntPair, NodeSet> IntPairToNodeSetMap;
typedef std::pair<int, IntPair> IntTriple;

typedef lemon::SubDigraph<const Digraph> SubDigraph;
typedef SubDigraph::ArcIt SubArcIt;
typedef SubDigraph::NodeIt SubNodeIt;
typedef SubDigraph::OutArcIt SubOutArcIt;
typedef SubDigraph::InArcIt SubInArcIt;

/// Verbosity level
typedef enum
{
  /// No output
  VERBOSE_NONE,
  /// Only essential information is output
  VERBOSE_ESSENTIAL,
  /// Essential and non-essential information is output
  VERBOSE_NON_ESSENTIAL,
  /// Essential, non-essential and debug information is output
  VERBOSE_DEBUG
} VerbosityLevel;

/// Verbosity level
extern VerbosityLevel g_verbosity;

/// Get line from stream in a platform independent manner
std::istream& getline(std::istream& is, std::string& t);

/// Random number generator
extern std::mt19937 g_rng;

/// Tolerance for floating point comparisons
extern lemon::Tolerance<double> g_tol;

/// Current line number
extern int g_lineNumber;

/// Get current line number as a string
std::string getLineNumber();

/// Write a matrix of doubles to the specified output stream
///
/// @param out Output stream
/// @param M Matrix of doubles
std::ostream& operator<<(std::ostream& out, const DoubleMatrix& M);

/// Read a matrix of doubles from the specified input stream
///
/// @param in Input stream
/// @param M Matrix of doubles
std::istream& operator>>(std::istream& in, DoubleMatrix& M);

/// Write a tensor of doubles to the specified output stream
///
/// @param out Output stream
/// @param M Tensor of doubles
std::ostream& operator<<(std::ostream& out, const DoubleTensor& M);

/// Read a tensor of doubles from the specified output stream
///
/// @param in Input stream
/// @param M Tensor of doubles
std::istream& operator>>(std::istream& in, DoubleTensor& M);

#endif // UTILS_H
