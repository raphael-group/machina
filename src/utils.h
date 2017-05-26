#ifndef UTILS_H
#define UTILS_H

#include <lemon/list_graph.h>
#include <lemon/adaptors.h>
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
DIGRAPH_TYPEDEFS(Digraph);
typedef Digraph::NodeMap<Digraph::Node> NodeNodeMap;
typedef Digraph::NodeMap<std::string> StringNodeMap;
typedef std::map<std::string, Node> StringToNodeMap;
typedef std::map<std::string, int> StringToIntMap;
typedef std::set<Node> NodeSet;
typedef std::set<Arc> ArcSet;
typedef Digraph::NodeMap<NodeSet> NodeNodeSetMap;
typedef std::pair<int, Node> IntNodePair;
typedef std::list<IntNodePair> IntNodePairList;
typedef std::vector<Node> NodeVector;
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
typedef std::vector<int> IntVector;
typedef std::vector<bool> BoolVector;
typedef std::vector<double> DoubleVector;
typedef std::vector<DoubleVector> DoubleMatrix;
typedef std::set<int> IntSet;
typedef std::list<std::string> StringList;
typedef std::set<std::string> StringSet;
typedef std::pair<std::string, std::string> StringPair;
typedef std::list<StringPair> StringPairList;
typedef Digraph::NodeMap<StringList> StringListNodeMap;
typedef Digraph::NodeMap<BoolVector> BoolVectorNodeMap;
typedef Digraph::NodeMap<IntVector> IntVectorNodeMap;
typedef Digraph::NodeMap<DoubleVector> DoubleVectorNodeMap;
typedef Digraph::NodeMap<IntSet> IntSetNodeMap;

typedef lemon::SubDigraph<const Digraph> SubDigraph;
typedef SubDigraph::ArcIt SubArcIt;
typedef SubDigraph::NodeIt SubNodeIt;
typedef SubDigraph::OutArcIt SubOutArcIt;
typedef SubDigraph::InArcIt SubInArcIt;

std::istream& getline(std::istream& is, std::string& t);
bool parseMigrationGraph(const std::string& migrationGraphFile,
                         const StringSet& samples,
                         StringPairList& forcedComigrations);

extern std::mt19937 g_rng;

#endif // UTILS_H
