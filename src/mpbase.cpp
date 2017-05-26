/*
 * mpbase.cpp
 *
 *  Created on: 10-jan-2017
 *      Author: M. El-Kebir
 */

#include "mpbase.h"

MPBase::MPBase(const CharacterTree& T)
  : _T(T)
  , _solutionIndex(_T.getNrCharacters())
  , _stateVector(_T.getNrCharacters())
  , _homoplasyVector(_T.getNrCharacters())
{
}

void MPBase::updateHomoplasy()
{
  const int nrCharacters = _T.getNrCharacters();
  
  for (int c = 0; c < nrCharacters; ++c)
  {
    int solIdx = _solutionIndex[c];
    _homoplasyVector[c][solIdx] = false;
    
    // update homoplasy
    const int nrStates = _T.getNrStates(c);
    for (int s = 0; s < nrStates; ++s)
    {
      NodeSet nodes;
      for (NodeIt u(_T.tree()); u != lemon::INVALID; ++u)
      {
        if ((*_stateVector[c][solIdx])[u] == s)
        {
          nodes.insert(u);
        }
      }
      if (!nodes.empty())
      {
        _homoplasyVector[c][solIdx] = _homoplasyVector[c][solIdx] || (!_T.isConnected(nodes));
      }
    }
  }
}

void MPBase::writeDOT(std::ostream& out) const
{
  const Digraph& tree = _T.tree();
  const int n = _T.getNrCharacters();
  
  out << "digraph T {" << std::endl;
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt u(tree); u != lemon::INVALID; ++u)
  {
    if (_T.isLeaf(u))
    {
      out << "\t\t" << tree.id(u) << " [label=<" << _T.label(u) << "<BR/>";
      
      for (int c = 0; c < n; ++c)
      {
        int solIdx = _solutionIndex[c];
        if (_homoplasyVector[c][solIdx])
        {
//          out << "<font color=\"red\">" << _stateVector[u][i] << " </font>";
          out << (*_stateVector[c][solIdx])[u] << " ";
        }
        else
        {
//          out << "<font color=\"black\">" << _stateVector[u][i] << " </font>";
          out << (*_stateVector[c][solIdx])[u] << " ";
        }
      }
      
      out << ">]" << std::endl;
    }
  }
  out << "\t}" << std::endl;
  
  for (NodeIt u(tree); u != lemon::INVALID; ++u)
  {
    if (!_T.isLeaf(u))
    {
      out << "\t" << tree.id(u) << " [label=<" << _T.label(u) << "<BR/>";
      
      for (int c = 0; c < n; ++c)
      {
        int solIdx = _solutionIndex[c];
        if (_homoplasyVector[c][solIdx])
        {
//          out << "<font color=\"red\">" << _stateVector[u][i] << " </font>";
          out << (*_stateVector[c][solIdx])[u] << " ";
        }
        else
        {
//          out << "<font color=\"black\">" << _stateVector[u][i] << " </font>";
          out << (*_stateVector[c][solIdx])[u] << " ";
        }
      }
      
      out << ">]" << std::endl;
    }
  }
  
  for (ArcIt a(tree); a != lemon::INVALID; ++a)
  {
    Node u = tree.source(a);
    Node v = tree.target(a);
    
    out << "\t" << tree.id(tree.source(a)) << " -> " << tree.id(tree.target(a)) << " [label=<";
    
    bool first = true;
    for (int c = 0; c < n; ++c)
    {
      if (_T.getNrStates(c) == 2)
      {
        int solIdx = _solutionIndex[c];
        if (_homoplasyVector[c][solIdx] && (*_stateVector[c][solIdx])[u] != (*_stateVector[c][solIdx])[v])
        {
          if (first)
            first = false;
          else
            out << "<BR/>";
          
          out << "<font color=\"";
          
          if ((*_stateVector[c][solIdx])[v] == 1)
          {
            out << "blue\">";
          }
          else
          {
            out << "red\">";
          }
          out << _T.characterLabel(c) << (((*_stateVector[c][solIdx])[v] == 1) ? "+" : "-") << "</font>";
        }
        else
        {
          out << _T.characterLabel(c);
        }
      }
    }
    out << ">]" << std::endl;
  }
  
  out << "}" << std::endl;
}

bool MPBase::nextSolution()
{
  const int nrCharacters = _T.getNrCharacters();
  
  for (int c = 0; c < nrCharacters; ++c)
  {
    if (_solutionIndex[c] != _stateVector[c].size() - 1)
    {
      ++_solutionIndex[c];
      for (int d = 0; d < c; ++d)
      {
        _solutionIndex[d] = 0;
      }
      return true;
    }
  }
  
  return false;
}
