/*
 * sankofflabeling.cpp
 *
 *  Created on: 12-jan-2017
 *      Author: M. El-Kebir
 */

#include "sankofflabeling.h"

void SankoffLabeling::run()
{
  assert(_stateToSample[0] == _primary);
  Sankoff sankoff(_charT);
  sankoff.run(0);
  
  do
  {
    StringNodeMap* pLabeling = new StringNodeMap(_T.tree());
    for (NodeIt v(_charT.tree()); v != lemon::INVALID; ++v)
    {
      const Node vv = _T.getNodeByLabel(_charT.label(v));
      assert(vv != lemon::INVALID);
      
      const int s = sankoff.state(v, 0);
      assert(0 <= s && s < _stateToSample.size());
      
      const std::string& sample = _stateToSample[s];
      pLabeling->set(vv, sample);
    }
    _mpLabelings.push_back(pLabeling);
  } while (sankoff.nextSolution());
}
