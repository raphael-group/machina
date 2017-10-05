/*
 * basematrix.cpp
 *
 *  Created on: 5-sep-2017
 *      Author: M. El-Kebir
 */

#include "basematrix.h"

BaseMatrix::BaseMatrix()
  : _m(0)
  , _k(0)
  , _n(0)
  , _indexToAnatomicalSite()
  , _anatomicalSiteToIndex()
  , _indexToSample()
  , _sampleToIndex()
  , _indexToCharacter()
  , _characterToIndex()
  , _sampleIndexToAnatomicalSiteIndex()
  , _anatomicalSiteIndexToSampleIndices()
{
}
