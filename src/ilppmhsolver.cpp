/*
 * ilppmhsolver.cpp
 *
 *  Created on: 25-aug-2017
 *      Author: M. El-Kebir
 */

#include "ilppmhsolver.h"
#include <lemon/time_measure.h>

IlpPmhSolver::IlpPmhSolver(const CloneTree& T,
                           const std::string& primary,
                           MigrationGraph::Pattern pattern,
                           const std::string& gurobiLogFilename,
                           const StringPairList& forcedComigrations)
  : _T(T)
  , _primary(primary)
  , _pattern(pattern)
  , _forcedComigrations(forcedComigrations)
  , _lca()
  , _anatomicalSiteToIndex()
  , _indexToAnatomicalSite()
  , _indexToNode()
  , _pNodeToIndex(NULL)
  , _indexToArc()
  , _pArcToIndex(NULL)
  , _primaryIndex(-1)
  , _L()
  , _env(gurobiLogFilename)
  , _model(_env)
  , _x()
  , _y()
  , _z()
  , _w()
  , _gamma()
  , _sigma()
  , _LB(-1)
  , _UB(-1)
  , _G()
  , _rootG(lemon::INVALID)
  , _nodeGToSubLabel(_G)
  , _subLabelToNodeG()
  , _pLPlus(NULL)
{
}

IlpPmhSolver::~IlpPmhSolver()
{
  delete _pLPlus;
  delete _pNodeToIndex;
  delete _pArcToIndex;
}

void IlpPmhSolver::init(const IntTriple& bounds)
{
  initIndices();
  initVariables();
  initLeafVariables();
  initVertexLabelingConstraints();
  initForcedComigrations();
  initConstraintsG();
  initConstraintsNonEdgesG();
  initConstraints();
  initLeafConstraints();
  initObjective(bounds);
  initCallbacks();
  initWarmStart();
}

void IlpPmhSolver::initWarmStart()
{
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    if (isLeaf(v_i)) continue;
    
    const int i = (*_pNodeToIndex)[v_i];
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      const int size_L_s = _L[s].size();
      for (int c = 0; c < size_L_s; ++c)
      {
        if (s == _primaryIndex && c == 0)
        {
          _x[i][s][c].set(GRB_DoubleAttr_Start, 1);
        }
        else
        {
          _x[i][s][c].set(GRB_DoubleAttr_Start, 0);
        }
      }
    }
  }
}

void IlpPmhSolver::initIndices()
{
  // 1. Initialize _primaryIndex, _anatomicalSiteToIndex
  // and _indexToAnatomicalSite
  _primaryIndex = -1;
  _anatomicalSiteToIndex.clear();
  _indexToAnatomicalSite.clear();

  const StringSet Sigma = getAnatomicalSites();
  for (const std::string& s : Sigma)
  {
    _anatomicalSiteToIndex[s] = _indexToAnatomicalSite.size();
    _indexToAnatomicalSite.push_back(s);
    if (_primary == s)
    {
      _primaryIndex = _indexToAnatomicalSite.size() - 1;
    }
  }
  assert(_primaryIndex != -1);
  
  // 2. Initialize _L and _lca
  _L = NodeSetVector(_indexToAnatomicalSite.size());
  for (NodeIt u(getTree()); u != lemon::INVALID; ++u)
  {
    if (isLeaf(u))
    {
      const std::string& sStr = getLeafAnatomicalSiteLabel(u);
      assert(_anatomicalSiteToIndex.count(sStr) == 1);
      
      int s = _anatomicalSiteToIndex[sStr];
      _L[s].insert(u);
    }
  }
  
  // add dummy vertex
  _L[_primaryIndex].insert(lemon::INVALID);
  
  const int nrAnatomicalSites = _indexToAnatomicalSite.size();
  _lca = NodeVector(nrAnatomicalSites, lemon::INVALID);
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    if (s == _primaryIndex)
      _lca[s] = getRoot();
    else
      _lca[s] = getLCA(_L[s]);
  }
  
  // 3. Initialize _indexToNode and _nodeToIndex
  const int nrNodes = lemon::countNodes(getTree());
  _indexToNode = NodeVector(nrNodes, lemon::INVALID);
  delete _pNodeToIndex;
  _pNodeToIndex = new IntNodeMap(getTree(), -1);
  
  int idx = 0;
  for (NodeIt v(getTree()); v != lemon::INVALID; ++v)
  {
    _indexToNode[idx] = v;
    _pNodeToIndex->set(v, idx);
    ++idx;
  }
  
  // 4. Initialize _indexToArc and _arcToIndex
  const int nrArcs = lemon::countArcs(getTree());
  _indexToArc = ArcVector(nrArcs, lemon::INVALID);
  delete _pArcToIndex;
  _pArcToIndex = new IntArcMap(getTree(), -1);
  
  idx = 0;
  for (ArcIt a(getTree()); a != lemon::INVALID; ++a)
  {
    _indexToArc[idx] = a;
    _pArcToIndex->set(a, idx);
    ++idx;
  }
}

IntTriple IlpPmhSolver::run(const CloneTree& T,
                            const std::string& primary,
                            const std::string& outputDirectory,
                            const std::string& outputPrefix,
                            const StringToIntMap& colorMap,
                            MigrationGraph::Pattern pattern,
                            int nrThreads,
                            bool outputILP,
                            bool outputSearchGraph,
                            int timeLimit,
                            const IntTriple& bounds,
                            const StringPairList& forcedComigrations)
{
  std::string filenameGurobiLog;
  if (!outputDirectory.empty())
  {
    char buf[1024];
    snprintf(buf, 1024, "%s/%slog-%s-%s.txt",
             outputDirectory.c_str(),
             outputPrefix.c_str(),
             primary.c_str(),
             MigrationGraph::getPatternString(pattern).c_str());
    
    filenameGurobiLog = buf;
  }
  
  IlpPmhSolver solver(T,
                      primary,
                      pattern,
                      filenameGurobiLog,
                      forcedComigrations);
  return run(solver,
             T,
             primary,
             outputDirectory,
             outputPrefix,
             colorMap,
             pattern,
             nrThreads,
             outputILP,
             outputSearchGraph,
             timeLimit,
             bounds,
             forcedComigrations);
}

IntTriple IlpPmhSolver::run(IlpPmhSolver& solver,
                            const CloneTree& T,
                            const std::string& primary,
                            const std::string& outputDirectory,
                            const std::string& outputPrefix,
                            const StringToIntMap& colorMap,
                            MigrationGraph::Pattern pattern,
                            int nrThreads,
                            bool outputILP,
                            bool outputSearchGraph,
                            int timeLimit,
                            const IntTriple& bounds,
                            const StringPairList& forcedComigrations)
{
  char buf[1024];
  std::string filenameSearchGraph;
  
  solver.init(bounds);
  
  if (!outputDirectory.empty())
  {
    snprintf(buf, 1024, "%s/%sGG-%s-%s.dot",
             outputDirectory.c_str(),
             outputPrefix.c_str(),
             primary.c_str(),
             MigrationGraph::getPatternString(pattern).c_str());
    
    filenameSearchGraph = buf;
    
    if (outputILP)
    {
      snprintf(buf, 1024, "%s/%s%s-%s.lp",
               outputDirectory.c_str(),
               outputPrefix.c_str(),
               primary.c_str(),
               MigrationGraph::getPatternString(pattern).c_str());
      
      solver.exportModel(buf);
    }
  }
  
  lemon::Timer timer;
  bool solved = solver.solve(nrThreads, timeLimit);
  if (!solved)
  {
    std::cout << outputPrefix << "\t"
              << "(" << MigrationGraph::getAllowedPatternsString(pattern) << ")\t"
              << "-" << "\t"
              << "-" << "\t"
              << "-" << "\t"
              << "-" << "\t"
              << "-" << "\t"
              << "-" << "\t"
              << timer.realTime()
              << std::endl;

    IntTriple res;
    res.first = -1;
    res.second.first = -1;
    res.second.second = -1;
    
    return res;
  }

  MigrationGraph G(solver.T(), solver.lPlus());
  
  int mu = G.getNrMigrations();
  int gamma = G.getNrComigrations(solver.T(), solver.lPlus());
  int sigma = G.getNrSeedingSites();
  
  std::cout << outputPrefix << "\t"
            << "(" << MigrationGraph::getAllowedPatternsString(pattern) << ")\t"
            << mu << "\t"
            << gamma << "\t"
            << sigma << "\t"
            << G.getPatternString(G.getPattern(), G.isMonoclonal()) << "\t"
            << solver.LB() << "\t"
            << solver.UB() << "\t"
            << timer.realTime()
            << std::endl;
  
  if (!outputDirectory.empty())
  {
    snprintf(buf, 1024, "%s/%sT-%s-%s.dot",
             outputDirectory.c_str(),
             outputPrefix.c_str(),
             primary.c_str(),
             MigrationGraph::getPatternString(pattern).c_str());
    std::ofstream outT(buf);
    solver.writeCloneTree(outT, colorMap);
    outT.close();
    
    snprintf(buf, 1024, "%s/%sG-%s-%s.dot",
             outputDirectory.c_str(),
             outputPrefix.c_str(),
             primary.c_str(),
             MigrationGraph::getPatternString(pattern).c_str());
    std::ofstream outG(buf);
    G.writeDOT(outG, colorMap);
    outG.close();

    snprintf(buf, 1024, "%s/%sG-%s-%s.tree",
             outputDirectory.c_str(),
             outputPrefix.c_str(),
             primary.c_str(),
             MigrationGraph::getPatternString(pattern).c_str());
    std::ofstream outGraph(buf);
    G.write(outGraph);
    outGraph.close();
    
    snprintf(buf, 1024, "%s/%sT-%s-%s.tree",
             outputDirectory.c_str(),
             outputPrefix.c_str(),
             primary.c_str(),
             MigrationGraph::getPatternString(pattern).c_str());
    std::ofstream outTree(buf);
    solver.T().write(outTree);
    outTree.close();
    
    snprintf(buf, 1024, "%s/%sT-%s-%s.labeling",
             outputDirectory.c_str(),
             outputPrefix.c_str(),
             primary.c_str(),
             MigrationGraph::getPatternString(pattern).c_str());
    std::ofstream outLabeling(buf);
    solver.T().writeVertexLabeling(outLabeling, solver.lPlus());
    outLabeling.close();
    
    std::ofstream outGG(filenameSearchGraph.c_str());
    solver.writeSolutionGraphDOT(outGG, colorMap);
    outGG.close();
  }
  
  IntTriple res;
  res.first = mu;
  res.second.first = gamma;
  res.second.second = sigma;
  
  return res;
}

void IlpPmhSolver::initCallbacks()
{
}

void IlpPmhSolver::initConstraintsNonEdgesG()
{
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  
  GRBLinExpr sum;
  
  for (ArcIt a_ij(getTree()); a_ij != lemon::INVALID; ++a_ij)
  {
    const int ij = (*_pArcToIndex)[a_ij];
    Node v_i = getTree().source(a_ij);
    Node v_j = getTree().target(a_ij);
    const int i = (*_pNodeToIndex)[v_i];
    const int j = (*_pNodeToIndex)[v_j];
    
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      const int size_L_s = _L[s].size();
      for (int c = 0; c < size_L_s; ++c)
      {
        for (int t = 0; t < nrAnatomicalSites; ++t)
        {
          const int size_L_t = _L[t].size();
          for (int d = 0; d < size_L_t; ++d)
          {
            _model.addConstr(_xx[ij][s][c][t][d] <= _x[i][s][c]);
            _model.addConstr(_xx[ij][s][c][t][d] <= _x[j][t][d]);
            _model.addConstr(_xx[ij][s][c][t][d] >= _x[i][s][c] + _x[j][t][d] - 1);
          }
        }
      }
    }
  }
  
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    const int size_L_s = _L[s].size();
    for (int c = 0; c < size_L_s; ++c)
    {
      for (int t = 0; t < nrAnatomicalSites; ++t)
      {
        const int size_L_t = _L[t].size();
        for (int d = 0; d < size_L_t; ++d)
        {
          if (s == t && c == d) continue;
          
          for (ArcIt a_ij(getTree()); a_ij != lemon::INVALID; ++a_ij)
          {
            const int ij = (*_pArcToIndex)[a_ij];
            sum += _xx[ij][s][c][t][d];
          }
          _model.addConstr(_z[s][c][t][d] <= sum);
          sum.clear();
        }
      }
    }
  }
}

void IlpPmhSolver::initConstraintsG()
{
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();

  GRBLinExpr sum, sum2;
  
  // Every (t,d) != (P,0) has at most one incoming edge
  for (int t = 0; t < nrAnatomicalSites; ++t)
  {
    const int size_L_t = _L[t].size();
    for (int d = 0; d < size_L_t; ++d)
    {
      if (t == _primaryIndex && d == 0) continue;
      for (int s = 0; s < nrAnatomicalSites; ++s)
      {
        const int size_L_s = _L[s].size();
        for (int c = 0; c < size_L_s; ++c)
        {
          sum += _z[s][c][t][d];
        }
      }
      
      _model.addConstr(sum == _y[t][d]);
      
      sum.clear();
    }
  }
  
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    const int size_L_s = _L[s].size();
    for (int c = 0; c < size_L_s; ++c)
    {
      _model.addConstr(_z[s][c][s][c] == 0);
    }
  }
  
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    const int size_L_s = _L[s].size();
    for (int c = 0; c < size_L_s; ++c)
    {
      for (int t = 0; t < nrAnatomicalSites; ++t)
      {
        const int size_L_t = _L[t].size();
        for (int d = 0; d < size_L_t; ++d)
        {
          if (s == t && c == d) continue;
          _model.addConstr(_z[s][c][t][d] <= _y[s][c]);
          _model.addConstr(_z[s][c][t][d] <= _y[t][d]);
        }
      }
    }
  }
  
  // No edge from (s,c) to (s,d)
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    const int size_L_s = _L[s].size();
    for (int c = 0; c < size_L_s; ++c)
    {
      for (int d = 0; d < size_L_s; ++d)
      {
        sum += _z[s][c][s][d];
      }
      _model.addConstr(sum == 0);
      sum.clear();
    }
  }
  
  
  // (P,0). I don't think this is necessary
  /*for (int t = 0; t < nrAnatomicalSites; ++t)
  {
    const int size_L_t = _L[t].size();
    for (int d = 0; d < size_L_t; ++d)
    {
      sum += _z[_primaryIndex][0][t][d];
    }
  }
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    const int i = (*_pNodeToIndex)[v_i];
    _model.addConstr(sum >= _x[i][_primaryIndex][0]);
  }
  sum.clear();*/
}

bool IlpPmhSolver::solve(int nrThreads, int timeLimit)
{
  try
  {
    if (nrThreads > 0)
    {
      _model.getEnv().set(GRB_IntParam_Threads, nrThreads);
    }
    if (timeLimit > 0)
    {
      _model.getEnv().set(GRB_DoubleParam_TimeLimit, timeLimit);
    }
    
    _model.getEnv().set(GRB_IntParam_LogToConsole, 0);
    _model.optimize();
    
    int status = _model.get(GRB_IntAttr_Status);
    if (status == GRB_OPTIMAL || status == GRB_SUBOPTIMAL)
    {
      _LB = _model.get(GRB_DoubleAttr_ObjBound);
      _UB = _model.get(GRB_DoubleAttr_ObjVal);
      processSolution();
      return true;
    }
    else if (status == GRB_INF_OR_UNBD)
    {
      std::cerr << "Model is infeasible or unbounded" << std::endl;
      return false;
    }
    else if (status == GRB_INFEASIBLE)
    {
      //std::cerr << "Model is infeasible" << std::endl;
//      _model.computeIIS();
//      _model.write("/tmp/model_IIS.ilp");
      return false;
    }
    else if (status == GRB_UNBOUNDED)
    {
      std::cerr << "Model is unbounded" << std::endl;
      return false;
    }
    else if (status == GRB_TIME_LIMIT)
    {
      _LB = _model.get(GRB_DoubleAttr_ObjBound);
      _UB = _model.get(GRB_DoubleAttr_ObjVal);
      if (_UB < std::numeric_limits<double>::max() && _UB != NAN)
      {
        processSolution();
        return true;
      }
      return false;
    }
  }
  catch (const GRBException& e)
  {
    std::cerr << "Error code = " << e.getErrorCode() << std::endl;
    std::cerr << e.getMessage() << std::endl;
    
    return false;
  }
  catch (...)
  {
    return false;
  }
  
  return true;
}

void IlpPmhSolver::initVertexLabelingConstraints()
{
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  GRBLinExpr sum;
  
  // Every vertex is labeled by exactly one (s,c)
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    const int i = (*_pNodeToIndex)[v_i];
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      const int size_L_s = _L[s].size();
      for (int c = 0; c < size_L_s; ++c)
      {
        sum += _x[i][s][c];
      }
    }
    
    _model.addConstr(sum == 1);
    sum.clear();
  }
}

void IlpPmhSolver::initLeafVariables()
{
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  const int nrNodes = _indexToNode.size();
  
  char buf[1024];
  
  for (int i = 0; i < nrNodes; ++i)
  {
    Node v_i = _indexToNode[i];
    if (!isLeaf(v_i)) continue;
    
    _x[i] = VarMatrix(nrAnatomicalSites);
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      int size_L_s = _L[s].size();
      _x[i][s] = VarArray(size_L_s);
      for (int c = 0; c < size_L_s; ++c)
      {
        snprintf(buf, 1024, "x;%s;%s;%d",
                 getLabel(v_i).c_str(),
                 _indexToAnatomicalSite[s].c_str(), c);
        _x[i][s][c] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
      }
    }
  }
  
  _model.update();
}

void IlpPmhSolver::initVariables()
{
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  const int nrNodes = _indexToNode.size();
  const int nrArcs = _indexToArc.size();

  char buf[1024];
  
  _x = Var3Matrix(nrNodes);
  for (int i = 0; i < nrNodes; ++i)
  {
    Node v_i = _indexToNode[i];
    if (isLeaf(v_i)) continue;
    _x[i] = VarMatrix(nrAnatomicalSites);
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      int size_L_s = _L[s].size();
      _x[i][s] = VarArray(size_L_s);
      for (int c = 0; c < size_L_s; ++c)
      {
        snprintf(buf, 1024, "x;%s;%s;%d",
                 getLabel(v_i).c_str(),
                 _indexToAnatomicalSite[s].c_str(), c);
        _x[i][s][c] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
      }
    }
  }
  
  _xx = Var5Matrix(nrArcs);
  for (int ij = 0; ij < nrArcs; ++ij)
  {
    Arc a_ij = _indexToArc[ij];
    Node v_i = getTree().source(a_ij);
    Node v_j = getTree().target(a_ij);
    
    _xx[ij] = Var4Matrix(nrAnatomicalSites);
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      const int size_L_s = _L[s].size();
      _xx[ij][s] = Var3Matrix(size_L_s);
      for (int c = 0; c < size_L_s; ++c)
      {
        _xx[ij][s][c] = VarMatrix(nrAnatomicalSites);
        for (int t = 0; t < nrAnatomicalSites; ++t)
        {
          const int size_L_t = _L[t].size();
          _xx[ij][s][c][t] = VarArray(size_L_t);
          for (int d = 0; d < size_L_t; ++d)
          {
            snprintf(buf, 1024, "xx;(%s,%s);%s;%d;%s;%d",
                     getLabel(v_i).c_str(),
                     getLabel(v_j).c_str(),
                     _indexToAnatomicalSite[s].c_str(), c,
                     _indexToAnatomicalSite[t].c_str(), d);
            _xx[ij][s][c][t][d] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
          }
        }
      }
    }
  }
  
  _y = VarMatrix(nrAnatomicalSites);
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    int size_L_s = _L[s].size();
    _y[s] = VarArray(size_L_s);
    for (int c = 0; c < size_L_s; ++c)
    {
      snprintf(buf, 1024, "y;%s;%d", _indexToAnatomicalSite[s].c_str(), c);
      _y[s][c] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
    }
  }
  
  _maxNrEdgesInG = 0;
  _z = Var4Matrix(nrAnatomicalSites);
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    const int size_L_s = _L[s].size();
    _z[s] = Var3Matrix(size_L_s);
    for (int c = 0; c < size_L_s; ++c)
    {
      _z[s][c] = VarMatrix(nrAnatomicalSites);
      for (int t = 0; t < nrAnatomicalSites; ++t)
      {
        const int size_L_t = _L[t].size();
        _z[s][c][t] = VarArray(size_L_t);
        for (int d = 0; d < size_L_t; ++d)
        {
          snprintf(buf, 1024, "z;%s;%d;%s;%d",
                   _indexToAnatomicalSite[s].c_str(), c,
                   _indexToAnatomicalSite[t].c_str(), d);
          _z[s][c][t][d] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
          ++_maxNrEdgesInG;
        }
      }
    }
  }
  
  _w = Var4Matrix(nrAnatomicalSites);
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    _w[s] = Var3Matrix(nrAnatomicalSites);
    for (int t = 0; t < nrAnatomicalSites; ++t)
    {
      const int size_L_t = _L[t].size();
      _w[s][t] = VarMatrix(size_L_t);
      for (int d = 0; d < size_L_t; ++d)
      {
        _w[s][t][d] = VarArray(size_L_t);
        for (int e = 0; e < size_L_t; ++e)
        {
          snprintf(buf, 1024, "w;%s;%s;%d;%d",
                   _indexToAnatomicalSite[s].c_str(),
                   _indexToAnatomicalSite[t].c_str(),
                   d, e);
          _w[s][t][d][e] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
        }
      }
    }
  }
  
  _gamma = VarMatrix(nrAnatomicalSites);
  _sigma = VarArray(nrAnatomicalSites);
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    _gamma[s] = VarArray(nrAnatomicalSites);
    for (int t = 0; t < nrAnatomicalSites; ++t)
    {
      snprintf(buf, 1024, "gamma;%s;%s",
               _indexToAnatomicalSite[s].c_str(),
               _indexToAnatomicalSite[t].c_str());
      _gamma[s][t] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
    }
    
    snprintf(buf, 1024, "sigma;%s",
             _indexToAnatomicalSite[s].c_str());
    _sigma[s] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
  }
  
  _model.update();
}

void IlpPmhSolver::initLeafConstraints()
{
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  
  GRBLinExpr sum, sum2;

  // Ensure that each leaf v_i receives a unique label
  // satisfying the definition of a sublabeling
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    if (!isLeaf(v_i)) continue;
    
    const int i = (*_pNodeToIndex)[v_i];
    const std::string& sStr = getLeafAnatomicalSiteLabel(v_i);
    
    assert(_anatomicalSiteToIndex.count(sStr) == 1);
    const int s = _anatomicalSiteToIndex[sStr];
    
    const int size_L_s = _L[s].size();
    for (int c = 0; c < size_L_s; ++c)
    {
      sum += _x[i][s][c];
    }
    _model.addConstr(sum == 1);
    sum.clear();
    
    for (int t = 0; t < nrAnatomicalSites; ++t)
    {
      if (t == s) continue;
      
      const int size_L_t = _L[t].size();
      for (int c = 0; c < size_L_t; ++c)
      {
        sum += _x[i][t][c];
      }
    }
    _model.addConstr(sum == 0);
    sum.clear();
  }
  _model.update();
}

void IlpPmhSolver::initConstraints()
{
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  
  GRBLinExpr sum, sum2;
  
  // symmetry breaking constraints
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    const NodeSet& L_s = _L[s];
    const int size_L_s = _L[s].size();
    for (int c = 1; c < size_L_s; ++c)
    {
      for (Node v_i : L_s)
      {
        if (v_i == lemon::INVALID)
        {
          // dummy vertex
          if (c == 1)
          {
            sum += 1;
          }
        }
        else
        {
          const int i = (*_pNodeToIndex)[v_i];
          sum += _x[i][s][c - 1];
          sum2 += _x[i][s][c];
          //std::cout << _x[i][s][c - 1].get(GRB_StringAttr_VarName) << "    " << _x[i][s][c].get(GRB_StringAttr_VarName) << std::endl;
        }
      }
      _model.addConstr(sum2 <= sum);
      sum.clear();
      sum2.clear();
    }
  }
  
  // root vertex is labeled by (P,1)
  _model.addConstr(_x[(*_pNodeToIndex)[getRoot()]][_primaryIndex][0] == 1);
  
  // Intuitively, we want that x[i][s][c] = 1 if vertex v_i in X(s,c)
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    const NodeSet& L_s = _L[s];
    for (Node v_i : L_s)
    {
      for (Node v_j : L_s)
      {
        if (getTree().id(v_i) > getTree().id(v_j))
        {
          assert(v_i != lemon::INVALID);
          
          Node v_k;
          if (v_j == lemon::INVALID)
          {
            v_k = getRoot();
          }
          else
          {
            NodeSet set_ij;
            set_ij.insert(v_i);
            set_ij.insert(v_j);
            
            v_k = getLCA(set_ij);
            assert(!isLeaf(v_k));
          }

          initConstraintsLCA(s, v_k, v_i, v_j);
          
          Node v_kk = getParent(v_i);
          while (v_kk != v_k)
          {
            initConstraintsLCA(s, v_kk, v_i, v_j);
            v_kk = getParent(v_kk);
          }
          
          if (v_j != lemon::INVALID)
          {
            v_kk = getParent(v_j);
            while (v_kk != v_k)
            {
              initConstraintsLCA(s, v_kk, v_i, v_j);
              v_kk = getParent(v_kk);
            }
          }
        }
      }
    }
  }
  
  // if _x[i][s][c] = 1 then at least one of the children v_j
  // must have _x[j][s][c] = 1
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    const int size_L_s = _L[s].size();
    for (int c = 0; c < size_L_s; ++c)
    {
      for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
      {
        if (isLeaf(v_i)) continue;
        if (s == _primaryIndex && c == 0) continue;
        
        const int i = (*_pNodeToIndex)[v_i];
        for (OutArcIt a_ij(getTree(), v_i); a_ij != lemon::INVALID; ++a_ij)
        {
          Node v_j = getTree().target(a_ij);
          const int j = (*_pNodeToIndex)[v_j];
          sum += _x[j][s][c];
        }
        
        _model.addConstr(_x[i][s][c] <= sum);
        sum.clear();
      }
    }
  }

  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    if (isLeaf(v_i)) continue;
    if (getRoot() == v_i) continue;
    
    const int i = (*_pNodeToIndex)[v_i];

    Node v_k = getParent(v_i);
    const int k = (*_pNodeToIndex)[v_k];

    _model.addConstr(_x[i][_primaryIndex][0] <= _x[k][_primaryIndex][0]);
  }
  
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    const int size_L_s = _L[s].size();
    for (int c = 0; c < size_L_s; ++c)
    {
      for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
      {
        const int i = (*_pNodeToIndex)[v_i];
        _model.addConstr(_y[s][c] >= _x[i][s][c]);
        if (isLeaf(v_i))
        {
          sum += _x[i][s][c];
        }
      }
      _model.addConstr(_y[s][c] <= sum);
      sum.clear();
    }
  }
  
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    for (int t = 0; t < nrAnatomicalSites; ++t)
    {
      const int size_L_t = _L[t].size();
      for (int d = 0; d < size_L_t; ++d)
      {
        _model.addConstr(_w[s][t][d][d] == 0);
      }
    }
  }

  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    const int size_L_s = _L[s].size();
    for (int t = 0; t < nrAnatomicalSites; ++t)
    {
      const int size_L_t = _L[t].size();
      for (int d = 0; d < size_L_t; ++d)
      {
        for (int e = 0; e < size_L_t; ++e)
        {
          if (d == e) continue;
          
          for (int c = 0; c < size_L_s; ++c)
          {
            for (int cc = 0; cc < size_L_s; ++cc)
            {
              _model.addConstr(_w[s][t][d][e] >= _z[s][c][t][d] + _z[s][cc][t][e] - 1);
            }
          }
          
          for (int c = 0; c < size_L_s; ++c)
          {
            sum += _z[s][c][t][d];
            sum2 += _z[s][c][t][e];
          }
          _model.addConstr(_w[s][t][d][e] <= sum);
          _model.addConstr(_w[s][t][d][e] <= sum2);
          sum.clear();
          sum2.clear();
        }
      }
    }
  }
  
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    for (int t = 0; t < nrAnatomicalSites; ++t)
    {
      if (s == t)
      {
        _model.addConstr(_gamma[s][t] == 0);
      }
      else
      {
        const int size_L_s = _L[s].size();
        const int size_L_t = _L[t].size();
        for (int c = 0; c < size_L_s; ++c)
        {
          for (int d = 0; d < size_L_t; ++d)
          {
            _model.addConstr(_gamma[s][t] >= _z[s][c][t][d]);
            sum += _z[s][c][t][d];
          }
        }
        
        _model.addConstr(_gamma[s][t] <= sum);
        sum.clear();
      }
    }
  }
  
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    for (int t = 0; t < nrAnatomicalSites; ++t)
    {
      if (s == t) continue;
      _model.addConstr(_sigma[s] >= _gamma[s][t]);
      sum += _gamma[s][t];
    }
    _model.addConstr(_sigma[s] <= sum);
    sum.clear();
  }

  if (_pattern == MigrationGraph::S
      || _pattern == MigrationGraph::PS)
  {
    initSingleSourceSeedingConstraints();
  }
  
  if (_pattern == MigrationGraph::PS)
  {
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      if (s == _primaryIndex) continue;
      
      const int size_L_s = _L[s].size();
      for (int c = 0; c < size_L_s; ++c)
      {
        for (Node v_i : _L[s])
        {
          const int i = (*_pNodeToIndex)[v_i];
          _model.addConstr(_z[_primaryIndex][0][s][c] >= _x[i][s][c]);
        }
      }
    }
  }
  
  if (_pattern == MigrationGraph::M)
  {
    initMultiSourceSeedingConstraints();
  }
  
  _model.update();
}

void IlpPmhSolver::initSingleSourceSeedingConstraints()
{
  const int nrAnatomicalSites = _indexToAnatomicalSite.size();
  GRBLinExpr sum;
  
  for (Node v_i : _L[_primaryIndex])
  {
    if (v_i == lemon::INVALID) continue;
    
    const int i = (*_pNodeToIndex)[v_i];
    _model.addConstr(_x[i][_primaryIndex][0] == 1);
  }

  for (int t = 0; t < nrAnatomicalSites; ++t)
  {
    if (t == _primaryIndex) continue;
    
    const int size_L_t = _L[t].size();
    for (int d = 0; d < size_L_t; ++d)
    {
      for (int e = 0; e < size_L_t; ++e)
      {
        if (d == e) continue;
        for (int s = 0; s < nrAnatomicalSites; ++s)
        {
          sum += _w[s][t][d][e];
        }
        
        _model.addConstr(sum >= _y[t][d] + _y[t][e] - 1);
        _model.addConstr(sum <= _y[t][d]);
        _model.addConstr(sum <= _y[t][e]);
        sum.clear();
      }
    }
  }
}

void IlpPmhSolver::initMultiSourceSeedingConstraints()
{
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    const int i = (*_pNodeToIndex)[v_i];
    for (NodeIt v_j(getTree()); v_j != lemon::INVALID; ++v_j)
    {
      const int j = (*_pNodeToIndex)[v_j];
      if (v_i == v_j) continue;
      if (isAncestor(v_i, v_j))
      {
        for (int s = 0; s < nrAnatomicalSites; ++s)
        {
          const int size_L_s = _L[s].size();
          for (int c = 0; c < size_L_s; ++c)
          {
            for (int d = 0; d < size_L_s; ++d)
            {
              if (c == d) continue;
              _model.addConstr(_x[i][s][c] + _x[j][s][d] <= 1);
            }
          }
        }
      }
    }
  }
}

void IlpPmhSolver::initConstraintsLCA(int s,
                                      Node v_k,
                                      Node v_i,
                                      Node v_j)
{
  assert(0 <= s && s < _anatomicalSiteToIndex.size());
  assert(v_i != lemon::INVALID);
  assert(v_k != lemon::INVALID);
  
  const int i = (*_pNodeToIndex)[v_i];
  const int k = (*_pNodeToIndex)[v_k];
  const int size_L_s = _L[s].size();
  
  for (int c = 0; c < size_L_s; ++c)
  {
    if (v_j == lemon::INVALID)
    {
      assert(s == _primaryIndex);
      if (c == 0)
      {
        _model.addConstr(_x[k][s][c] >= _x[i][s][c]);
      }
      else
      {
        // TODO: can be removed
        _model.addConstr(_x[k][s][c] >= _x[i][s][c] - 1);
      }
    }
    else
    {
      int j = (*_pNodeToIndex)[v_j];
      _model.addConstr(_x[k][s][c] >= _x[i][s][c] + _x[j][s][c] - 1);
    }
  }
}

GRBLinExpr IlpPmhSolver::initObjective(const IntTriple& bounds)
{
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  
  // migration number
  GRBLinExpr migrationNumber;
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    const int size_L_s = _L[s].size();
    for (int c = 0; c < size_L_s; ++c)
    {
      if (s == _primaryIndex && c == 0) continue;
      migrationNumber += _y[s][c];
    }
  }
  
  if (bounds.first != -1)
  {
    _model.addConstr(migrationNumber <= bounds.first);
  }

  // comigration number
  GRBLinExpr comigrationNumber;
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    for (int t = 0; t < nrAnatomicalSites; ++t)
    {
      comigrationNumber += _gamma[s][t];
    }
  }
  
  if (bounds.second.first != -1)
  {
    _model.addConstr(comigrationNumber <= bounds.second.first);
  }
  
  GRBLinExpr seedingSiteNumber;
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    seedingSiteNumber += _sigma[s];
  }
  
  if (bounds.second.second != -1)
  {
    _model.addConstr(seedingSiteNumber <= bounds.second.second);
  }
  
  GRBLinExpr obj;
  obj += migrationNumber
      + (1. / (nrAnatomicalSites * nrAnatomicalSites)) * comigrationNumber
      + (1. / (nrAnatomicalSites * nrAnatomicalSites)) * (1. / (nrAnatomicalSites + 1)) * seedingSiteNumber;

  _model.setObjective(obj, GRB_MINIMIZE);
  _model.update();
  
  return obj;
}

void IlpPmhSolver::initForcedComigrations()
{
  for (const StringPair& st : _forcedComigrations)
  {
    assert(_anatomicalSiteToIndex.count(st.first) == 1);
    assert(_anatomicalSiteToIndex.count(st.second) == 1);
    
    int s = _anatomicalSiteToIndex[st.first];
    int t = _anatomicalSiteToIndex[st.second];
    
    _model.addConstr(_gamma[s][t] == 1);
  }
}

void IlpPmhSolver::exportModel(const std::string& filename)
{
  _model.write(filename);
}

void IlpPmhSolver::constructGraph()
{
  const int nrAnatomicalSites = _indexToAnatomicalSite.size();
  
  // Initialize _subLabelToNodeG
  _subLabelToNodeG = NodeMatrix(nrAnatomicalSites);
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    _subLabelToNodeG[s] = NodeVector(_L[s].size(), lemon::INVALID);
  }
  
  // Initialize edges and vertices
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    const int size_L_s = _L[s].size();
    for (int c = 0; c < size_L_s; ++c)
    {
      for (int t = 0; t < nrAnatomicalSites; ++t)
      {
        const int size_L_t = _L[t].size();
        for (int d = 0; d < size_L_t; ++d)
        {
          if (_z[s][c][t][d].get(GRB_DoubleAttr_X) >= 0.4)
          {
            Node v_sc = _subLabelToNodeG[s][c];
            if (v_sc == lemon::INVALID)
            {
              v_sc = _G.addNode();
              if (s == _primaryIndex && c == 0)
              {
                _rootG = v_sc;
              }
              _nodeGToSubLabel[v_sc] = std::make_pair(s, c);
              _subLabelToNodeG[s][c] = v_sc;
            }
            
            Node v_td = _subLabelToNodeG[t][d];
            if (v_td == lemon::INVALID)
            {
              v_td = _G.addNode();
              _nodeGToSubLabel[v_td] = std::make_pair(t, d);
              _subLabelToNodeG[t][d] = v_td;
            }
            
//            std::cout << "Edge from (" << _indexToAnatomicalSite[s] << "," << c << ") to (" << _indexToAnatomicalSite[t] << "," << d << ")" << std::endl;
            
            _G.addArc(v_sc, v_td);
          }
        }
      }
    }
  }
  
  assert(lemon::dag(_G));
}

void IlpPmhSolver::writeSolutionGraphDOT(std::ostream& out,
                                         const StringToIntMap& colorMap) const
{
  out << "digraph G {" << std::endl;
  
  for (NodeIt v_sc(_G); v_sc != lemon::INVALID; ++v_sc)
  {
    const IntPair& sc = _nodeGToSubLabel[v_sc];
    const std::string& sStr = _indexToAnatomicalSite[sc.first];
    
    out << "\t" << _G.id(v_sc) << " [penwidth=3,colorscheme=set19,color="
        << colorMap.find(sStr)->second << ",label=\"("
        << sStr << "," << sc.second << ")\"]" << std::endl;
  }
  
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    Node v_sc = _G.source(a);
    Node v_td = _G.target(a);
    
    const IntPair& sc = _nodeGToSubLabel[v_sc];
    const IntPair& td = _nodeGToSubLabel[v_td];
    
    const std::string& sStr = _indexToAnatomicalSite[sc.first];
    const std::string& tStr = _indexToAnatomicalSite[td.first];
    
    out << "\t" << _G.id(v_sc) << " -> " << _G.id(v_td);
    if (sc.first == td.first)
    {
      out << " [penwidth=3,colorscheme=set19,color=" << colorMap.find(sStr)->second << "]";
    }
    else
    {
      out << " [penwidth=3,colorscheme=set19,color=\"" << colorMap.find(sStr)->second << ";0.5:" << colorMap.find(tStr)->second << "\"]";
    }
    
    out << std::endl;
  }
  
  out << "}";
}

void IlpPmhSolver::processSolution()
{
  const int nrNodes = _indexToNode.size();
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  
  _pLPlus = new StringNodeMap(getTree());

  for (int i = 0; i < nrNodes; ++i)
  {
    Node v_i = _indexToNode[i];
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      const std::string& sStr = _indexToAnatomicalSite[s];

      int size_L_s = _L[s].size();
      for (int c = 0; c < size_L_s; ++c)
      {
        if (_x[i][s][c].get(GRB_DoubleAttr_X) >= 0.4)
        {
          //std::cout << _x[i][s][c].get(GRB_StringAttr_VarName)
          //          << " = " << _x[i][s][c].get(GRB_DoubleAttr_X) << std::endl;
          _pLPlus->set(v_i, sStr);
        }
      }
    }
  }
  
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    // march up until you find a labeled ancestor
    Node v_j = v_i;
    while (v_j != lemon::INVALID && (*_pLPlus)[v_j].empty())
    {
      v_j = getParent(v_j);
    }
    assert(v_j != lemon::INVALID);
    
    if (v_j != v_i)
    {
      _pLPlus->set(v_i, (*_pLPlus)[v_j]);
    }
  }
  
//
//  std::cout << std::endl;
//  
//  for (int s = 0; s < nrAnatomicalSites; ++s)
//  {
//    int size_L_s = _L[s].size();
//    for (int c = 0; c < size_L_s; ++c)
//    {
//      if (_y[s][c].get(GRB_DoubleAttr_X) >= 0.4)
//      {
//        std::cout << _y[s][c].get(GRB_StringAttr_VarName)
//                  << " = " << _y[s][c].get(GRB_DoubleAttr_X) << std::endl;
//      }
//    }
//  }
//
//  std::cout << std::endl;
//  
//  for (int s = 0; s < nrAnatomicalSites; ++s)
//  {
//    const int size_L_s = _L[s].size();
//    for (int c = 0; c < size_L_s; ++c)
//    {
//      for (int t = 0; t < nrAnatomicalSites; ++t)
//      {
//        const int size_L_t = _L[t].size();
//        for (int d = 0; d < size_L_t; ++d)
//        {
//          if (_z[s][c][t][d].get(GRB_DoubleAttr_X) >= 0.4)
//          {
//            std::cout << _z[s][c][t][d].get(GRB_StringAttr_VarName)
//                      << " = " << _z[s][c][t][d].get(GRB_DoubleAttr_X) << std::endl;
//          }
//        }
//      }
//    }
//  }
  
  constructGraph();
}

void IlpPmhSolver::writeCloneTree(std::ostream& out,
                                  const StringToIntMap& colorMap) const
{
  T().writeDOT(out, lPlus(), colorMap);
}
