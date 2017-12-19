/*
 * simulationmain.cpp
 *
 *  Created on: 25-aug-2017
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "simulation/simulation.h"
#include <lemon/arg_parser.h>
#include <fstream>

void output(const Simulation& simulation,
            const std::string& outputDirectory,
            const StringToIntMap& colorMap,
            int seed)
{
  char buf[1024];
  
  const CloneTree& T = simulation.getCloneTree();
  const CloneTree& TT = simulation.getObservedCloneTree();
  const StringNodeMap& lPlus = simulation.getObservedVertexLabeling();
//  MigrationGraph G = cellTree.getMigrationGraph();
  MigrationGraph GG = simulation.getObservedMigrationGraph();
  
  snprintf(buf, 1024, "%s/T_all_seed%d.dot", outputDirectory.c_str(), seed);
  std::ofstream outT_dot(buf);
  T.writeDOT(outT_dot, colorMap, simulation.getAnatomicalSiteProportions());
  outT_dot.close();
  
  snprintf(buf, 1024, "%s/clustering_observed_seed%d.txt", outputDirectory.c_str(), seed);
  std::ofstream outClustering(buf);
  simulation.writeObservedClustering(outClustering);
  outClustering.close();
  
  snprintf(buf, 1024, "%s/T_seed%d.dot", outputDirectory.c_str(), seed);
  std::ofstream outTT_dot(buf);
  TT.writeDOT(outTT_dot, lPlus, colorMap, simulation.getObservedSampleProportions());
  outTT_dot.close();
  
  snprintf(buf, 1024, "%s/T_seed%d.tree", outputDirectory.c_str(), seed);
  std::ofstream outTT_tree(buf);
  TT.write(outTT_tree);
  outTT_tree.close();
  
  snprintf(buf, 1024, "%s/T_seed%d.labeling", outputDirectory.c_str(), seed);
  std::ofstream outTT_labeling(buf);
  TT.writeLeafLabeling(outTT_labeling);
  outTT_labeling.close();
  
  snprintf(buf, 1024, "%s/T_seed%d.vertex.labeling", outputDirectory.c_str(), seed);
  std::ofstream outTT_vertex_labeling(buf);
  TT.writeVertexLabeling(outTT_vertex_labeling, lPlus);
  outTT_vertex_labeling.close();
  
  snprintf(buf, 1024, "%s/drivers_seed%d.txt", outputDirectory.c_str(), seed);
  std::ofstream outDrivers(buf);
  simulation.writeDrivers(outDrivers);
  outDrivers.close();
  
//  snprintf(buf, 1024, "%s/G_all_seed%d.dot", outputDirectory.c_str(), seed);
//  std::ofstream outG_dot(buf);
//  G.writeDOT(outG_dot, colorMap);
//  outG_dot.close();
  
  snprintf(buf, 1024, "%s/G_seed%d.dot", outputDirectory.c_str(), seed);
  std::ofstream outGG_dot(buf);
  GG.writeDOT(outGG_dot, colorMap);
  outGG_dot.close();
  
  snprintf(buf, 1024, "%s/G_seed%d.tree", outputDirectory.c_str(), seed);
  std::ofstream outGG_tree(buf);
  GG.write(outGG_tree);
  outGG_tree.close();
  
  snprintf(buf, 1024, "%s/reads_seed%d.tsv", outputDirectory.c_str(), seed);
  std::ofstream out_reads(buf);
  simulation.writeReadCounts(out_reads);
  out_reads.close();
}

bool success(const MigrationGraph& G,
             Simulation::Pattern pattern)
{
  switch (pattern)
  {
    case Simulation::PATTERN_mS:
      return G.getNrMigrations() == G.getNrAnatomicalSites() - 1;
    case Simulation::PATTERN_S:
      return G.isPolyclonal() && G.isSingleSourceSeeded();
    case Simulation::PATTERN_M:
      return G.isPolyclonal() && G.getPattern() == MigrationGraph::M;
    case Simulation::PATTERN_R:
      return G.isPolyclonal() && G.hasReseeding();
  }
}

int main(int argc, char** argv)
{
  int seed = 0;
  double K = 5e4;
  double mutFreqThreshold = 0.05;
  double migrationRate = 1e-6;
  double mutationRate = 0.1;
  double driverProb = 1e-7;
  int maxNrAnatomicalSites = 8;
  std::string filenameColorMap;
  std::string outputDirectory = ".";
  int pattern = 0;
  int N = -1;
  int coverage = 200;
  int nrSamplesPerAnatomicalSite = 2;
  int nrSamplesPrimary = 2;
  double seqErrorRate = 0;
  double purity = 1;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("s", "Random number generator seed (default: 0)", seed, false)
    .refOption("c", "Color map file", filenameColorMap, true)
    .refOption("C", "Target coverage", coverage, true)
    .refOption("D", "Driver probability (default: 1e-7)", driverProb)
    .refOption("E", "Per base sequencing error rate (default: 0)", seqErrorRate)
    .refOption("P", "Purity (default: 1)", purity)
    .refOption("K", "Carrying capacity (default: 5e4)", K)
    .refOption("k", "Number of samples per anatomical site (default: 2)", nrSamplesPerAnatomicalSite)
    .refOption("kP", "Number of samples for the primary tumor (default: 2)", nrSamplesPrimary)
    .refOption("f", "Mutation frequency threshold (default: 0.05)", mutFreqThreshold)
    .refOption("mig", "Migration rate (default: 1e-6)", migrationRate)
    .refOption("mut", "Mutation rate (default: 0.1)", mutationRate)
    .refOption("N", "Number of successful simulations (default: -1)", N)
    .refOption("m", "Maximum number of detectable anatomical sites (default: 8)", maxNrAnatomicalSites)
    .refOption("p", "Allowed migration patterns:\n"\
               "       0 : mS (default)\n"\
               "       1 : mS, S\n" \
               "       2 : mS, S, M\n" \
               "       3 : mS, S, M, R", pattern)
    .refOption("o", "Output directory (default: '.')", outputDirectory);
  ap.parse();
  
  // Read color map
  StringToIntMap colorMap;
  std::ifstream inColorMap(filenameColorMap.c_str());
  if (!inColorMap.good())
  {
    std::cerr << "Could not open '" << filenameColorMap << "' for reading" << std::endl;
    return 1;
  }
  if (!BaseTree::readColorMap(inColorMap, colorMap))
  {
    return 1;
  }
  
  Simulation::Pattern migrationPattern = static_cast<Simulation::Pattern>(pattern);
  
  if (N == -1)
  {
    g_rng = std::mt19937(seed);
    Simulation simulation(K,
                          migrationRate,
                          mutationRate,
                          driverProb,
                          mutFreqThreshold,
                          maxNrAnatomicalSites,
                          nrSamplesPerAnatomicalSite,
                          nrSamplesPrimary,
                          coverage,
                          migrationPattern,
                          seqErrorRate,
                          purity);
    if (simulation.simulate(true))
    {
      output(simulation, outputDirectory, colorMap, seed);
    }
  }
  else
  {
    int n = 0;
    while (n < N)
    {
      std::cout << "Seed: " << seed << "... " << std::flush;
      g_rng = std::mt19937(seed);
      Simulation simulation(K,
                            migrationRate,
                            mutationRate,
                            driverProb,
                            mutFreqThreshold,
                            maxNrAnatomicalSites,
                            nrSamplesPerAnatomicalSite,
                            nrSamplesPrimary,
                            coverage,
                            migrationPattern,
                            seqErrorRate,
                            purity);
      if (simulation.simulate(true))
      {
        MigrationGraph GG = simulation.getObservedMigrationGraph();
        if (success(GG, migrationPattern))
        {
          ++n;
          output(simulation, outputDirectory, colorMap, seed);
          std::cout << "SUCCESS" << std::endl;
        }
        else
        {
          std::cout << "FAILURE" << std::endl;
        }
      }
      else
      {
        std::cout << "FAILURE" << std::endl;
      }
      ++seed;
    }
  }
  
  return 0;
}
