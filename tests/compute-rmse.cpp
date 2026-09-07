/**
 * @file compute-rmse.cpp
 * @brief Tool CLI per il calcolo della Root Mean Square Error (RMSE) e Spread
 *        tra due popolazioni (es. baseline sequenziale vs parallela).
 */

#include "fmo/metrics/metrics.hpp"
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cerr << "Uso: " << argv[0]
              << " <pop_seq.csv> <pop_par.csv> [variant] [mode] [workers] [rep]\n";
    return 1;
  }

  const std::string file_seq = argv[1];
  const std::string file_par = argv[2];
  const std::string variant = (argc > 3) ? argv[3] : "unknown";
  const std::string mode = (argc > 4) ? argv[4] : "unknown";
  const int workers = (argc > 5) ? std::atoi(argv[5]) : 1;
  const int rep = (argc > 6) ? std::atoi(argv[6]) : 1;

  try {
    Population pop_seq = loadPopulation(file_seq);
    Population pop_par = loadPopulation(file_par);

    if (pop_seq.size() == 0 || pop_par.size() == 0) {
      std::cerr << "Errore: una o entrambe le popolazioni sono vuote.\n";
      return 1;
    }
    if (pop_seq.size() != pop_par.size()) {
      std::cerr << "Errore: dimensione delle popolazioni non coerente ("
                << pop_seq.size() << " vs " << pop_par.size() << ")\n";
      return 1;
    }

    double spread_seq = computeSpreadMetric(pop_seq);
    double spread_par = computeSpreadMetric(pop_par);
    PopulationDiscrepancy disc = computePopulationRMSE(pop_seq, pop_par);

    // Formato machine-readable (per grep/awk nello script bash e python)
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "RMSE variant=" << variant
              << " mode=" << mode
              << " workers=" << workers
              << " rep=" << rep
              << " rmse_ptv=" << disc.rmse_ptv
              << " rmse_rectum=" << disc.rmse_rectum
              << " rmse_bladder=" << disc.rmse_bladder
              << " rmse_total_fitness=" << disc.rmse_total_fitness
              << " rmse_genes=" << disc.rmse_genes
              << " spread_seq=" << spread_seq
              << " spread_par=" << spread_par
              << std::endl;

    return 0;
  } catch (const std::exception &e) {
    std::cerr << "Errore durante il calcolo delle metriche: " << e.what()
              << std::endl;
    return 1;
  }
}
