/**
 * @file nsga2-seq.cpp
 * @brief Implementazione dell'algoritmo NSGA-II (Deb et al., 2002) in modalità
 * sequenziale per la risoluzione del problema FMO.
 * @details L'algoritmo NSGA-II (Non-dominated Sorting Genetic Algorithm II) è
 * un algoritmo genetico multi-obiettivo che utilizza la selezione a torneo
 * binario, il crossover SBX (Simulated Binary Crossover) e la mutazione
 * polinomiale per generare una popolazione di soluzioni ottimali. Questa
 * implementazione sequenziale esegue tutte le operazioni in un singolo thread,
 * senza parallelizzazione.
 */

#include "fmo/core/individual.hpp"
#include "fmo/evaluation/evaluator.hpp"
#include "fmo/nsga2/gen-op.hpp"
#include "fmo/nsga2/nsga2-steps.hpp"
#include "fmo/nsga2/offspring.hpp"
#include "fmo/utilities/hpc_helpers.hpp"
#include "fmo/utilities/utimer.hpp"

#include "fmo/preprocessing/manager.hpp"

#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>
#include <nsga-utils.hpp>
#define GANTRIES                                                               \
  {0, 72, 144, 216, 288}                                                       \
  // Angoli di gantry da caricare per la matrice globale D
#define PATH "../data/prostate" // Percorso base per i file di input

FMODataManager getDataFromPathAndAngles(const std::string &base_dir,
                                        const std::vector<int> &gantry_angles) {
  FMODataManager manager;
  manager.loadAllData(base_dir, gantry_angles);
  return manager;
}

void nsga2seq(Population &pop, int num_generations, int population_size, double eta_c, double eta_m,
              Evaluator &evaluator, std::mt19937 &rng,
              std::map<std::string, std::vector<long>> &samples) {
  auto measure = [&](const std::string &phase, auto operation) {
    long elapsed_ms = 0;
    {
      utimer timer(phase, &elapsed_ms, true);
      operation();
    }
    samples[phase].push_back(elapsed_ms);
  };

  // 1. Valutazione iniziale di P_0
  measure("initial_evaluation", [&] { evaluatePopulation(pop, evaluator); });

  // Classificazione iniziale di P_0
  std::vector<std::vector<int>> fronts;
  measure("initial_sorting", [&] { fronts = sortPopulation(pop); });

  measure("initial_crowding", [&] {
    assignPopulationCrowding(pop, fronts);
  });

  // 2. Loop Generazionale
  for (int gen = 0; gen < num_generations; ++gen) {
    long generation_ms = 0;
    {
      utimer generation_timer("generation_total", &generation_ms, true);

      // A. Generazione discendenza Q_t (taglia N) tramite Torneo, SBX e Mutazione
      Population offspring(pop.size());
      measure("offspring_generation", [&] {
        offspring = generatePopulationOffspring(pop, eta_c, eta_m, rng);
      });

      // B. Valutazione della discendenza Q_t (calcolo delle fitness)
      measure("population_evaluation", [&] {
        evaluatePopulation(offspring, evaluator);
      });

      // C. Fusione R_t = P_t U Q_t (taglia 2N)
      Population combined_pop;
      measure("population_merge", [&] {
        combined_pop = mergePopulations(pop, offspring);
      });

      // D. Non-dominated sorting ed estrazione dei fronti su R_t
      std::vector<std::vector<int>> combined_fronts;
      measure("population_sorting", [&] {
        combined_fronts = sortPopulation(combined_pop);
      });

      measure("population_crowding", [&] {
        assignPopulationCrowding(combined_pop, combined_fronts);
      });

      // E. Elitismo e troncamento: R_t -> P_{t+1} (taglia N)
      measure("population_truncation", [&] {
        pop = truncatePopulationByFronts(combined_pop, combined_fronts,
                                         population_size);
      });
    }
    if (gen > 0) {
      samples["generation_total"].push_back(generation_ms);
    }

  }
}

int main(int argc, char *argv[]) {
  constexpr int population_size = 100;
  constexpr int num_generations = 50;
  std::map<std::string, std::vector<long>> samples;
  std::cout << "RUN variant=seq mode=sequential workers=1 generations="
            << num_generations << " population=" << population_size
            << std::endl;

  long data_loading_ms = 0;
  FMODataManager manager;
  {
    utimer timer("data_loading", &data_loading_ms, true);
    manager = getDataFromPathAndAngles(PATH, GANTRIES);
  }
  samples["data_loading"].push_back(data_loading_ms);
  manager.printSummary();

  // generazione della popolazione iniziale
  long initial_population_ms = 0;
  std::mt19937 rng(42); // Inizializza il generatore di numeri casuali con un
                        // seed fisso per la riproducibilità
  Population start;
  {
    utimer timer("initial_population", &initial_population_ms, true);
    start = generateRandomPopulation(population_size,
                                     manager.getTotalBeamlets(), rng);
  }
  samples["initial_population"].push_back(initial_population_ms);

  // Esecuzione dell'algoritmo NSGA-II per un numero prefissato di generazioni
  long nsga2_total_ms = 0;
  Evaluator evaluator(manager.getData());
  {
    utimer timer("nsga2_total", &nsga2_total_ms, true);
    nsga2seq(start, num_generations, population_size, 20.0, 20.0, evaluator,
             rng, samples);
  }
  samples["nsga2_total"].push_back(nsga2_total_ms);

  double spread = computeSpreadMetric(start);
  std::cout << "QUALITY variant=seq mode=sequential spread_metric=" 
            << spread << std::endl;

  for (const auto &[phase, values] : samples) {
    std::vector<long> filtered(values.begin(), values.end());
    if (phase != "data_loading" && phase != "initial_population" &&
        phase != "nsga2_total" && filtered.size() > 1) {
      filtered.erase(filtered.begin());
    }
    if (filtered.empty()) {
      continue;
    }
    std::sort(filtered.begin(), filtered.end());
    const double mean = static_cast<double>(
        std::accumulate(filtered.begin(), filtered.end(), 0L)) /
        filtered.size();
    const double median = filtered.size() % 2 == 0
        ? (filtered[filtered.size() / 2 - 1] +
           filtered[filtered.size() / 2]) / 2.0
        : filtered[filtered.size() / 2];
    std::cout << "TIMING variant=seq mode=sequential phase=" << phase
              << " samples=" << filtered.size() << " mean_ms=" << mean
              << " median_ms=" << median << std::endl;
  }

  return 0;
}
