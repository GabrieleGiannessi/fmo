/**
 * @file nsga2-omp.cpp
 * @brief Implementazione dell'algoritmo NSGA-II (Deb et al., 2002) con
 * parallelizzazione OpenMP per la generazione della prole.
 * @details
 */

#include "fmo/core/individual.hpp"
#include "fmo/evaluation/evaluator.hpp"
#include "fmo/nsga2/gen-op.hpp"
#include "fmo/nsga2/nsga2-steps-omp.hpp"
#include "fmo/nsga2/nsga2-steps.hpp"
#include "fmo/nsga2/offspring.hpp"
#include "fmo/preprocessing/manager.hpp"
#include "fmo/utilities/hpc_helpers.hpp"

#include <cstdlib>
#include <iostream>
#include <vector>

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

void nsga2Omp(Population &pop, int num_generations, int population_size,
              double crossover_probability, double mutation_probability,
              double eta_c, double eta_m, Evaluator &evaluator,
              std::mt19937 &rng, int nw) {

  //   std::cout << "Esecuzione dell'algoritmo NSGA-II per " << num_generations
  //             << " generazioni..." << std::endl;

  //   std::cout << "Valutazione iniziale degli individui" << std::endl;

  // 1. Valutazione iniziale di P_0 usando il metodo parallelizzato tramite la
  // libreria OpenMP
  evaluatePopulationOmp(pop, evaluator, nw);

  std::cout << "Classificazione degli individui" << std::endl;

  // Classificazione iniziale di P_0
  auto fronts = sortPopulation(pop);
  assignPopulationCrowding(pop, fronts);

  // 2. Loop Generazionale
  for (int gen = 0; gen < num_generations; ++gen) {

    std::cout << "Generazione " << gen << std::endl;
    // A. Generazione discendenza Q_t (taglia N) tramite Torneo, SBX e Mutazione
    Population offspring = generatePopulationOffspring(
        pop, crossover_probability, mutation_probability, eta_c, eta_m, rng);

    // B. Valutazione della discendenza Q_t (calcolo delle fitness)
    evaluatePopulationOmp(offspring, evaluator, nw);

    // C. Fusione R_t = P_t U Q_t (taglia 2N)
    Population combined_pop = mergePopulations(pop, offspring);

    // D. Non-dominated sorting ed estrazione dei fronti su R_t
    auto combined_fronts = sortPopulation(combined_pop);
    assignPopulationCrowding(combined_pop, combined_fronts);

    // E. Elitismo e troncamento: R_t -> P_{t+1} (taglia N)
    pop = truncatePopulationByFronts(combined_pop, combined_fronts,
                                     population_size);

    std::cout << "Generazione " << gen + 1 << "/" << num_generations
              << " completata. Fronti di Pareto: " << combined_fronts.size()
              << std::endl;
  }
}

int main(int argc, char *argv[]) {

  if (argc != 2) {
    std::cout << "Inserisci il numero degli workers" << std::endl;
    exit(1);
  }

  int nw = std::atoi(argv[1]);

  TIMERSTART(data_loading);
  FMODataManager manager = getDataFromPathAndAngles(PATH, GANTRIES);
  // manager.printSummary();
  TIMERSTOP(data_loading);

  // generazione della popolazione iniziale
  TIMERSTART(inizial_population);
  std::mt19937 rng(42); // Inizializza il generatore di numeri casuali con un
                        // seed fisso per la riproducibilità
  Population start =
      generateRandomPopulation(100, manager.getTotalBeamlets(), rng);
  TIMERSTOP(inizial_population);

  // Esecuzione dell'algoritmo NSGA-II per un numero prefissato di generazioni
  TIMERSTART(nsga2seq);
  Evaluator evaluator(manager.getData());
  nsga2Omp(start, 50, 100, 0.9, 0.1, 20.0, 20.0, evaluator, rng, nw);
  TIMERSTOP(nsga2seq);

  return 0;
}