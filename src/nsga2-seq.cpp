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

#include "include/gen-op.hpp"
#include "include/hpc_helpers.hpp"
#include "include/individual.hpp"
#include "include/nsga-utils.hpp"
#include "include/offspring.hpp"
#include "include/population.hpp"

#include "preprocessing/manager.hpp"

#include <iostream>
#include <vector>
#define GANTRIES                                                               \
  {0, 72, 144, 216, 288}                                                       \
  // Angoli di gantry da caricare per la matrice globale D
#define PATH "../data/phantom" // Percorso base per i file di input

FMODataManager getDataFromPathAndAngles(const std::string &base_dir,
                                        const std::vector<int> &gantry_angles) {
  FMODataManager manager;
  manager.loadGlobalDMatrixFromAngles(base_dir, gantry_angles);
  return manager;
}

void nsga2seq(Population &pop, int num_generations, int population_size,
              double crossover_probability, double mutation_probability,
              double eta_c, double eta_m) {
  std::cout << "Esecuzione dell'algoritmo NSGA-II per " << num_generations
            << " generazioni..." << std::endl;
  // Implementazione dell'algoritmo NSGA-II sequenziale
  // Si divide l'intera popolazione in fronti di pareto che indicano il
  // livello di dominanza delle soluzioni.
}

int main(int argc, char *argv[]) {
  std::cout << "NSGA-II sequenziale" << std::endl;

  TIMERSTART(data_loading);
  FMODataManager manager = getDataFromPathAndAngles(PATH, GANTRIES);
  manager.printSummary();
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
  nsga2seq(start, 50, 100, 0.9, 0.1, 20.0, 20.0);
  TIMERSTOP(nsga2seq);

  return 0;
}
