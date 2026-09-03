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
#include "include/evaluator.hpp"

#include "preprocessing/manager.hpp"

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

void nsga2seq(Population &pop, int num_generations, int population_size,
              double crossover_probability, double mutation_probability,
              double eta_c, double eta_m, Evaluator &evaluator, std::mt19937 &rng) {
  
  std::cout << "Esecuzione dell'algoritmo NSGA-II per " << num_generations
            << " generazioni..." << std::endl;

  std::cout << "Valutazione iniziale degli individui" << std::endl; 

  // 1. Valutazione iniziale di P_0
  for (size_t i = 0; i < pop.size(); ++i) {
    pop.getIndividual(i).setFitness(evaluator.evaluate(pop.getIndividual(i)));
  }

  std::cout << "Classificazione degli individui" << std::endl; 
  // Classificazione iniziale di P_0
  auto fronts = fastNondominatedSort(pop);
  for (const auto &front : fronts) {
    assignCrowdingDistance(front, pop);
  }

  // 2. Loop Generazionale
  for (int gen = 0; gen < num_generations; ++gen) {

    std::cout << "Generazione " << gen << std::endl; 
    // A. Generazione discendenza Q_t (taglia N) tramite Torneo, SBX e Mutazione
    Population offspring = generateOffSpring(pop, crossover_probability,
                                             mutation_probability, eta_c, eta_m, rng);

    // B. Valutazione della discendenza Q_t
    for (size_t i = 0; i < offspring.size(); ++i) {
      offspring.getIndividual(i).setFitness(evaluator.evaluate(offspring.getIndividual(i)));
    }

    // C. Fusione R_t = P_t U Q_t (taglia 2N)
    Population combined_pop(population_size * 2);
    for (size_t i = 0; i < pop.size(); ++i) {
      combined_pop.addIndividual(pop.getIndividual(i));
    }
    for (size_t i = 0; i < offspring.size(); ++i) {
      combined_pop.addIndividual(offspring.getIndividual(i));
    }

    // D. Non-dominated sorting ed estrazione dei fronti su R_t
    auto combined_fronts = fastNondominatedSort(combined_pop);
    for (const auto &front : combined_fronts) {
      assignCrowdingDistance(front, combined_pop);
    }

    // E. Elitismo e troncamento: R_t -> P_{t+1} (taglia N)
    pop = truncatePopulation(combined_pop, combined_fronts, population_size);

    std::cout << "Generazione " << gen + 1 << "/" << num_generations 
              << " completata. Fronti di Pareto: " << combined_fronts.size() << std::endl;
  }
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
  Evaluator evaluator(manager.getData());
  nsga2seq(start, 50, 100, 0.9, 0.1, 20.0, 20.0, evaluator, rng);
  TIMERSTOP(nsga2seq);

  return 0;
}
