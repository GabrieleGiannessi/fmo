/**
 * @file nsga2-ff.cpp
 * @brief Implementazione dell'algoritmo NSGA-II (Deb et al., 2002) con
 * parallelizzazione tramite la libreria FastFlow per la generazione della prole
 * e le valutazioni delle fitness degli individui della popolazione.
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

/**
 * @brief esegue l'algoritmo NSGA-II (Deb et al., 2002) usando la libreria
 * FastFlow per parallelizzare delle parti di calcolo come la generazione delle
 * popolazioni successive e la valutazione delle fitness degli individui della
 * popolazione.
 *
 * @param pop: popolazione iniziale da cui parte l'algoritmo
 * @param num_generations numero di generazioni su cui iterare
 * @param population_size dimensione della popolazione iniziare e di quelle da
 * generare
 * @param eta_c Indice di distribuzione del crossover, controlla quanto i
 * figli si "allontanano" dai genitori. Valori più alti producono figli più
 * vicini ai genitori
 * @param eta_m Indice di distribuzione della mutazione
 * @param evaluator oggetto Evaluator, usato per calcolare le Fitness degli
 * individui
 * @param rng pseudo-random number generator
 * @param nw numero di workers
 *
 */
void nsga2ff(Population &pop, int num_generations, int population_size,
             double eta_c, double eta_m, Evaluator &evaluator,
             std::mt19937 &rng, int nw);

int main(int argc, char *argv[]) {

  if (argc != 2) {
    std::cout << "Inserisci il numero degli workers" << std::endl;
    exit(1);
  }

  int nw = std::atoi(argv[1]);
  if (nw < 1) {
    std::cout << "Inserisci un numero valido di workers" << std::endl;
    exit(1);
  }

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
}