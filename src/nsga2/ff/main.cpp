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
#include "fmo/nsga2/nsga2-steps-ff.hpp"
#include "fmo/nsga2/nsga2-steps.hpp"
#include "fmo/nsga2/offspring.hpp"
#include "fmo/preprocessing/manager.hpp"
#include "fmo/utilities/hpc_helpers.hpp"
#include "fmo/utilities/utimer.hpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
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
             std::mt19937 &rng, int nw, int mod,
             std::map<std::string, std::vector<long>> &samples) {
  auto measure = [&](const std::string &phase, auto operation) {
    long elapsed_ms = 0;
    {
      utimer timer(phase, &elapsed_ms, true);
      operation();
    }
    samples[phase].push_back(elapsed_ms);
  };

  if (mod == 0) {
    measure("initial_evaluation", [&] {
      evaluatePopulationFFParFor(pop, evaluator, nw);
    });
  } else {
    measure("initial_evaluation", [&] {
      evaluatePopulationFFFarm(pop, evaluator, nw);
    });
  }

  // Classificazione iniziale di P_0
  std::vector<std::vector<int>> fronts;
  measure("initial_sorting", [&] { fronts = sortPopulation(pop); });

  // assegnazione distance crowding iniziale alla P_0
  measure("initial_crowding", [&] {
    assignPopulationCrowding(pop, fronts);
  });

  // 2. Loop Generazionale
  for (int gen = 0; gen < num_generations; ++gen) {
    long generation_ms = 0;
    {
      utimer generation_timer("generation_total", &generation_ms, true);

      // A. Generazione discendenza Q_t (taglia N) tramite Torneo, SBX e Mutazione
      Population offspring;
      measure("offspring_generation", [&] {
        offspring = generatePopulationOffspring(pop, eta_c, eta_m, rng);
      });

      // B. Valutazione della discendenza Q_t (calcolo delle fitness)
      if (mod == 0) {
        measure("population_evaluation", [&] {
          evaluatePopulationFFParFor(offspring, evaluator, nw);
        });
      } else {
        measure("population_evaluation", [&] {
          evaluatePopulationFFFarm(offspring, evaluator, nw);
        });
      }

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

  if (argc != 3) {
    std::cout << "Inserisci il numero degli workers e la mod. di esecuzione"
              << std::endl;
    exit(1);
  }

  int nw = std::atoi(argv[1]);
  if (nw < 1) {
    std::cout << "Inserisci un numero valido di workers" << std::endl;
    exit(1);
  }

  int mod = std::atoi(argv[2]);
  if (mod != 0 && mod != 1) {
    std::cout
        << "Solo due modalità di esecuzione disponibili (0: parfor; 1: farm)"
        << std::endl;
    exit(1);
  }

  constexpr int population_size = 100;
  constexpr int num_generations = 50;
  const char *mode_name = mod == 0 ? "parfor" : "farm";
  std::map<std::string, std::vector<long>> samples;

  std::cout << "RUN variant=ff mode=" << mode_name << " workers=" << nw
            << " generations=" << num_generations << " population="
            << population_size << std::endl;

  long data_loading_ms = 0;
  FMODataManager manager;
  {
    utimer timer("data_loading", &data_loading_ms, true);
    manager = getDataFromPathAndAngles(PATH, GANTRIES);
  }
  samples["data_loading"].push_back(data_loading_ms);
  // manager.printSummary();

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

  long nsga2_total_ms = 0;
  Evaluator evaluator(manager.getData());
  {
    utimer timer("nsga2_total", &nsga2_total_ms, true);
    nsga2ff(start, num_generations, population_size, 20.0, 20.0, evaluator,
            rng, nw, mod, samples);
  }
  samples["nsga2_total"].push_back(nsga2_total_ms);

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
    std::cout << "TIMING variant=ff mode=" << mode_name << " phase="
              << phase << " samples=" << filtered.size()
              << " mean_ms=" << mean << " median_ms=" << median
              << std::endl;
  }

  return 0;

  return 0;
}