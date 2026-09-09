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
#include "fmo/metrics/metrics.hpp"
#include "fmo/nsga2/nsga-utils.hpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>
#include <memory>

#include <ff/parallel_for.hpp>

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

  std::unique_ptr<ff::ParallelFor> pf;
  if (mod == 0) {
    pf = std::make_unique<ff::ParallelFor>(nw, true);
  }

  if (mod == 0) {
    measure("initial_evaluation",
            [&] { evaluatePopulationFFParFor(pop, evaluator, *pf); });
  } else {
    measure("initial_evaluation",
            [&] { evaluatePopulationFFFarm(pop, evaluator, nw); });
  }

  // Classificazione iniziale di P_0
  std::vector<std::vector<int>> fronts;
  measure("initial_sorting", [&] { fronts = sortPopulation(pop); });

  // assegnazione distance crowding iniziale alla P_0
  measure("initial_crowding", [&] { assignPopulationCrowding(pop, fronts); });

  // 2. Loop Generazionale
  for (int gen = 0; gen < num_generations; ++gen) {
    long generation_ms = 0;
    {
      utimer generation_timer("generation_total", &generation_ms, true);

      // A. Generazione discendenza Q_t (taglia N) tramite Torneo, SBX e
      // Mutazione
      Population offspring;
      measure("offspring_generation", [&] {
        offspring = generatePopulationOffspring(pop, eta_c, eta_m, rng);
      });

      // B. Valutazione della discendenza Q_t (calcolo delle fitness)
      if (mod == 0) {
        measure("population_evaluation",
                [&] { evaluatePopulationFFParFor(offspring, evaluator, *pf); });
      } else {
        measure("population_evaluation",
                [&] { evaluatePopulationFFFarm(offspring, evaluator, nw); });
      }

      // C. Fusione R_t = P_t U Q_t (taglia 2N)
      Population combined_pop;
      measure("population_merge",
              [&] { combined_pop = mergePopulations(pop, offspring); });

      // D. Non-dominated sorting ed estrazione dei fronti su R_t
      std::vector<std::vector<int>> combined_fronts;
      measure("population_sorting",
              [&] { combined_fronts = sortPopulation(combined_pop); });

      measure("population_crowding",
              [&] { assignPopulationCrowding(combined_pop, combined_fronts); });

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

void printUsage(const char *prog_name) {
  std::cout << "Uso: " << prog_name << " <workers> <mode> [pop_size] [generations] [pop_out_path]\n"
            << "   o: " << prog_name << " [OPZIONI]\n"
            << "Esegue l'algoritmo NSGA-II con parallelizzazione FastFlow (parfor o farm) per il problema FMO.\n\n"
            << "Opzioni:\n"
            << "  -w, --workers <NUM>      Numero di workers FastFlow\n"
            << "  -m, --mode <0|1>         Modalità di esecuzione (0: parfor, 1: farm)\n"
            << "  -p, --pop-size <N>       Dimensione della popolazione (default: 256)\n"
            << "  -g, --generations <G>    Numero di generazioni (default: 50)\n"
            << "  -o, --output <FILE>      Percorso per salvare la popolazione finale in CSV\n"
            << "  -h, --help               Mostra questo messaggio di aiuto ed esce\n\n"
            << "Variabili d'ambiente (usate se non specificate da CLI):\n"
            << "  FMO_POPULATION_SIZE, FMO_GENERATIONS, FMO_POPULATION_OUT\n";
}

int main(int argc, char *argv[]) {
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      printUsage(argv[0]);
      return 0;
    }
  }

  int nw = -1;
  int mod = -1;

  int population_size = 256;
  if (const char *env_p = std::getenv("FMO_POPULATION_SIZE")) {
    int p = std::atoi(env_p);
    if (p > 0) population_size = p;
  }
  int num_generations = 50;
  if (const char *env_g = std::getenv("FMO_GENERATIONS")) {
    int g = std::atoi(env_g);
    if (g > 0) num_generations = g;
  }

  std::string pop_out_path;
  if (const char *env_o = std::getenv("FMO_POPULATION_OUT")) {
    pop_out_path = env_o;
  }

  std::vector<std::string> pos_args;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if ((arg == "-w" || arg == "--workers") && i + 1 < argc) {
      nw = std::atoi(argv[++i]);
    } else if ((arg == "-m" || arg == "--mode") && i + 1 < argc) {
      std::string m_str = argv[++i];
      if (m_str == "parfor" || m_str == "0") mod = 0;
      else if (m_str == "farm" || m_str == "1") mod = 1;
      else mod = std::atoi(m_str.c_str());
    } else if ((arg == "-p" || arg == "--pop-size" || arg == "--population") && i + 1 < argc) {
      population_size = std::atoi(argv[++i]);
    } else if ((arg == "-g" || arg == "--generations") && i + 1 < argc) {
      num_generations = std::atoi(argv[++i]);
    } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
      pop_out_path = argv[++i];
    } else if (!arg.empty() && arg[0] == '-') {
      std::cerr << "Opzione sconosciuta: " << arg << "\n";
      printUsage(argv[0]);
      return 1;
    } else {
      pos_args.push_back(arg);
    }
  }

  auto is_number = [](const std::string &s) {
    return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
  };

  int pos_idx = 0;
  if (nw < 1 && !pos_args.empty() && is_number(pos_args[pos_idx])) {
    nw = std::atoi(pos_args[pos_idx++].c_str());
  }
  if (mod < 0 && pos_idx < static_cast<int>(pos_args.size())) {
    std::string m_str = pos_args[pos_idx++];
    if (m_str == "parfor" || m_str == "0") mod = 0;
    else if (m_str == "farm" || m_str == "1") mod = 1;
    else mod = std::atoi(m_str.c_str());
  }

  if (nw < 1) {
    std::cerr << "Errore: inserisci un numero valido di workers.\n";
    printUsage(argv[0]);
    return 1;
  }
  if (mod != 0 && mod != 1) {
    std::cerr << "Errore: solo due modalità di esecuzione disponibili (0: parfor, 1: farm).\n";
    printUsage(argv[0]);
    return 1;
  }

  if (pos_idx < static_cast<int>(pos_args.size())) {
    if (is_number(pos_args[pos_idx])) {
      population_size = std::atoi(pos_args[pos_idx++].c_str());
      if (pos_idx < static_cast<int>(pos_args.size()) && is_number(pos_args[pos_idx])) {
        num_generations = std::atoi(pos_args[pos_idx++].c_str());
      }
      if (pos_idx < static_cast<int>(pos_args.size())) {
        pop_out_path = pos_args[pos_idx++];
      }
    } else {
      pop_out_path = pos_args[pos_idx++];
    }
  }

  if (population_size <= 0) {
    std::cerr << "Errore: la dimensione della popolazione deve essere > 0 (fornita: " << population_size << ")\n";
    return 1;
  }
  if (num_generations <= 0) {
    std::cerr << "Errore: il numero di generazioni deve essere > 0 (fornito: " << num_generations << ")\n";
    return 1;
  }

  const char *mode_name = mod == 0 ? "parfor" : "farm";
  std::map<std::string, std::vector<long>> samples;

  std::cout << "RUN variant=ff mode=" << mode_name << " workers=" << nw
            << " generations=" << num_generations
            << " population=" << population_size << std::endl;

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
    nsga2ff(start, num_generations, population_size, 20.0, 20.0, evaluator, rng,
            nw, mod, samples);
  }
  samples["nsga2_total"].push_back(nsga2_total_ms);

  double spread = computeSpreadMetric(start);
  std::cout << "QUALITY variant=ff mode=" << mode_name
            << " spread_metric=" << spread << std::endl;

  if (!pop_out_path.empty()) {
    savePopulation(start, pop_out_path);
    std::cout << "SAVED_POPULATION variant=ff mode=" << mode_name << " path="
              << pop_out_path << std::endl;
  }

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
    const double mean = static_cast<double>(std::accumulate(
                            filtered.begin(), filtered.end(), 0L)) /
                        filtered.size();
    const double median = filtered.size() % 2 == 0
                              ? (filtered[filtered.size() / 2 - 1] +
                                 filtered[filtered.size() / 2]) /
                                    2.0
                              : filtered[filtered.size() / 2];
    std::cout << "TIMING variant=ff mode=" << mode_name << " phase=" << phase
              << " samples=" << filtered.size() << " mean_ms=" << mean
              << " median_ms=" << median << std::endl;
  }

  return 0;
}