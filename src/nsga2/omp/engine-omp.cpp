
#include "fmo/nsga2/nsga-utils.hpp"
#include "fmo/nsga2/nsga2-steps-omp.hpp"
#include "fmo/nsga2/offspring.hpp"
#include "fmo/utilities/hpc_helpers.hpp"

#ifdef _OPENMP
#include "omp.h"
#endif

/**
 * @brief valutazione delle fitness per un insieme di individui (popolazione) in
 * parallelo (embarassingly parallel). 
 */
inline void evaluatePopulationOmp(Population &population, Evaluator &evaluator,
                                  int nw) {
  TIMERSTART(population_evaluation);
#pragma omp parallel for num_threads(nw)
  for (int i = 0; i < static_cast<int>(population.size()); ++i) {
    auto &individual = population.getIndividual(i);
    individual.setFitness(evaluator.evaluate(individual));
  }
  TIMERSTOP(population_evaluation);
}

/**
 * @brief generazione della nuova prole distribuendo il calcolo dei nuovi
 * individui su più cores (nw). In particolare, ogni thread riceve il seme base
 * e calcola il suo PRNG per effettuare poi le chiamate ai metodi di selezione a
 * torneo, crossover e mutazione polinomiale.
 */
inline Population generatePopulationOffspringOmp(const Population &population,
                                                 double eta_c, double eta_m,
                                                 uint64_t base_seed, int nw) {
  TIMERSTART(offspring_generation);
  const size_t N = population.size();
  Population offspring(N);

  // Pre-allochiamo lo spazio per evitare reallocazioni concorrenti
  for (size_t i = 0; i < N; ++i) {
    offspring.addIndividual(population.getIndividual(i));
  }

// Creazione della regione parallela
#pragma omp parallel num_threads(nw)
  {
    int tid = 0;
#ifdef _OPENMP
    int tid = omp_get_thread_num();
#endif
    // PRNG privato per thread con seed deterministico derivato
    std::mt19937 thread_rng(base_seed + tid * 10007);

// Partizionamento statico del loop
#pragma omp for schedule(static)
    for (int i = 0; i < static_cast<int>(N); i += 2) {

      const Individual &parent1 = GeneticOperator::tournament_selection(
          population.individuals, thread_rng);
      const Individual &parent2 = GeneticOperator::tournament_selection(
          population.individuals, thread_rng);

      auto [child1, child2] = GeneticOperator::crossover(
          parent1, parent2, CROSSOVER_PROB, eta_c, thread_rng);

      GeneticOperator::mutate(child1, 1.0 / static_cast<double>(child1.size()),
                              eta_m, thread_rng);
      GeneticOperator::mutate(child2, 1.0 / static_cast<double>(child2.size()),
                              eta_m, thread_rng);

      // Scrittura sicura su indici di memoria esclusivi per il task corrente
      offspring.getIndividual(i) = child1;
      if (i + 1 < static_cast<int>(N)) {
        offspring.getIndividual(i + 1) = child2;
      }
    }
  }
  TIMERSTOP(offspring_generation);
  return offspring;
}

std::vector<std::vector<int>> sortPopulation(Population &population) {
  TIMERSTART(population_sorting);
  auto fronts = fastNondominatedSort(population);
  TIMERSTOP(population_sorting);
  return fronts;
}

void assignPopulationCrowding(Population &population,
                              const std::vector<std::vector<int>> &fronts) {
  TIMERSTART(population_crowding);
  for (const auto &front : fronts) {
    assignCrowdingDistance(front, population);
  }
  TIMERSTOP(population_crowding);
}

Population mergePopulations(const Population &population,
                            const Population &offspring) {
  TIMERSTART(population_merge);
  Population combined_population(population.size() + offspring.size());
  combined_population.addIndividuals(population.individuals);
  combined_population.addIndividuals(offspring.individuals);
  TIMERSTOP(population_merge);
  return combined_population;
}

Population
truncatePopulationByFronts(Population &population,
                           const std::vector<std::vector<int>> &fronts,
                           int population_size) {
  TIMERSTART(population_truncation);
  Population truncated_population =
      truncatePopulation(population, fronts, population_size);
  TIMERSTOP(population_truncation);
  return truncated_population;
}