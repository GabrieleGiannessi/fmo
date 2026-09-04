
#include "fmo/nsga2/nsga-utils.hpp"
#include "fmo/nsga2/nsga2-steps-omp.hpp"
#include "fmo/nsga2/offspring.hpp"
#include "fmo/utilities/hpc_helpers.hpp"
#include "omp.h"

// calcolo delle fitness in parallelo (embarassingly parallel)
void evaluatePopulationOmp(Population &population, Evaluator &evaluator,
                           int nw) {
  TIMERSTART(population_evaluation);
#pragma omp parallel for num_threads(nw)
  for (int i = 0; i < static_cast<int>(population.size()); ++i) {
    auto &individual = population.getIndividual(i);
    individual.setFitness(evaluator.evaluate(individual));
  }
  TIMERSTOP(population_evaluation);
}

// generazione della prole parallela
// Population generatePopulationOffspringOmp(Population &population,
//                                           double crossover_probability,
//                                           double mutation_probability,
//                                           double eta_c, double eta_m,
//                                           std::mt19937 &rng) {
//   TIMERSTART(offspring_generation);
//   Population offspring =
//       generateOffSpring(population, crossover_probability,
//       mutation_probability,
//                         eta_c, eta_m, rng);
//   TIMERSTOP(offspring_generation);
//   return offspring;
// }

// assegnazione crowding distance parallelo
// void assignPopulationCrowdingOmp(Population &population,
//                                  const std::vector<std::vector<int>> &fronts)
//                                  {
//   TIMERSTART(crowding_assignment);
// #pragma omp parallel for schedule(dynamic)
//   for (int i = 0; i < static_cast<int>(fronts.size()); ++i) {
//     assignCrowdingDistance(fronts[i], population);
//   }
//   TIMERSTOP(crowding_assignment);
// }

// non-dominated sorting in parallelo
// std::vector<std::vector<int>> sortPopulationOmp(Population &population){

// }

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

Population generatePopulationOffspring(Population &population, double eta_c,
                                       double eta_m, std::mt19937 &rng) {
  TIMERSTART(offspring_generation);
  Population offspring = generateOffSpring(population, eta_c, eta_m, rng);
  TIMERSTOP(offspring_generation);
  return offspring;
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