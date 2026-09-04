#include "fmo/nsga2/nsga2-steps.hpp"

#include "fmo/nsga2/nsga-utils.hpp"
#include "fmo/nsga2/offspring.hpp"
#include "fmo/utilities/hpc_helpers.hpp"

void evaluatePopulation(Population &population, Evaluator &evaluator) {
  TIMERSTART(population_evaluation);
  for (size_t i = 0; i < population.size(); ++i) {
    Individual &individual = population.getIndividual(i);
    individual.setFitness(evaluator.evaluate(individual));
  }
  TIMERSTOP(population_evaluation);
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
