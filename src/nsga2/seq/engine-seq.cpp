#include "fmo/nsga2/nsga2-steps.hpp"

#include "fmo/nsga2/nsga-utils.hpp"
#include "fmo/nsga2/offspring.hpp"
#include "fmo/utilities/hpc_helpers.hpp"

void evaluatePopulation(Population &population, Evaluator &evaluator) {
  for (size_t i = 0; i < population.size(); ++i) {
    Individual &individual = population.getIndividual(i);
    individual.setFitness(evaluator.evaluate(individual));
  }
}

std::vector<std::vector<int>> sortPopulation(Population &population) {
  auto fronts = fastNondominatedSort(population);
  return fronts;
}

void assignPopulationCrowding(Population &population,
                              const std::vector<std::vector<int>> &fronts) {

  for (const auto &front : fronts) {
    assignCrowdingDistance(front, population);
  }
}

Population generatePopulationOffspring(Population &population, double eta_c,
                                       double eta_m, std::mt19937 &rng) {

  Population offspring = generateOffSpring(population, eta_c, eta_m, rng);
  return offspring;
}

Population mergePopulations(const Population &population,
                            const Population &offspring) {

  Population combined_population(population.size() + offspring.size());
  combined_population.addIndividuals(population.individuals);
  combined_population.addIndividuals(offspring.individuals);
  return combined_population;
}

Population
truncatePopulationByFronts(Population &population,
                           const std::vector<std::vector<int>> &fronts,
                           int population_size) {

  Population truncated_population =
      truncatePopulation(population, fronts, population_size);
  return truncated_population;
}
