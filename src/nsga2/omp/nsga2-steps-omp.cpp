
#include "fmo/nsga2/nsga-utils.hpp"
#include "fmo/nsga2/offspring.hpp"
#include "fmo/nsga2/nsga2-steps-omp.hpp"
#include "fmo/utilities/hpc_helpers.hpp"

// calcolo delle fitness in parallelo (embarassingly parallel)
void evaluatePopulationOmp(Population &population, Evaluator &evaluator) {
  TIMERSTART(population_evaluation);
#pragma omp parallel for
  for (int i = 0; i < static_cast<int>(population.size()); ++i) {
    auto &individual = population.getIndividual(i);
    individual.setFitness(evaluator.evaluate(individual));
  }
  TIMERSTOP(population_evaluation);
}

// generazione della prole parallela
Population generatePopulationOffspringOmp(Population &population,
                                          double crossover_probability,
                                          double mutation_probability,
                                          double eta_c, double eta_m,
                                          std::mt19937 &rng);

// assegnazione crowding distance parallelo
void assignPopulationCrowdingOmp(Population &population,
                                 const std::vector<std::vector<int>> &fronts) {
  TIMERSTART(crowding_assignment);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < static_cast<int>(fronts.size()); ++i) {
    assignCrowdingDistance(fronts[i], population);
  }
  TIMERSTOP(crowding_assignment);
}

// non-dominated sorting in parallelo
std::vector<std::vector<int>> sortPopulationOmp(Population &population);