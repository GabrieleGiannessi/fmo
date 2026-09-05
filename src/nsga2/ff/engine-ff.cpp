#include <ff/ff.hpp>
#include <ff/node.hpp>
#include <ff/parallel_for.hpp>

#include "fmo/evaluation/evaluator.hpp"
#include "fmo/nsga2/nsga-utils.hpp"
#include "fmo/nsga2/offspring.hpp"
#include "fmo/utilities/hpc_helpers.hpp"

using namespace ff;

// devo creare i nodi per l'elaborazione della fitness dato un dato in input
// (Individuo)
class FitnessWorker : public ff_node_t<Individual, Individual> {
private:
  Evaluator &evaluator;

public:
  FitnessWorker(Evaluator &e) : evaluator(e) {}

  Individual *svc(Individual *ind) {
    if (ind == nullptr) {
      return EOS;
    }
    ind->setFitness(evaluator.evaluate(*ind));
    return ind;
  }
};

ParallelFor map_fitness(
    8); // map pattern per il calcolo delle fitness in parallelo sugli individui

void evaluatePopulationFF(Population &population, Evaluator &evaluator, int nw);

/**
 * @brief generazione della nuova prole distribuendo il calcolo dei nuovi
 * individui su più cores (nw). In particolare, ogni thread riceve il seme base
 * e calcola il suo PRNG per effettuare poi le chiamate ai metodi di selezione a
 * torneo, crossover e mutazione polinomiale.
 */
Population generatePopulationOffspringFF(const Population &population,
                                         double eta_c, double eta_m,
                                         uint64_t base_seed, int nw);

// void assignPopulationCrowding(Population &population,
//                               const std::vector<std::vector<int>> &fronts) {
//   TIMERSTART(population_crowding);
//   for (const auto &front : fronts) {
//     assignCrowdingDistance(front, population);
//   }
//   TIMERSTOP(population_crowding);
// }
//
// Population mergePopulations(const Population &population,
//                             const Population &offspring) {
//   TIMERSTART(population_merge);
//   Population combined_population(population.size() + offspring.size());
//   combined_population.addIndividuals(population.individuals);
//   combined_population.addIndividuals(offspring.individuals);
//   TIMERSTOP(population_merge);
//   return combined_population;
// }
//
// Population
// truncatePopulationByFronts(Population &population,
//                            const std::vector<std::vector<int>> &fronts,
//                            int population_size) {
//   TIMERSTART(population_truncation);
//   Population truncated_population =
//       truncatePopulation(population, fronts, population_size);
//   TIMERSTOP(population_truncation);
//   return truncated_population;
// }
