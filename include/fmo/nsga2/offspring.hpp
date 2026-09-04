/**
 * @file offspring.hpp
 * @brief definizione dell'algoritmo funzionale per la generazione della prole a
 * partire da una popolazione di individui.
 * @details La funzione generateOffspring prende in input una popolazione di
 * individui e genera una nuova popolazione di figli utilizzando le operazioni
 * genetiche di crossover e mutazione presenti in GeneticOperator (gen-op.hpp).
 * La funzione seleziona casualmente due genitori dalla popolazione, applica
 * l'operatore di crossover per generare due figli e successivamente applica
 * l'operatore di mutazione su ciascun figlio. La nuova popolazione di figli
 * viene restituita come output della funzione.
 */
#pragma once
#include "fmo/core/individual.hpp"
#include "fmo/core/population.hpp"
#include "fmo/nsga2/gen-op.hpp"
#include <vector>

#define CROSSOVER_PROB 0.9 // probabilità definita da paper (Deb et al., 2002)

inline Population generateRandomPopulation(int population_size, int num_bixels,
                                           std::mt19937 &rng) {
  Population population(population_size);
  std::uniform_real_distribution<double> dist(
      0.0, 10.0); // la popolazione iniziale contiene individui con geni i cui
                  // valori reali appartengono a [0, 10]

  for (int i = 0; i < population_size; ++i) {
    Individual individual(num_bixels);
    for (int j = 0; j < num_bixels; ++j) {
      individual.genes[j] = dist(rng);
    }
    population.addIndividual(individual);
  }

  return population;
}

inline Population generateOffSpring(Population &population, double eta_c,
                                    double eta_m, std::mt19937 &rng) {
  Population offspring(population.size());
  for (int i = 0; i < population.size() / 2; ++i) {
    // Selezione dei genitori tramite torneo binario
    const Individual &parent1 =
        GeneticOperator::tournament_selection(population.individuals, rng);
    const Individual &parent2 =
        GeneticOperator::tournament_selection(population.individuals, rng);

    // Crossover per generare due figli
    auto [child1, child2] = GeneticOperator::crossover(
        parent1, parent2, CROSSOVER_PROB, eta_c, rng);

    // Mutazione dei figli
    GeneticOperator::mutate(child1, 1.0 / static_cast<double>(child1.size()), eta_m, rng);
    GeneticOperator::mutate(child2, 1.0 / static_cast<double>(child2.size()), eta_m, rng);

    // Aggiunta dei figli alla nuova popolazione
    offspring.addIndividual(child1);
    offspring.addIndividual(child2);
  }
  return offspring;
}