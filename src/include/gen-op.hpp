/**
 * @file gen-op.hpp
 * @brief Definizione della classe GeneticOperator per le operazioni genetiche
 * utilizzate nell'algoritmo NSGA-II (Deb et al., 2002) applicato a FMO.
 * @details La classe espone metodi statici e stateless per:
 * - Crowded-Comparison Operator (<_n) per il confronto secondo Pareto e
 * densità.
 * - Selezione a torneo binario affollato (Binary Tournament Selection).
 * - Crossover continuo SBX (Simulated Binary Crossover).
 * - Mutazione Polinomiale a valori reali.
 * - Correzione fisica dei vincoli di non-negatività (clipping).
 */

#pragma once

#include "individual.hpp"
#include <algorithm>
#include <random>
#include <utility>
#include <vector>

class GeneticOperator {
public:
  /**
   * @brief Crowded-Comparison Operator (<_n) definito da Deb et al. (2002).
   * @details Restituisce true se l'individuo 'a' domina o è preferibile a 'b':
   * 1. Priorità al rango di dominanza inferiore (rank minore batte rank
   * maggiore).
   * 2. A parità di rango, vince chi ha crowding distance maggiore (regione meno
   * densa).
   * @param a Primo individuo
   * @param b Secondo individuo
   * @return true se a è migliore di b, false altrimenti
   */
  static inline bool crowded_compare(const Individual &a, const Individual &b) {
    if (a.rank != b.rank) {
      return a.rank < b.rank;
    }
    return a.crowding_distance > b.crowding_distance;
  }

  /**
   * @brief Selezione a torneo binario per la scelta di un genitore.
   * @details Estrae uniformemente a caso due candidati dalla popolazione e
   * seleziona il vincitore tramite il Crowded-Comparison Operator
   * (crowded_compare).
   * @param population Popolazione da cui estrarre i candidati
   * @param rng Generatore di numeri pseudo-casuali (thread-safe)
   * @return Riferimento costante all'individuo genitore selezionato
   */
  static const Individual &
  tournament_selection(const std::vector<Individual> &population,
                       std::mt19937 &rng) {
    std::uniform_int_distribution<size_t> dist(0, population.size() - 1);

    size_t idx1 = dist(rng);
    size_t idx2 = dist(rng);

    const Individual &cand1 = population[idx1];
    const Individual &cand2 = population[idx2];

    return crowded_compare(cand1, cand2) ? cand1 : cand2;
  }

  /**
   * @brief Crossover SBX (Simulated Binary Crossover) per variabili continue.
   * @param parent1 Primo genitore
   * @param parent2 Secondo genitore
   * @param crossover_probability Probabilità di applicare il crossover (es.
   * 0.9)
   * @param eta_c Indice di distribuzione del crossover
   * @param rng Generatore pseudo-casuale
   * @return Coppia di nuovi individui figli
   */
  static std::pair<Individual, Individual>
  crossover(const Individual &parent1, const Individual &parent2,
            double crossover_probability, double eta_c, std::mt19937 &rng);

  /**
   * @brief Mutazione Polinomiale dei geni dell'individuo.
   * @param individual Riferimento all'individuo da mutare in-place
   * @param mutation_probability Probabilità di mutazione per ciascun gene
   * (es. 1.0 / num_bixels)
   * @param eta_m Indice di distribuzione della mutazione (default: 20.0)
   * @param rng Generatore pseudo-casuale
   * @param min_val Limite inferiore fisico per gene (default: 0.0)
   * @param max_val Limite superiore fisico per gene (default: 100.0)
   */
  static void mutate(Individual &individual, double mutation_probability,
                     double eta_m, std::mt19937 &rng, double min_val = 0.0,
                     double max_val = 100.0);

  /**
   * @brief Correzione fisica dei vincoli di non-negatività (clipping).
   * @details Garantisce che nessun bixel assuma intensità negative dopo
   * mutazione o crossover.
   * @param ind Individuo da correggere
   * @param min_val Valore minimo consentito (default: 0.0)
   * @param max_val Valore massimo consentito (default: 100.0)
   */
  static void clip(Individual &ind, double min_val = 0.0,
                   double max_val = 100.0) {
    for (double &gene : ind.genes) {
      gene = std::max(min_val, std::min(gene, max_val));
    }
  }
};