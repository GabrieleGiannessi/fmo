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

#include "fmo/core/individual.hpp"
#include <algorithm>
#include <random>
#include <stdexcept>
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
   * @param eta_c Indice di distribuzione del crossover, controlla quanto i
   * figli si "allontanano" dai genitori. Valori più alti producono figli più
   * vicini ai genitori
   * @param rng Generatore pseudo-casuale
   * @return Coppia di nuovi individui figli
   */
  static std::pair<Individual, Individual>
  crossover(const Individual &parent1, const Individual &parent2,
            double crossover_probability, double eta_c, std::mt19937 &rng) {
    if (parent1.genes.size() != parent2.genes.size()) {
      throw std::invalid_argument(
          "I genitori devono avere lo stesso numero di geni");
    }

    Individual child1 = parent1;
    Individual child2 = parent2;

    std::uniform_real_distribution<double> dist(
        0.0, 1.0); // distribuzione uniforme per la probabilità di crossover e
                   // il calcolo di beta

    if (dist(rng) <= crossover_probability) {
      for (size_t i = 0; i < parent1.genes.size(); ++i) {
        // Swap probabilistico per gene (50%)
        if (dist(rng) <= 0.5) {
          double u = dist(rng);
          double beta; // fattore di diffusione per il crossover SBX
          if (u <= 0.5) {
            beta = std::pow(2.0 * u, 1.0 / (eta_c + 1.0));
          } else {
            beta = std::pow(1.0 / (2.0 * (1.0 - u)), 1.0 / (eta_c + 1.0));
          }

          double x1 = parent1.genes[i];
          double x2 = parent2.genes[i];

          child1.genes[i] = 0.5 * ((1.0 + beta) * x1 + (1.0 - beta) * x2);
          child2.genes[i] = 0.5 * ((1.0 - beta) * x1 + (1.0 + beta) * x2);
        }
      }
      // Applicazione del clipping sui figli generati
      clip(child1);
      clip(child2);
    }

    return std::make_pair(child1, child2);
  }

  /**
   * @brief Mutazione Polinomiale dei geni dell'individuo.
   * @details Per ogni gene, con probabilità mutation_probability, viene
   * applicata la mutazione polinomiale secondo l'indice di distribuzione eta_m.
   * La mutazione polinomiale genera un nuovo valore del gene vicino al valore
   * originale, con una distribuzione controllata da eta_m. Valori più alti
   * di eta_m producono mutazioni più vicine al gene originale, mentre valori
   * più bassi producono mutazioni più lontane.
   * La mutazione polinomiale è definita come segue:
   * - Genera un numero casuale u uniformemente distribuito tra 0 e 1.
   * - Calcola il fattore di mutazione delta in base a u e eta_m.
   * - Aggiorna il gene con il nuovo valore calcolato.
   * @note Dopo la mutazione, viene applicato il clipping per garantire che i
   * valori dei geni rimangano entro i limiti fisici specificati (min_val e
   * max_val).
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
                     double max_val = 100.0) {

    std::uniform_real_distribution<double> dist(0.0, 1.0);
    if (dist(rng) <= mutation_probability) {
      double u = dist(rng);
      double delta;
      if (u < 0.5) {
        delta = std::pow(2.0 * u, 1.0 / (eta_m + 1.0)) - 1.0;
      } else {
        delta = 1.0 - std::pow(2.0 * (1.0 - u), 1.0 / (eta_m + 1.0));
      }
      for (double &gene : individual.genes) {
        gene += delta * (max_val - min_val);
      }
    }
    // Applicazione del clipping dopo la mutazione
    clip(individual, min_val, max_val);
  }

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