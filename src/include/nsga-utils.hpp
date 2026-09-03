/**
 * @file nsga-utils.hpp
 * @brief Funzioni di utilità per l'algoritmo NSGA-II
 * @details Contiene strutture e predicati di dominanza per il fast
 * non-dominated sorting.
 */
#pragma once

#include "fitness.hpp"
#include "population.hpp"

#include <vector>

struct SortNode {
  int n_p =
      0; // Numero di individui della popolazione che dominano questa soluzione
  std::vector<int> s_p; // Indici degli individui della popolazione dominati da
                        // questa soluzione
};

/**
 * @brief Verifica se l'individuo 'a' domina 'b' in senso di Pareto.
 * @details Semantica:
 *  - PTV: Massimizzazione
 *  - Retto (OAR): Minimizzazione
 *  - Vescica (OAR): Minimizzazione
 * @return true se 'a' domina 'b', false altrimenti.
 */
inline bool dominates(const Fitness &a, const Fitness &b) {
  // 1. Condizione di "non peggiore": se 'a' è peggiore di 'b' in anche solo
  // uno, NON può dominare b
  if (a.value_target_ptv < b.value_target_ptv)
    return false; // PTV peggiore (minore)
  if (a.value_oar_rectal > b.value_oar_rectal)
    return false; // Retto peggiore (maggiore)
  if (a.value_oar_bladder > b.value_oar_bladder)
    return false; // Vescica peggiore (maggiore)

  // 2. Condizione di "strettamente migliore": 'a' deve superare 'b' in almeno
  // un obiettivo
  bool at_least_one_better = (a.value_target_ptv > b.value_target_ptv) ||
                             (a.value_oar_rectal < b.value_oar_rectal) ||
                             (a.value_oar_bladder < b.value_oar_bladder);

  return at_least_one_better;
}

/**
 * @brief Esegue il Fast Non-dominated Sorting
 * @details Questa funzione implementa l'algoritmo di ordinamento non dominato
 * veloce per classificare gli individui della popolazione in fronti di Pareto.
 * Ogni fronte contiene individui che non dominano tra loro, e il primo fronte
 * contiene gli individui non dominati da nessun altro.
 * @param population Riferimento alla popolazione da ordinare
 * @return Un vettore di fronti di Pareto, dove ogni fronte è un vettore di
 * indici degli individui nella popolazione.
 */
inline std::vector<std::vector<int>>
fastNondominatedSort(Population &population) {

  size_t N = population.size();
  std::vector<SortNode> nodes(N);

  std::vector<std::vector<int>> pareto_fronts;
  pareto_fronts.emplace_back();

  // Inizializzo i nodi di ordinamento non dominato
  for (size_t p = 0; p < N; ++p) {
    nodes[p].n_p = 0;
    nodes[p].s_p.clear();
  }

  // Riempio il vettore di SortNodes per verificare le dominanze tra gli
  // individui della popolazione O(MN^2)
  for (size_t p = 0; p < N; ++p) {
    for (size_t q = 0; q < N; ++q) {
      if (p == q)
        continue;

      if (dominates(population.getIndividual(p).getFitness(),
                    population.getIndividual(q).getFitness())) {
        // p domina q: aggiungo q alla lista delle soluzioni dominate da p
        nodes[p].s_p.push_back(q);
      } else if (dominates(population.getIndividual(q).getFitness(),
                           population.getIndividual(p).getFitness())) {
        // q domina p: incremento il contatore di quante soluzioni dominano p
        nodes[p].n_p++;
      }
    }
  }

  // Costruisco il fronte di Pareto F_1
  for (size_t p = 0; p < N; ++p) {
    if (nodes[p].n_p == 0) {
      // Se nessuno domina p, allora p appartiene al primo fronte
      pareto_fronts[0].push_back(p);
      population.getIndividual(p).setRank(1);
    }
  }

  int current_front_idx = 0;

  // Costruisco i fronti successivi F_2, F_3, ... partendo da F_1
  while (current_front_idx < static_cast<int>(pareto_fronts.size()) &&
         !pareto_fronts[current_front_idx].empty()) {

    std::vector<int> next_front; // prossimo fronte di Pareto da costruire

    for (int p : pareto_fronts[current_front_idx]) {
      for (int q : nodes[p].s_p) {
        nodes[q].n_p--;
        if (nodes[q].n_p == 0) {
          // Se nessuno domina più q, allora q appartiene al prossimo fronte
          next_front.push_back(q);
          population.getIndividual(q).setRank(current_front_idx + 2);
        }
      }
    }

    if (!next_front.empty()) {
      pareto_fronts.push_back(next_front);
    }
    current_front_idx++; // passo al prossimo fronte (quello appena generato)
  }

  return pareto_fronts;
}

/**
 * @brief Calcola e assegna la crowding distance agli individui di un fronte di
 * Pareto.
 *
 * @details Implementa la procedura di Deb et al. (2002) per stimare la densità
 * delle soluzioni:
 * 1. Inizializza le distanze a 0.0 (o a infinito se il fronte ha dimensione <=
 * 2).
 * 2. Per ciascuno dei 3 obiettivi clinici (PTV, Retto, Vescica):
 *    - Ordina gli individui in ordine crescente di obiettivo.
 *    - Assegna distanza infinita agli estremi per preservare i confini della
 * frontiera.
 *    - Per i punti interni, accumula la differenza normalizzata rispetto al
 * range (f_max - f_min).
 *
 * I valori calcolati vengono salvati in-place negli individui e usati da
 * crowded_compare.
 *
 * @param[in] front Vettore di indici delle soluzioni appartenenti al
 * fronte.
 * @param[in,out] population Popolazione in cui aggiornare il membro
 * crowding_distance.
 */
inline void assignCrowdingDistance(const std::vector<int> &front,
                                   Population &population) {
  size_t l = front.size();
  if (l == 0)
    return;

  if (l <= 2) {
    for (int idx : front) {
      population.getIndividual(idx).setCrowdingDistance(
          std::numeric_limits<double>::infinity());
    }
    return;
  }

  for (int idx : front) {
    population.getIndividual(idx).setCrowdingDistance(0.0);
  }

  std::vector<int> sorted_indices = front;

  for (int m = 0; m < 3; ++m) {
    // m = 0 (PTV), m = 1 (Rectum), m=2 (Bladder)
    auto get_obj = [m](const Fitness &f) -> double {
      if (m == 0)
        return f.value_target_ptv;
      if (m == 1)
        return f.value_oar_rectal;
      return f.value_oar_bladder;
    };

    std::sort(sorted_indices.begin(), sorted_indices.end(), [&](int a, int b) {
      return get_obj(population.getIndividual(a).getFitness()) <
             get_obj(population.getIndividual(b).getFitness());
    });

    // Assegnazione di infinito ai limiti della frontiera
    population.getIndividual(sorted_indices.front()).crowding_distance =
        std::numeric_limits<double>::infinity();
    population.getIndividual(sorted_indices.back()).crowding_distance =
        std::numeric_limits<double>::infinity();

    double f_min =
        get_obj(population.getIndividual(sorted_indices.front()).getFitness());
    double f_max =
        get_obj(population.getIndividual(sorted_indices.back()).getFitness());
    double denom = f_max - f_min;

    if (denom <= 1e-9)
      continue;

    for (size_t i = 1; i < l - 1; ++i) {
      int prev_idx = sorted_indices[i - 1];
      int next_idx = sorted_indices[i + 1];
      int curr_idx = sorted_indices[i];

      double dist = (get_obj(population.getIndividual(next_idx).getFitness()) -
                     get_obj(population.getIndividual(prev_idx).getFitness())) /
                    denom;

      if (population.getIndividual(curr_idx).getCrowdingDistance() !=
          std::numeric_limits<double>::infinity()) {
        population.getIndividual(curr_idx).setCrowdingDistance(
            population.getIndividual(curr_idx).getCrowdingDistance() + dist);
      }
    }
  }
}