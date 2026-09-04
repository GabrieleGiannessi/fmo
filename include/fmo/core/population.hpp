/**
 * @file population.hpp
 * @brief Definizione della classe Population per rappresentare la popolazione
 * genetica
 * @details La classe Population rappresenta la popolazione genetica, contenente
 * un vettore di individui.
 * Contiene metodi per l'aggiunta e l'accesso agli individui nella popolazione.
 */
#pragma once
#include "fmo/core/individual.hpp"
#include <stdexcept>
#include <vector>

class Population {
public:
  std::vector<Individual> individuals;

  Population() = default;

  // Riserva la memoria iniziale senza imporre un limite rigido
  explicit Population(size_t reserve_count) {
    individuals.reserve(reserve_count);
  }

  // Inserimento elementi
  inline void addIndividual(const Individual &individual) {
    individuals.push_back(individual);
  }

  inline void addIndividuals(const std::vector<Individual> &new_individuals) {
    individuals.insert(individuals.end(), new_individuals.begin(),
                       new_individuals.end());
  }

  // Accesso sicuro (per riferimento per evitare copie)
  inline Individual &getIndividual(size_t index) {
    if (index >= individuals.size()) {
      throw std::out_of_range("Index out of range");
    }
    return individuals[index];
  }

  inline const Individual &getIndividual(size_t index) const {
    if (index >= individuals.size()) {
      throw std::out_of_range("Index out of range");
    }
    return individuals[index];
  }

  inline void removeIndividual(size_t index) {
    if (index >= individuals.size()) {
      throw std::out_of_range("Index out of range");
    }
    individuals.erase(individuals.begin() + index);
  }

  inline size_t size() const { return individuals.size(); }
  inline bool empty() const { return individuals.empty(); }
  inline void clear() { individuals.clear(); }
  inline void reserve(size_t n) { individuals.reserve(n); }
};