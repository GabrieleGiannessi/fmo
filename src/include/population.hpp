/**
 * @file population.hpp
 * @brief Definizione della classe Population per rappresentare la popolazione
 * genetica
 * @details La classe Population rappresenta la popolazione genetica, contenente
 * un vettore di individui.
 * Contiene metodi per l'aggiunta e l'accesso agli individui nella popolazione.
 */

#pragma once
#include "individual.hpp"
#include <algorithm>
#include <stdexcept>
#include <vector>

class Population {
public:
  std::vector<Individual> individuals;
  size_t capacity;

  explicit Population(size_t cap) : capacity(cap) { individuals.reserve(cap); }

  // Aggiunta di un individuo
  inline void addIndividual(const Individual &individual) {
    if (individuals.size() >= capacity) {
      throw std::runtime_error("Population is full");
    }
    individuals.push_back(individual);
  }

  // Aggiunta di n individui
  inline void addIndividuals(const std::vector<Individual> &new_individuals) {
    if (individuals.size() + new_individuals.size() > capacity) {
      throw std::runtime_error(
          "Adding these individuals would exceed population capacity");
    }
    individuals.insert(individuals.end(), new_individuals.begin(),
                       new_individuals.end());
  }

  // Rimozione di un individuo
  inline void removeIndividual(size_t index) {
    if (index >= individuals.size()) {
      throw std::out_of_range("Index out of range");
    }
    individuals.erase(individuals.begin() + index);
  }

  // Accesso sicuro a un individuo
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

  inline size_t size() const { return individuals.size(); }
  inline bool empty() const { return individuals.empty(); }
  inline void clear() { individuals.clear(); }
};