/**
 * @file population.hpp
 * @brief Definizione della classe Population per rappresentare la popolazione
 * genetica
 * @details La classe Population rappresenta la popolazione genetica, contenente
 * un vettore di individui.
 * Contiene metodi per l'aggiunta e l'accesso agli individui nella popolazione.
 */

#include "individual.hpp"
#include <vector>
#include <algorithm>

class Population {
public:
  std::vector<Individual> individuals;
  int size;

  Population(int size) : size(size) { individuals.reserve(size); }

  // aggiunta di un individuo alla popolazione
  void addIndividual(const Individual &individual) {
    if (individuals.size() >= size) {
      throw std::runtime_error("Population is full");
    }
    individuals.push_back(individual);
  }

  // aggiunta di n individui alla popolazione
  void addIndividuals(const std::vector<Individual> &new_individuals) {
    if (individuals.size() + new_individuals.size() > size) {
      throw std::runtime_error(
          "Adding these individuals would exceed population size");
    }
    individuals.insert(individuals.end(), new_individuals.begin(),
                       new_individuals.end());
  }

  // rimozione di un individuo dalla popolazione
  void removeIndividual(int index) {
    if (index < 0 || index >= individuals.size()) {
      throw std::out_of_range("Index out of range");
    }
    individuals.erase(individuals.begin() + index);
  }

  // ottenere un individuo dalla popolazione
  Individual getIndividual(int index) const {
    if (index < 0 || index >= individuals.size()) {
      throw std::out_of_range("Index out of range");
    }
    return individuals[index];
  }

  // ottenere i migliori individui dalla popolazione in base al punteggio di
  // fitness
  // std::vector<Individual> getBestIndividuals(int n) const {
  //   if (n <= 0 || n > individuals.size()) {
  //     throw std::out_of_range("Invalid number of individuals requested");
  //   }
  //   std::vector<Individual> sorted_individuals = individuals;
  //   std::sort(sorted_individuals.begin(), sorted_individuals.end(),
  //             [](const Individual &a, const Individual &b) {
  //               return a.fitness > b.fitness; // Ordinamento decrescente
  //             });
  //   return std::vector<Individual>(sorted_individuals.begin(),
  //                                     sorted_individuals.begin() + n);
  // }
};
