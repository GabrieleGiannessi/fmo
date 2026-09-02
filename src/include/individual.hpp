/**
 * @file individual.hpp
 * @brief Definizione della classe Individual per rappresentare un individuo
 * nella popolazione genetica
 * @details La classe Individual rappresenta un individuo nella popolazione
 * genetica (cromosoma), contenente un vettore di intensità dei bixels e un
 * punteggio di fitness. Fornisce metodi per accedere e modificare i geni e il
 * punteggio di fitness.
 */

#include <type_traits>
#include <vector>
#include "fitness.hpp"
#include <stdexcept>

class Individual {
public:
  /**
   * @brief Costruttore della classe Individual
   * @param genes Vettore di intensità dei bixels
   * @param fitness Punteggio di fitness dell'individuo
   */
  std::vector<double> genes;
  int num_bixels;
  Fitness fitness;

  Individual(Fitness fitness, int num_bixels)
      : genes(std::vector<double>(num_bixels)), fitness(fitness),
        num_bixels(num_bixels) {}

  Individual(const std::vector<double> &genes)
      : genes(genes), fitness(Fitness()), num_bixels(genes.size()) {}

private:
  // metodo di modifica di un gene specifico
  void setGene(int index, double value) {
    if (index < 0 || index >= num_bixels) {
      throw std::out_of_range("Index out of range");
    }
    this->genes[index] = value;
  }

  // metodo di modifica del punteggio di fitness
  void setFitness(Fitness new_fitness) { this->fitness = new_fitness; }
};
