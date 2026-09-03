/**
 * @file individual.hpp
 * @brief Definizione della classe Individual per rappresentare un individuo
 * nella popolazione genetica
 * @details La classe Individual rappresenta un individuo nella popolazione
 * genetica (cromosoma), contenente un vettore di intensità dei bixels e un
 * punteggio di fitness. Fornisce metodi per accedere e modificare i geni e il
 * punteggio di fitness.
 */
#pragma once
#include "fitness.hpp"
#include <stdexcept>
#include <type_traits>
#include <vector>

class Individual {
public:
  /**
   * @brief Costruttore della classe Individual
   * @param genes Vettore di intensità dei bixels
   * @param fitness Punteggio di fitness dell'individuo
   * @param rank Rango di dominanza di Pareto (1 = ottimo)
   * @param crowding_distance Distanza di affollamento (stima della densità
   * attorno alla soluzione)
   */
  std::vector<double> genes;
  Fitness fitness;
  int rank = 0;
  double crowding_distance = 0.0;

  // Costruttore di default
  Individual() = default;

  /**
   * @brief Costruttore con dimensione prefissata
   * @param num_bixels Numero totale di bixels del piano di trattamento
   */
  explicit Individual(int num_bixels)
      : genes(num_bixels, 0.0), fitness(), rank(0), crowding_distance(0.0) {}

  /**
   * @brief Costruttore con geni preesistenti
   * @param genes Vettore di intensità dei bixels
   */
  explicit Individual(const std::vector<double> &genes)
      : genes(genes), fitness(), rank(0), crowding_distance(0.0) {}

  /**
   * @brief Costruttore completo
   */
  Individual(const std::vector<double> &genes, const Fitness &fitness)
      : genes(genes), fitness(fitness), rank(0), crowding_distance(0.0) {}

  // Utility per ottenere la dimensione del cromosoma
  inline size_t size() const { return genes.size(); }

  // metodo di modifica di un gene specifico
  inline void setGene(int index, double value) {
    if (index >= genes.size() || index < 0) {
      throw std::out_of_range("Index out of range");
    }
    this->genes[index] = value;
  }

  // Metodo di accesso al gene specifico
  inline double getGene(size_t index) const {
    if (index >= genes.size()) {
      throw std::out_of_range("Index out of range in getGene");
    }
    return genes[index];
  }

  // metodo di accesso alle fitness dell'individuo
  inline Fitness getFitness() const { return fitness; }

  // metodo di modifica del punteggio di fitness
  inline void setFitness(Fitness new_fitness) { this->fitness = new_fitness; }
  
  inline int getRank() const { return rank; }
  inline void setRank(int new_rank) { this->rank = new_rank; }
  inline double getCrowdingDistance() const { return crowding_distance; }
  inline void setCrowdingDistance(double new_distance) {
    this->crowding_distance = new_distance;
  }
};
