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

template <typename T> class Individual {
  static_assert(std::is_floating_point<T>::value, "T must be float or double");

public:
  /**
   * @brief Costruttore della classe Individual
   * @param genes Vettore di intensità dei bixels
   * @param fitness Punteggio di fitness dell'individuo
   */
  std::vector<T> genes;
  int num_bixels;
  double fitness;

  Individual(double fitness, int num_bixels)
      : genes(std::vector<T>(num_bixels)), fitness(fitness),
        num_bixels(num_bixels) {}

  Individual(const std::vector<T> &genes)
      : genes(genes), fitness(0.0), num_bixels(genes.size()) {}

private:
  // metodo di modifica di un gene specifico
  void setGene(int index, T value) {
    if (index < 0 || index >= num_bixels) {
      throw std::out_of_range("Index out of range");
    }
    genes[index] = value;
  }

  // metodo di modifica del punteggio di fitness
  void setFitness(double value) { this->fitness = value; }
};
