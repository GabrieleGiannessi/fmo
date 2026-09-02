/**
 * @file gen-op.hpp
 * @brief Definizione della classe GenOp per rappresentare le operazioni
 * genetiche utilizzate nell'algoritmo genetico. La classe implementa le
 * operazioni di crossover e mutazione per generare nuovi individui a partire da
 * quelli esistenti.
 * @details La classe GenOp rappresenta le operazioni genetiche utilizzate
 * nell'algoritmo genetico per generare nuovi individui a partire da quelli
 * esistenti. La classe implementa le operazioni di crossover (SBX) e mutazione,
 * che consentono di combinare le caratteristiche di due individui genitori per
 * creare un nuovo individuo figlio, e di introdurre variazioni casuali nei geni
 * di un individuo per esplorare nuove soluzioni. La classe GenOp è progettata
 * per essere utilizzata in combinazione con la classe Individual e Population,
 * che rappresenta un individuo nella popolazione genetica.
 *
 */

#include "individual.hpp"
#include <utility>

class GeneticOperator {
public:
  /** @brief Implementazione della selezione a torneo per la scelta dei genitori
   * @param population Popolazione di individui da cui selezionare i genitori
   * @param tournament_size Dimensione del torneo (numero di individui da
   * confrontare)
   * @return Vettore di individui selezionati come genitori
   */
  static std::vector<Individual>
  tournament_selection(const std::vector<Individual> &population,
                       size_t tournament_size) {}

  /** @brief Implementazione del crossover SBX (Simulated Binary Crossover)
   * @param parent1 Primo genitore
   * @param parent2 Secondo genitore
   * @param crossover_probability Probabilità di crossover
   * @return Coppia di figli generati dai genitori
   */
  static std::pair<Individual, Individual>
  crossover(const Individual &parent1,
            const Individual &parent2, double crossover_probability) {}

  /**
   * @brief Implementazione della mutazione polinomiale dei geni di un individuo
   * @param individual Individuo da mutare
   * @param mutation_probability Probabilità di mutazione per ciascun gene
   */
  static void mutate(Individual &individual,
                     double mutation_probability) {}

    /**
     * @brief Metodo statico per limitare i valori dei geni di un individuo (le intensità devono rientrare nei valori positivi consentiti).
     * @param ind Riferimento all'individuo da modificare.
     * @param min_val Valore minimo consentito per i geni (default: 0.0).
     * @param max_val Valore massimo consentito per i geni (default: 100.0).
     */
    static void clampBounds(
        Individual& ind, 
        double min_val = 0.0, 
        double max_val = 100.0
    );
};