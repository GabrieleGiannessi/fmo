/**
 * @file offspring.hpp
 * @brief definizione dell'algoritmo funzionale per la generazione della prole a
 * partire da una popolazione di individui.
 * @details La funzione generateOffspring prende in input una popolazione di
 * individui e genera una nuova popolazione di figli utilizzando le operazioni
 * genetiche di crossover e mutazione presenti in GeneticOperator (gen-op.hpp).
 * La funzione seleziona casualmente due genitori dalla popolazione, applica
 * l'operatore di crossover per generare due figli e successivamente applica
 * l'operatore di mutazione su ciascun figlio. La nuova popolazione di figli
 * viene restituita come output della funzione.
 */
#include <vector>
#include "individual.hpp"
#include "gen-op.hpp"

