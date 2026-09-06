/**
 * @file nsga2-steps-ff.hpp
 * @brief file che espone le firme delle procedure parallelizzate utilizzate
 * nella loro versione FastFlow per l'implementazione dell'algoritmo NSGA2.
 */

#pragma once

#include "fmo/core/population.hpp"
#include "fmo/evaluation/evaluator.hpp"

#include <random>
#include <vector>

void evaluatePopulationFFFarm(Population &population, Evaluator &evaluator,
                              int nw);

void evaluatePopulationFFParFor(Population &population, Evaluator &evaluator,
                                int nw); 