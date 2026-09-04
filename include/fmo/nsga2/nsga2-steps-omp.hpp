/**
 * @file nsga2-steps-omp.hpp
 * @brief file che espone le firme delle procedure parallelizzate utilizzate
 * nella versione parallelizzata implementata attraverso OpenMP dell'algoritmo
 * NSGA-II.
 */

#pragma once

#include "fmo/core/population.hpp"
#include "fmo/evaluation/evaluator.hpp"

#include <random>
#include <vector>

// calcolo delle fitness in parallelo (embarassingly parallel)
void evaluatePopulationOmp(Population &population, Evaluator &evaluator);

// generazione della prole parallela
Population generatePopulationOffspringOmp(Population &population,
                                          double crossover_probability,
                                          double mutation_probability,
                                          double eta_c, double eta_m,
                                          std::mt19937 &rng);

// assegnazione crowding distance parallelo
void assignPopulationCrowdingOmp(Population &population,
                                 const std::vector<std::vector<int>> &fronts);

// non-dominated sorting in parallelo
std::vector<std::vector<int>> sortPopulationOmp(Population &population);