/**
 * @file nsga2-steps.hpp
 * @brief Operazioni timerizzate componibili dell'algoritmo NSGA-II.
 */
#pragma once

#include "fmo/evaluation/evaluator.hpp"
#include "fmo/core/population.hpp"

#include <random>
#include <vector>

void evaluatePopulation(Population &population, Evaluator &evaluator);

std::vector<std::vector<int>> sortPopulation(Population &population);

void assignPopulationCrowding(
    Population &population, const std::vector<std::vector<int>> &fronts);

Population generatePopulationOffspring(
    Population &population, double eta_c, double eta_m,
    std::mt19937 &rng);

Population mergePopulations(const Population &population,
                            const Population &offspring);

Population truncatePopulationByFronts(
    Population &population, const std::vector<std::vector<int>> &fronts,
    int population_size);
