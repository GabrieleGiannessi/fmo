/**
 * @file test-metrics.cpp
 * @brief Unit tests per le metriche di prestazione e stabilità (RMSE, Spread)
 *        e per la serializzazione delle popolazioni.
 */

#include "fmo/metrics/metrics.hpp"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>

void test_population_save_load() {
  Population pop;
  Individual ind1({1.5, 2.5, 3.5}, Fitness(10.0, 20.0, 30.0));
  ind1.setRank(1);
  ind1.setCrowdingDistance(0.75);

  Individual ind2({4.0, 5.0, 6.0}, Fitness(12.0, 22.0, 32.0));
  ind2.setRank(2);
  ind2.setCrowdingDistance(1.5);

  pop.addIndividual(ind1);
  pop.addIndividual(ind2);

  const std::string tmp_file = "/tmp/test_pop_save_load.csv";
  savePopulation(pop, tmp_file);

  Population loaded = loadPopulation(tmp_file);
  std::remove(tmp_file.c_str());

  assert(loaded.size() == 2);
  assert(loaded.getIndividual(0).getRank() == 1);
  assert(std::abs(loaded.getIndividual(0).getCrowdingDistance() - 0.75) < 1e-9);
  assert(std::abs(loaded.getIndividual(0).getFitness().getPTVFitness() - 10.0) < 1e-9);
  assert(std::abs(loaded.getIndividual(0).getFitness().getRectalFitness() - 20.0) < 1e-9);
  assert(std::abs(loaded.getIndividual(0).getFitness().getBladderFitness() - 30.0) < 1e-9);
  assert(loaded.getIndividual(0).genes.size() == 3);
  assert(std::abs(loaded.getIndividual(0).genes[0] - 1.5) < 1e-9);
  assert(std::abs(loaded.getIndividual(0).genes[1] - 2.5) < 1e-9);
  assert(std::abs(loaded.getIndividual(0).genes[2] - 3.5) < 1e-9);

  std::cout << "[PASS] Test Population Save & Load\n";
}

void test_rmse_identical_populations() {
  Population pop1;
  pop1.addIndividual(Individual({1.0, 2.0}, Fitness(5.0, 10.0, 15.0)));
  pop1.addIndividual(Individual({3.0, 4.0}, Fitness(7.0, 12.0, 17.0)));

  // Popolazione identica
  Population pop2 = pop1;

  PopulationDiscrepancy res = computePopulationRMSE(pop1, pop2);

  assert(std::abs(res.rmse_ptv) < 1e-12);
  assert(std::abs(res.rmse_rectum) < 1e-12);
  assert(std::abs(res.rmse_bladder) < 1e-12);
  assert(std::abs(res.rmse_total_fitness) < 1e-12);
  assert(std::abs(res.rmse_genes) < 1e-12);

  std::cout << "[PASS] Test RMSE Identical Populations (Zero Error)\n";
}

void test_rmse_perturbed_populations() {
  Population pop1;
  pop1.addIndividual(Individual({0.0, 0.0}, Fitness(10.0, 10.0, 10.0)));
  pop1.addIndividual(Individual({0.0, 0.0}, Fitness(20.0, 20.0, 20.0)));

  Population pop2;
  // Differenza di +2 su PTV, -2 su Rectum, 0 su Bladder, +1 su tutti i geni (2 individui x 2 geni)
  pop2.addIndividual(Individual({1.0, 1.0}, Fitness(12.0, 8.0, 10.0)));
  pop2.addIndividual(Individual({1.0, 1.0}, Fitness(22.0, 18.0, 20.0)));

  PopulationDiscrepancy res = computePopulationRMSE(pop1, pop2);

  // sqrt((2^2 + 2^2) / 2) = 2.0
  assert(std::abs(res.rmse_ptv - 2.0) < 1e-9);
  assert(std::abs(res.rmse_rectum - 2.0) < 1e-9);
  assert(std::abs(res.rmse_bladder - 0.0) < 1e-9);
  // total = sqrt((8 + 8 + 0) / 6) = sqrt(16 / 6) = sqrt(8 / 3)
  assert(std::abs(res.rmse_total_fitness - std::sqrt(16.0 / 6.0)) < 1e-9);
  // genes: 4 scarti di 1.0 -> sqrt(4 * 1^2 / 4) = 1.0
  assert(std::abs(res.rmse_genes - 1.0) < 1e-9);

  std::cout << "[PASS] Test RMSE Perturbed Populations\n";
}

void test_spread_metric() {
  Population pop;
  // Front 1 con punti perfettamente equi-spaziati: (1,1,1), (2,2,2), (3,3,3)
  // distanze: sqrt(3), sqrt(3). Tutte uguali -> spread = 0
  Individual i1({0.0}, Fitness(1.0, 1.0, 1.0)); i1.setRank(1);
  Individual i2({0.0}, Fitness(2.0, 2.0, 2.0)); i2.setRank(1);
  Individual i3({0.0}, Fitness(3.0, 3.0, 3.0)); i3.setRank(1);

  pop.addIndividual(i1);
  pop.addIndividual(i2);
  pop.addIndividual(i3);

  double spread = computeSpreadMetric(pop);
  assert(std::abs(spread) < 1e-9);

  std::cout << "[PASS] Test Spread Metric (Uniform Frontier = 0)\n";
}

int main() {
  std::cout << "--- ESECUZIONE TEST METRICHE E STABILITÀ ---\n";
  test_population_save_load();
  test_rmse_identical_populations();
  test_rmse_perturbed_populations();
  test_spread_metric();
  std::cout << "Tutti i test delle metriche sono passati con successo!\n";
  return 0;
}
