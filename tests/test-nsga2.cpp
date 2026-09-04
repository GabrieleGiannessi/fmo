#include "fmo/nsga2/gen-op.hpp"
#include "fmo/nsga2/nsga-utils.hpp"
#include "fmo/nsga2/offspring.hpp"
#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

namespace {
int failures = 0;
void check(bool condition, const std::string &message) {
  if (!condition) { ++failures; std::cerr << "[FAIL] " << message << '\n'; }
}
void check_throws(const std::function<void()> &operation, const std::string &message) {
  try { operation(); check(false, message); }
  catch (const std::exception &) { check(true, message); }
}
Population make_population(const std::vector<Fitness> &fitnesses) {
  Population population(fitnesses.size());
  for (const Fitness &fitness : fitnesses)
    population.addIndividual(Individual({0.0}, fitness));
  return population;
}

void test_dominance() {
  const Fitness best(1.0, 2.0, 3.0);
  check(dominates(best, Fitness(2.0, 2.0, 3.0)), "dominance accepts equal objectives");
  check(!dominates(best, best), "an individual does not dominate itself");
  check(!dominates(Fitness(1.0, 4.0, 3.0), Fitness(2.0, 2.0, 3.0)), "incomparable points do not dominate");
  check(!dominates(Fitness(2.0, 3.0, 4.0), best), "a worse point does not dominate");
}

void test_fast_non_dominated_sort() {
  Population population = make_population({
      Fitness(1, 3, 3), Fitness(3, 1, 3), Fitness(3, 3, 1),
      Fitness(3, 3, 3), Fitness(4, 4, 4)});
  const auto fronts = fastNondominatedSort(population);
  check(fronts.size() == 3, "sorting creates all expected fronts");
  check(fronts[0] == std::vector<int>({0, 1, 2}), "trade-off points are in front one");
  check(fronts[1] == std::vector<int>{3}, "first dominated point is in front two");
  check(fronts[2] == std::vector<int>{4}, "second dominated point is in front three");
  check(population.getIndividual(0).getRank() == 1 && population.getIndividual(4).getRank() == 3,
        "sorting writes Pareto ranks");
}

void test_crowding_distance() {
  Population population = make_population({Fitness(0, 3, 3), Fitness(1, 2, 2), Fitness(2, 1, 1)});
  assignCrowdingDistance({0, 1, 2}, population);
  check(std::isinf(population.getIndividual(0).getCrowdingDistance()), "crowding keeps first boundary");
  check(std::isinf(population.getIndividual(2).getCrowdingDistance()), "crowding keeps second boundary");
  check(std::abs(population.getIndividual(1).getCrowdingDistance() - 3.0) < 1e-12,
        "crowding accumulates normalized distances");
  Population small = make_population({Fitness(1, 1, 1), Fitness(2, 2, 2)});
  assignCrowdingDistance({0, 1}, small);
  check(std::isinf(small.getIndividual(0).getCrowdingDistance()) &&
            std::isinf(small.getIndividual(1).getCrowdingDistance()),
        "fronts of size two get infinite crowding");
}

void test_truncation_and_elitism() {
  Population combined = make_population({Fitness(0, 0, 0), Fitness(1, 1, 1), Fitness(2, 2, 2), Fitness(3, 3, 3)});
  combined.getIndividual(0).setRank(1); combined.getIndividual(1).setRank(1);
  combined.getIndividual(2).setRank(2); combined.getIndividual(3).setRank(2);
  combined.getIndividual(2).setCrowdingDistance(0.1); combined.getIndividual(3).setCrowdingDistance(0.9);
  Population next = truncatePopulation(combined, {{0, 1}, {2, 3}}, 3);
  check(next.size() == 3, "truncation returns requested size");
  check(next.getIndividual(0).getRank() == 1 && next.getIndividual(1).getRank() == 1,
        "elitism keeps complete better fronts");
  check(next.getIndividual(2).getFitness().value_target_ptv == 3.0,
        "truncation keeps least crowded split-front member");
}

void test_genetic_operators() {
  Individual parent1({-10.0, 20.0, 150.0});
  Individual parent2({10.0, 40.0, 50.0});
  std::mt19937 rng(7);
  auto unchanged = GeneticOperator::crossover(parent1, parent2, 0.0, 20.0, rng);
  check(unchanged.first.genes == parent1.genes && unchanged.second.genes == parent2.genes,
        "disabled crossover copies parents");
  check_throws([&] { GeneticOperator::crossover(Individual(std::vector<double>{1.0}), parent2, 1.0, 20.0, rng); },
               "crossover rejects different chromosome sizes");
  auto children = GeneticOperator::crossover(parent1, parent2, 1.0, 20.0, rng);
  for (const Individual &child : {children.first, children.second})
    for (double gene : child.genes)
      check(gene >= 0.0 && gene <= 100.0, "crossover clips child genes");
  Individual individual({10.0, 20.0, 30.0});
  GeneticOperator::mutate(individual, 0.0, 20.0, rng);
  check(individual.genes == std::vector<double>({10.0, 20.0, 30.0}), "disabled mutation preserves genes");
  Individual clipped({-5.0, 50.0, 120.0});
  GeneticOperator::mutate(clipped, 0.0, 20.0, rng, 0.0, 100.0);
  check(clipped.genes == std::vector<double>({0.0, 50.0, 100.0}), "mutation applies clipping");
  Individual preferred(std::vector<double>{0.0}); preferred.setRank(1); preferred.setCrowdingDistance(0.1);
  Individual fallback(std::vector<double>{0.0}); fallback.setRank(2); fallback.setCrowdingDistance(100.0);
  check(GeneticOperator::crowded_compare(preferred, fallback), "comparison prioritizes rank");
  check(!GeneticOperator::crowded_compare(fallback, preferred), "crowding cannot overcome worse rank");
}

void test_random_population_and_offspring() {
  std::mt19937 first_rng(42), second_rng(42);
  Population first = generateRandomPopulation(8, 4, first_rng);
  Population second = generateRandomPopulation(8, 4, second_rng);
  check(first.size() == 8 && first.getIndividual(0).size() == 4, "random population has requested dimensions");
  check(first.getIndividual(0).genes == second.getIndividual(0).genes, "random population is reproducible");
  for (const Individual &individual : first.individuals)
    for (double gene : individual.genes)
      check(gene >= 0.0 && gene <= 10.0, "initial genes respect range");
  std::mt19937 offspring_rng(9);
  Population offspring = generateOffSpring(first, 0.0, 0.0, 20.0, 20.0, offspring_rng);
  check(offspring.size() == first.size(), "offspring preserves even population size");
  for (const Individual &child : offspring.individuals) {
    check(child.size() == 4, "offspring preserves chromosome size");
    for (double gene : child.genes)
      check(gene >= 0.0 && gene <= 100.0, "offspring genes stay in bounds");
  }
}
}

int main() {
  test_dominance(); test_fast_non_dominated_sort(); test_crowding_distance();
  test_truncation_and_elitism(); test_genetic_operators(); test_random_population_and_offspring();
  if (failures != 0) { std::cerr << failures << " NSGA-II test(s) failed.\n"; return 1; }
  std::cout << "All NSGA-II tests passed.\n";
  return 0;
}