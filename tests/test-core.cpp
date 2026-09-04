#include <cassert>
#include <iostream>
#include <random>
#include "fmo/nsga2/gen-op.hpp"
#include "fmo/core/population.hpp"
#include "fmo/nsga2/offspring.hpp"

void test_clipping() {
    Individual ind({-10.0, 5.0, 150.0});
    GeneticOperator::clip(ind, 0.0, 100.0);
    assert(ind.genes[0] == 0.0);
    assert(ind.genes[1] == 5.0);
    assert(ind.genes[2] == 100.0);
    std::cout << "[PASS] Test Clipping\n";
}

void test_crowded_compare() {
    Individual best, worst;
    best.rank = 1;  best.crowding_distance = 0.5;
    worst.rank = 2; worst.crowding_distance = 1.0;
    assert(GeneticOperator::crowded_compare(best, worst) == true);

    Individual dense, sparse;
    dense.rank = 1;  dense.crowding_distance = 0.1;
    sparse.rank = 1; sparse.crowding_distance = 0.9;
    assert(GeneticOperator::crowded_compare(sparse, dense) == true);
    std::cout << "[PASS] Test Crowded Compare\n";
}

void test_population_limits() {
    Population pop(2);
    Individual a(5), b(5), c(5);
    pop.addIndividual(a);
    pop.addIndividual(b);
    try {
        pop.addIndividual(c);
        assert(false); // Non deve arrivare qui
    } catch (const std::runtime_error&) {
        std::cout << "[PASS] Test Population Overflow Guard\n";
    }
}

int main() {
    std::cout << "--- ESECUZIONE TEST FUNZIONALI ---\n";
    test_clipping();
    test_crowded_compare();
    test_population_limits();
    std::cout << "Tutti i test unitari di base hanno avuto esito positivo.\n";
    return 0;
}