#include <cassert>
#include <iostream>
#include <random>
#include "fmo/core/population.hpp"
#include "fmo/evaluation/evaluator.hpp"
#include "fmo/nsga2/gen-op.hpp"
#include "fmo/nsga2/offspring.hpp"
#include "fmo/preprocessing/manager.hpp"

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
    pop.addIndividual(c);
    assert(pop.size() == 3);
    std::cout << "[PASS] Test Population Dynamic Growth\n";
}

void test_evaluator_equivalence() {
    // Creazione di una matrice di influenza e ROI sintetiche
    const int total_voxels = 50;
    const int total_beamlets = 10;

    FMODataManager manager;
    // Creiamo una matrice D sparsa
    std::vector<Eigen::Triplet<double>> triplets;
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> val_dist(0.1, 5.0);
    std::uniform_int_distribution<int> row_dist(0, total_voxels - 1);

    for (int col = 0; col < total_beamlets; ++col) {
        for (int k = 0; k < 15; ++k) {
            triplets.emplace_back(row_dist(rng), col, val_dist(rng));
        }
    }

    Eigen::SparseMatrix<double> D(total_voxels, total_beamlets);
    D.setFromTriplets(triplets.begin(), triplets.end());
    D.makeCompressed();

    FMOData data_raw;
    data_raw.total_voxels = total_voxels;
    data_raw.total_beamlets = total_beamlets;
    data_raw.D = D;
    data_raw.ptv_indices = {1, 3, 5, 7, 11, 13, 17, 19};
    data_raw.rectum_indices = {2, 4, 6, 8, 10};
    data_raw.bladder_indices = {12, 14, 16, 18, 20, 22, 24};

    // Creiamo data_preprocessed estraendo le sottomatrici
    FMOData data_pre = data_raw;
    // Usiamo il manager per estrarre le sottomatrici
    FMODataManager tm(data_raw);
    tm.extractROISubmatrices();
    const FMOData &data_preprocessed = tm.getData();

    Evaluator evaluator_legacy(data_raw);
    Evaluator evaluator_optimized(data_preprocessed);

    assert(data_raw.is_preprocessed == false);
    assert(data_preprocessed.is_preprocessed == true);

    std::uniform_real_distribution<double> gene_dist(0.0, 50.0);
    for (int trial = 0; trial < 20; ++trial) {
        std::vector<double> genes(total_beamlets);
        for (int g = 0; g < total_beamlets; ++g) {
            genes[g] = gene_dist(rng);
        }
        Individual ind(genes);

        Fitness f_legacy = evaluator_legacy.evaluate(ind);
        Fitness f_opt = evaluator_optimized.evaluate(ind);

        double diff_ptv = std::abs(f_legacy.getPTVFitness() - f_opt.getPTVFitness());
        double diff_rec = std::abs(f_legacy.getRectalFitness() - f_opt.getRectalFitness());
        double diff_bla = std::abs(f_legacy.getBladderFitness() - f_opt.getBladderFitness());

        assert(diff_ptv < 1e-10);
        assert(diff_rec < 1e-10);
        assert(diff_bla < 1e-10);
    }
    std::cout << "[PASS] Test Evaluator Equivalence (Legacy vs Optimized)\n";
}

int main() {
    std::cout << "--- ESECUZIONE TEST FUNZIONALI ---\n";
    test_clipping();
    test_crowded_compare();
    test_population_limits();
    test_evaluator_equivalence();
    std::cout << "Tutti i test unitari di base hanno avuto esito positivo.\n";
    return 0;
}