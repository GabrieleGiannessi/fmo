/**
 * @file preprocessing/main.cpp
 * @brief Test per il caricamento delle matrici di influenza D da file .mat
 * @details Questo file contiene un test per verificare il corretto caricamento
 * delle matrici di influenza D da file .mat. Viene testato sia il caricamento
 * di una singola matrice D che la concatenazione di più matrici D in una
 * matrice globale.
 */

#include "src/preprocessing/utils.hpp"
#include <Eigen/Sparse>
#include <iostream>
#include <matio.h>
#include <stdexcept>
#include <string>
#include <vector>

int main() {
  try {
    std::cout << "\n=== Test caricamento matrici di influenza ===" << std::endl;

    // Test 1: Singola matrice D
    std::cout << "\n1. Caricamento singola matrice D (Gantry 0°)..."
              << std::endl;
    int out_rows, out_cols;
    std::string d_path = "../data/phantom/Gantry0_Couch0_D.mat";
    Eigen::SparseMatrix<double> D_single =
        loadDMatrix(d_path, out_rows, out_cols);

    std::cout << "✓ Dimensioni: " << out_rows << " x " << out_cols << std::endl;
    std::cout << "  Non-zero elements: " << D_single.nonZeros() << std::endl;

    // Test 2: Matrice globale da singoli fasci
    std::cout << "\n2. Caricamento matrice globale (configurazione standard 5 "
                 "fasci -> 0, 72, 144, 216, 288)..."
              << std::endl;
    std::vector<int> gantry_angles = {0, 72, 144, 216, 288};
    std::string base_dir = "../data/phantom";

    Eigen::SparseMatrix<double> D_global =
        loadGlobalDMatrixFromAngles(base_dir, gantry_angles);

    std::cout << "\n✓ Matrice globale caricata con successo!" << std::endl;
    std::cout << "  Dimensioni finali: " << D_global.rows() << " x "
              << D_global.cols() << std::endl;
    std::cout << "  Non-zero elements: " << D_global.nonZeros() << std::endl;
    std::cout << "  Densità: "
              << (100.0 * D_global.nonZeros() /
                  (D_global.rows() * D_global.cols()))
              << "%" << std::endl;

    // Test 3: Caricamento con liste esplicite di file
    std::cout
        << "\n3. Caricamento matrice globale da lista esplicita di file..."
        << std::endl;
    std::vector<std::string> filepaths = {
        "../data/phantom/Gantry0_Couch0_D.mat",
        "../data/phantom/Gantry72_Couch0_D.mat",
        "../data/phantom/Gantry144_Couch0_D.mat"};

    Eigen::SparseMatrix<double> D_partial = loadGlobalDMatrix(filepaths);
    std::cout << "\n✓ Matrice parziale (3 fasci) caricata!" << std::endl;
    std::cout << "  Dimensioni: " << D_partial.rows() << " x "
              << D_partial.cols() << std::endl;
    std::cout << "  Non-zero elements: " << D_partial.nonZeros() << std::endl;

    return 0;
  } catch (const std::exception &e) {
    std::cerr << "Errore critico: " << e.what() << std::endl;
    return 1;
  }
}