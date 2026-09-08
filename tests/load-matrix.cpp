/**
 * @file tests/load-matrix.cpp
 * @brief Test per il caricamento delle matrici di influenza D da file .mat
 * @details Questo file contiene un test per verificare il corretto caricamento
 * delle matrici di influenza D da file .mat usando la classe FMODataManager.
 * Dimostra l'uso di un manager che mantiene internamente uno stato FMOData
 * e lo popola attraverso i metodi di caricamento.
 */

#include "fmo/preprocessing/manager.hpp"
#include <iostream>
#include <stdexcept>
#include <vector>

int main() {
  try {
    std::cout << "\n=== Test FMODataManager ===" << std::endl;

    // Test 1: Caricamento singola matrice D
    std::cout << "\n1. Creazione manager e caricamento singola matrice D..."
              << std::endl;
    FMODataManager manager1;
    manager1.loadDMatrixFromFile("../data/phantom/Gantry0_Couch0_D.mat");
    manager1.printSummary();

    // Test 2: Caricamento matrice globale da angoli
    std::cout << "\n2. Caricamento matrice globale da 5 fasci (0°, 72°, 144°, "
                 "216°, 288°)..."
              << std::endl;
    FMODataManager manager2;
    std::vector<int> gantry_angles = {0, 72, 144, 216, 288};
    manager2.loadGlobalDMatrixFromAngles("../data/phantom", gantry_angles);
    manager2.printSummary();

    // Test 3: Caricamento da lista esplicita di file
    std::cout << "\n3. Caricamento matrice globale da lista esplicita di 3 "
                 "fasci..."
              << std::endl;
    FMODataManager manager3;
    std::vector<std::string> filepaths = {
        "../data/phantom/Gantry0_Couch0_D.mat",
        "../data/phantom/Gantry72_Couch0_D.mat",
        "../data/phantom/Gantry144_Couch0_D.mat"};
    manager3.loadGlobalDMatrixFromFiles(filepaths);
    manager3.printSummary();

    // Test 4: Accesso ai dati tramite getter
    std::cout << "\n4. Test getter per accedere ai dati..." << std::endl;
    std::cout << "Dimensioni matrice: " << manager2.getTotalVoxels() << " x "
              << manager2.getTotalBeamlets() << std::endl;
    const auto &D_global = manager2.getInfluenceMatrix();
    std::cout << "Non-zero elements: " << D_global.nonZeros() << std::endl;
    std::cout << "Densità: "
              << (100.0 * D_global.nonZeros() /
                  (D_global.rows() * D_global.cols()))
              << "%" << std::endl;

    // Test 5: Accesso alla struttura completa FMOData
    std::cout << "\n5. Accesso alla struttura FMOData completa..." << std::endl;
    const FMOData &fmo_complete = manager2.getData();
    std::cout << "Total voxels (dal getData): " << fmo_complete.total_voxels
              << std::endl;
    std::cout << "Total beamlets (dal getData): "
              << fmo_complete.total_beamlets << std::endl;

    // Test 6: Caricamento VOILIST
    std::cout << "\n6. Caricamento VOILIST (ROI indices)..." << std::endl;
    manager2.loadVoiList("../data/phantom/BODY_VOILIST.mat", "v", "ptv");
    manager2.loadVoiList("../data/phantom/Core_VOILIST.mat", "v", "rectum");
    manager2.loadVoiList("../data/phantom/OuterTarget_VOILIST.mat", "v",
                         "bladder");
    manager2.printSummary();

    // Test 7: Estrazione sottomatrici ROI e verifica coerenza
    std::cout << "\n7. Estrazione sottomatrici ROI..." << std::endl;
    manager2.extractROISubmatrices();
    manager2.printSummary();

    if (!manager2.isPreprocessed()) {
      throw std::runtime_error("Errore: isPreprocessed() dovrebbe essere true.");
    }
    if (manager2.getDPTV().rows() != static_cast<int>(manager2.getPTVIndices().size())) {
      throw std::runtime_error("Errore: dimensione righe D_ptv non coerente.");
    }
    if (manager2.getDRectum().rows() != static_cast<int>(manager2.getRectumIndices().size())) {
      throw std::runtime_error("Errore: dimensione righe D_rectum non coerente.");
    }
    if (manager2.getDBladder().rows() != static_cast<int>(manager2.getBladderIndices().size())) {
      throw std::runtime_error("Errore: dimensione righe D_bladder non coerente.");
    }

    // Verifica coerenza tra D_rectum e rectum_weights
    const auto &D_rec = manager2.getDRectum();
    const auto &w_rec = manager2.getRectumWeights();
    double norm_rec = static_cast<double>(manager2.getRectumIndices().size());
    for (int col = 0; col < D_rec.cols(); ++col) {
      double sum_col = 0.0;
      for (Eigen::SparseMatrix<double>::InnerIterator it(D_rec, col); it; ++it) {
        sum_col += it.value();
      }
      double expected = sum_col / norm_rec;
      if (std::abs(expected - w_rec[col]) > 1e-12) {
        throw std::runtime_error("Errore: discrepanza tra D_rectum e rectum_weights alla colonna " + std::to_string(col));
      }
    }
    std::cout << "✓ Coerenza matematica tra sottomatrici e pesi verificata con successo!" << std::endl;

    std::cout << "\n✓ Test completati con successo!" << std::endl;
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "Errore critico: " << e.what() << std::endl;
    return 1;
  }
}