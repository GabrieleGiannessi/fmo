/**
 * @file tests/load-matrix.cpp
 * @brief Test per il caricamento delle matrici di influenza D da file .mat
 * @details Questo file contiene un test per verificare il corretto caricamento
 * delle matrici di influenza D da file .mat usando la classe FMODataManager.
 * Dimostra l'uso di un manager che mantiene internamente uno stato FMOData
 * e lo popola attraverso i metodi di caricamento.
 */

#include "src/preprocessing/manager.hpp"
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

    std::cout << "\n✓ Test completati con successo!" << std::endl;
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "Errore critico: " << e.what() << std::endl;
    return 1;
  }
}