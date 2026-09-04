/**
 * @file fmo.hpp
 * @brief Struttura per memorizzare il problema FMO pre-elaborato
 * @details Questa struttura contiene le informazioni necessarie per
 * rappresentare il problema FMO (Fluence Map Optimization) dopo la fase di
 * pre-elaborazione. Include il numero totale di voxel, il numero totale di
 * beamlet, la matrice globale concatenata D e i vettori di indici per le
 * regioni di interesse (PTV, retto e vescica).
 */
#pragma once
#include <Eigen/Sparse>
#include <vector>

struct FMOData {
  int total_voxels = 0;             // M voxels
  int total_beamlets = 0;           // N beamlets (o bixels)
  Eigen::SparseMatrix<double> D;    // Matrice di influenza [M x N]
  std::vector<int> ptv_indices;     // Indici PTV (Planning Target Volume)
  std::vector<int> rectum_indices;  // Indici retto
  std::vector<int> bladder_indices; // Indici vescica
};
