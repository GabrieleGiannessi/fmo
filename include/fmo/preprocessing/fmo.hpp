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
  int global_nonzeros = 0;          // Conteggio non-zeri prima della deallocazione
  Eigen::SparseMatrix<double> D;    // Matrice di influenza globale [M x N]
  std::vector<int> ptv_indices;     // Indici PTV (Planning Target Volume)
  std::vector<int> rectum_indices;  // Indici retto
  std::vector<int> bladder_indices; // Indici vescica

  // Sottomatrici per le singole Regioni di Interesse (ROI)
  Eigen::SparseMatrix<double> D_ptv;    // Matrice PTV [num_ptv_voxels x N]
  Eigen::SparseMatrix<double> D_rectum; // Matrice Retto [num_rectum_voxels x N]
  Eigen::SparseMatrix<double>
      D_bladder; // Matrice Vescica [num_bladder_voxels x N]

  // Vettori di pesi lineari medi precalcolati [dimensione N]
  std::vector<double> rectum_weights;
  std::vector<double> bladder_weights;

  bool is_preprocessed =
      false; // Flag per indicare se le sottomatrici sono state estratte
};
