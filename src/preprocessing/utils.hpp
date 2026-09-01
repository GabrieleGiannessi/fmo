/**
 * @file utils.hpp
 * @brief Funzioni di utilità per il caricamento delle matrici sparse D (matrici
 * di influenza) e concatenarle in una matrice globale. Inoltre, dispone anche
 * di metodi per ottenere dei vettori di indici VOILIST da file .mat
 *
 */

#include <Eigen/Sparse>
#include <iostream>
#include <matio.h>
#include <string>
#include <vector>
#include "fmo.hpp"

// Funzione per caricare una singola matrice sparsa D da un file .mat
Eigen::SparseMatrix<double> loadDMatrix(const std::string &filepath,
                                        int &out_rows, int &out_cols) {
  mat_t *matfp = Mat_Open(filepath.c_str(), MAT_ACC_RDONLY);
  if (!matfp) {
    throw std::runtime_error("Errore nell'apertura del file: " + filepath);
  }

  matvar_t *matvar = Mat_VarRead(matfp, "D");
  if (!matvar || matvar->class_type != MAT_C_SPARSE) {
    Mat_Close(matfp);
    throw std::runtime_error("Variabile D sparsa non trovata in: " + filepath);
  }

  mat_sparse_t *sparseData = static_cast<mat_sparse_t *>(matvar->data);
  out_rows = matvar->dims[0];
  out_cols = matvar->dims[1];

  // Creazione della matrice sparsa Eigen a partire dal formato CSC di matio
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(sparseData->ndata);

  double *data = static_cast<double *>(sparseData->data);
  for (int col = 0; col < out_cols; ++col) {
    for (int idx = sparseData->jc[col]; idx < sparseData->jc[col + 1]; ++idx) {
      int row = sparseData->ir[idx];
      double val = data[idx];
      triplets.emplace_back(row, col, val);
    }
  }

  Eigen::SparseMatrix<double> localD(out_rows, out_cols);
  localD.setFromTriplets(triplets.begin(), triplets.end());
  localD.makeCompressed();

  Mat_VarFree(matvar);
  Mat_Close(matfp);
  return localD;
}

// Funzione per caricare i vettori di indici (VOILIST)
std::vector<int> loadVoiList(const std::string &filepath,
                             const std::string &var_name) {
  mat_t *matfp = Mat_Open(filepath.c_str(), MAT_ACC_RDONLY);
  if (!matfp) {
    throw std::runtime_error("Errore nell'apertura di: " + filepath);
  }

  matvar_t *matvar = Mat_VarRead(matfp, var_name.c_str());
  if (!matvar) {
    Mat_Close(matfp);
    throw std::runtime_error("Variabile " + var_name + " non trovata.");
  }

  int num_elements = matvar->dims[0] * matvar->dims[1];
  std::vector<int> indices(num_elements);

  if (matvar->data_type == MAT_T_DOUBLE) {
    double *data = static_cast<double *>(matvar->data);
    for (int i = 0; i < num_elements; ++i) {
      // Conversione da 1-based (MATLAB) a 0-based (C++)
      indices[i] = static_cast<int>(data[i]) - 1;
    }
  } else if (matvar->data_type == MAT_T_INT32) {
    int32_t *data = static_cast<int32_t *>(matvar->data);
    for (int i = 0; i < num_elements; ++i) {
      indices[i] = data[i] - 1;
    }
  }

  Mat_VarFree(matvar);
  Mat_Close(matfp);
  return indices;
}

// Funzione per concatenare orizzontalmente le matrici D di più fasci
// D_global = [D^(1) | D^(2) | ... | D^(K)]
// Ritorna la matrice globale [M x N] dove N = somma di bixels dei fasci (i N_k)
Eigen::SparseMatrix<double>
loadGlobalDMatrix(const std::vector<std::string> &beam_filepaths) {
  if (beam_filepaths.empty()) {
    throw std::runtime_error("Lista di fasci vuota");
  }

  std::vector<Eigen::SparseMatrix<double>> D_matrices;
  int M = -1;      // numero di voxel (deve essere uguale per tutti i fasci)
  int total_N = 0; // numero totale di beamlet

  std::cout << "Caricamento matrici D per " << beam_filepaths.size()
            << " fascio/i..." << std::endl;

  // Carica tutte le matrici D
  for (size_t k = 0; k < beam_filepaths.size(); ++k) {
    int rows, cols;
    Eigen::SparseMatrix<double> D = loadDMatrix(beam_filepaths[k], rows, cols);

    if (M == -1) {
      M = rows; // Il primo file definisce il numero di voxel
    } else if (M != rows) {
      throw std::runtime_error(
          "Numero di voxel non coerente tra i fasci: " + std::to_string(rows) +
          " vs " + std::to_string(M));
    }

    D_matrices.push_back(D);
    total_N += cols;
    std::cout << "  Fascio " << (k + 1) << ": " << rows << " x " << cols
              << " (nnz: " << D.nonZeros() << ")" << std::endl;
  }

  // Concatena le matrici orizzontalmente usando Triplet
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(D_matrices[0].nonZeros() *
                   beam_filepaths.size()); // stima della capacità

  int col_offset = 0;
  for (const auto &D : D_matrices) {
    for (int col = 0; col < D.cols(); ++col) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(D, col); it; ++it) {
        triplets.emplace_back(it.row(), col_offset + col, it.value());
      }
    }
    col_offset += D.cols();
  }

  Eigen::SparseMatrix<double> D_global(M, total_N);
  D_global.setFromTriplets(triplets.begin(), triplets.end());
  D_global.makeCompressed();

  std::cout << "Matrice globale concatenata: " << M << " x " << total_N
            << " (nnz: " << D_global.nonZeros() << ")" << std::endl;

  return D_global;
}

// Funzione per caricare la matrice globale da angoli di gantry
// Accetta una directory base e un vettore di angoli (es. [0, 72, 144, 216,
// 288])
Eigen::SparseMatrix<double>
loadGlobalDMatrixFromAngles(const std::string &base_dir,
                            const std::vector<int> &gantry_angles) {
  std::vector<std::string> filepaths;

  for (int angle : gantry_angles) {
    // Formato: {base_dir}/Gantry{angle}_Couch0_D.mat
    std::string filepath =
        base_dir + "/Gantry" + std::to_string(angle) + "_Couch0_D.mat";
    filepaths.push_back(filepath);
  }

  return loadGlobalDMatrix(filepaths);
}