/**
 * @file utils.hpp
 * @brief Classe FMODataManager per il caricamento e la manipolazione dei dati FMO
 * @details La classe gestisce internamente un'istanza di FMOData e fornisce metodi
 * per caricare e manipolare le matrici sparse D (matrici di influenza) e i vettori
 * di indici VOILIST (Volume of Interest List) da file .mat. Gli indici VOILIST 
 * sono utilizzati per identificare le regioni di interesse (ROI) come PTV, retto e vescica.
 *
 */

#ifndef FMO_UTILS_HPP
#define FMO_UTILS_HPP

#include <Eigen/Sparse>
#include <iostream>
#include <matio.h>
#include <string>
#include <vector>
#include "fmo.hpp"

/**
 * @class FMODataManager
 * @brief Manager per il caricamento e la manipolazione dei dati FMO
 * @details Mantiene internamente un'istanza di FMOData e fornisce metodi per:
 *   - Caricamento di matrici D singole e multiple da file .mat
 *   - Concatenazione orizzontale di matrici D da multiple configurazioni di fasci
 *   - Caricamento di vettori di indici VOILIST (regioni di interesse)
 *   - Accesso ai dati FMO elaborati
 */
class FMODataManager {
private:
  // Stato interno
  FMOData fmo_data;

  /**
   * @brief Helper privato per leggere una matrice sparsa D da file matio
   * @param matvar Puntatore alla variabile MAT letto da matio
   * @param out_rows Numero di righe della matrice
   * @param out_cols Numero di colonne della matrice
   * @return Matrice sparsa Eigen convertita dal formato CSC di matio
   */
  static Eigen::SparseMatrix<double>
  readSparseMatrixFromMatio(matvar_t *matvar, int &out_rows, int &out_cols) {
    mat_sparse_t *sparseData = static_cast<mat_sparse_t *>(matvar->data);
    out_rows = matvar->dims[0];
    out_cols = matvar->dims[1];

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(sparseData->ndata);

    double *data = static_cast<double *>(sparseData->data);
    for (int col = 0; col < out_cols; ++col) {
      for (int idx = sparseData->jc[col]; idx < sparseData->jc[col + 1];
           ++idx) {
        int row = sparseData->ir[idx];
        double val = data[idx];
        triplets.emplace_back(row, col, val);
      }
    }

    Eigen::SparseMatrix<double> matrix(out_rows, out_cols);
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    matrix.makeCompressed();
    return matrix;
  }

  /**
   * @brief Helper privato per caricare una singola matrice D senza aggiornare
   * lo stato interno
   */
  static Eigen::SparseMatrix<double> loadDMatrixHelper(const std::string &filepath,
                                                       int &out_rows,
                                                       int &out_cols) {
    mat_t *matfp = Mat_Open(filepath.c_str(), MAT_ACC_RDONLY);
    if (!matfp) {
      throw std::runtime_error("Errore nell'apertura del file: " + filepath);
    }

    matvar_t *matvar = Mat_VarRead(matfp, "D");
    if (!matvar || matvar->class_type != MAT_C_SPARSE) {
      Mat_Close(matfp);
      throw std::runtime_error("Variabile D sparsa non trovata in: " +
                               filepath);
    }

    Eigen::SparseMatrix<double> result =
        readSparseMatrixFromMatio(matvar, out_rows, out_cols);

    Mat_VarFree(matvar);
    Mat_Close(matfp);
    return result;
  }

public:
  /**
   * @brief Costruttore di default
   */
  FMODataManager() : fmo_data() {}

  /**
   * @brief Carica una singola matrice sparsa D da file .mat nello stato interno
   * @param filepath Percorso completo del file .mat
   * @throw std::runtime_error Se il file non può essere aperto o la variabile
   * D non è trovata
   */
  void loadDMatrixFromFile(const std::string &filepath) {
    int rows, cols;
    Eigen::SparseMatrix<double> D = loadDMatrixHelper(filepath, rows, cols);

    fmo_data.D = D;
    fmo_data.total_voxels = rows;
    fmo_data.total_beamlets = cols;

    std::cout << "✓ Matrice D caricata: " << rows << " x " << cols
              << " (nnz: " << D.nonZeros() << ")" << std::endl;
  }

  /**
   * @brief Carica un vettore di indici VOILIST da file .mat nello stato interno
   * @param filepath Percorso del file .mat
   * @param var_name Nome della variabile nel file
   * @param roi_type Tipo di ROI: "ptv", "rectum", "bladder"
   * @throw std::runtime_error Se il file non può essere aperto o la variabile
   * non è trovata
   */
  void loadVoiList(const std::string &filepath, const std::string &var_name,
                   const std::string &roi_type) {
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

    // Assegna al tipo di ROI appropriato
    if (roi_type == "ptv") {
      fmo_data.ptv_indices = indices;
      std::cout << "✓ PTV_VOILIST caricato: " << indices.size() << " voxel"
                << std::endl;
    } else if (roi_type == "rectum") {
      fmo_data.rectum_indices = indices;
      std::cout << "✓ RECTUM_VOILIST caricato: " << indices.size() << " voxel"
                << std::endl;
    } else if (roi_type == "bladder") {
      fmo_data.bladder_indices = indices;
      std::cout << "✓ BLADDER_VOILIST caricato: " << indices.size() << " voxel"
                << std::endl;
    } else {
      throw std::runtime_error("Tipo di ROI non riconosciuto: " + roi_type);
    }
  }

  /**
   * @brief Concatena orizzontalmente matrici D di più fasci nello stato interno
   * @details Crea la matrice globale D = [D^(1) | D^(2) | ... | D^(K)]
   * @param beam_filepaths Vettore di percorsi ai file .mat contenenti le
   * matrici D
   * @throw std::runtime_error Se la lista è vuota o le matrici hanno righe
   * incoerenti
   */
  void loadGlobalDMatrixFromFiles(const std::vector<std::string> &beam_filepaths) {
    if (beam_filepaths.empty()) {
      throw std::runtime_error("Lista di fasci vuota");
    }

    std::vector<Eigen::SparseMatrix<double>> D_matrices;
    int M = -1;
    int total_N = 0;

    std::cout << "Caricamento matrici D per " << beam_filepaths.size()
              << " fascio/i..." << std::endl;

    // Carica tutte le matrici D
    for (size_t k = 0; k < beam_filepaths.size(); ++k) {
      int rows, cols;
      Eigen::SparseMatrix<double> D =
          loadDMatrixHelper(beam_filepaths[k], rows, cols);

      if (M == -1) {
        M = rows;
      } else if (M != rows) {
        throw std::runtime_error("Numero di voxel non coerente tra i fasci: " +
                                 std::to_string(rows) + " vs " +
                                 std::to_string(M));
      }

      D_matrices.push_back(D);
      total_N += cols;
      std::cout << "  Fascio " << (k + 1) << ": " << rows << " x " << cols
                << " (nnz: " << D.nonZeros() << ")" << std::endl;
    }

    // Concatena le matrici orizzontalmente
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(D_matrices[0].nonZeros() * beam_filepaths.size());

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

    // Aggiorna lo stato interno
    fmo_data.D = D_global;
    fmo_data.total_voxels = M;
    fmo_data.total_beamlets = total_N;

    std::cout << "Matrice globale concatenata: " << M << " x " << total_N
              << " (nnz: " << D_global.nonZeros() << ")" << std::endl;
  }

  /**
   * @brief Carica matrice globale partendo da una lista di angoli di gantry
   * @details Costruisce automaticamente i percorsi ai file basandosi sulla
   * convenzione di naming: {base_dir}/Gantry{angle}_Couch0_D.mat
   * @param base_dir Directory base contenente i file (es.
   * "../data/phantom")
   * @param gantry_angles Vettore di angoli in gradi (es. [0, 72, 144, 216,
   * 288])
   * @throw std::runtime_error Se uno dei file non può essere caricato
   */
  void loadGlobalDMatrixFromAngles(const std::string &base_dir,
                                    const std::vector<int> &gantry_angles) {
    std::vector<std::string> filepaths;

    for (int angle : gantry_angles) {
      std::string filepath =
          base_dir + "/Gantry" + std::to_string(angle) + "_Couch0_D.mat";
      filepaths.push_back(filepath);
    }

    loadGlobalDMatrixFromFiles(filepaths);
  }

  // ========================================================================
  // Getter per accedere ai dati
  // ========================================================================

  /**
   * @brief Ritorna il numero totale di voxel
   */
  int getTotalVoxels() const { return fmo_data.total_voxels; }

  /**
   * @brief Ritorna il numero totale di beamlet
   */
  int getTotalBeamlets() const { return fmo_data.total_beamlets; }

  /**
   * @brief Ritorna la matrice di influenza globale
   */
  const Eigen::SparseMatrix<double> &getInfluenceMatrix() const {
    return fmo_data.D;
  }

  /**
   * @brief Ritorna gli indici PTV
   */
  const std::vector<int> &getPTVIndices() const {
    return fmo_data.ptv_indices;
  }

  /**
   * @brief Ritorna gli indici retto
   */
  const std::vector<int> &getRectumIndices() const {
    return fmo_data.rectum_indices;
  }

  /**
   * @brief Ritorna gli indici vescica
   */
  const std::vector<int> &getBladderIndices() const {
    return fmo_data.bladder_indices;
  }

  /**
   * @brief Ritorna l'intera struttura FMOData
   */
  const FMOData &getData() const { return fmo_data; }

  /**
   * @brief Stampa un riassunto dello stato attuale
   */
  void printSummary() const {
    std::cout << "\n=== FMOData Summary ===" << std::endl;
    std::cout << "Total voxels: " << fmo_data.total_voxels << std::endl;
    std::cout << "Total beamlets: " << fmo_data.total_beamlets << std::endl;
    std::cout << "Influence matrix size: " << fmo_data.D.rows() << " x "
              << fmo_data.D.cols() << std::endl;
    std::cout << "Non-zero elements: " << fmo_data.D.nonZeros() << std::endl;
    if (fmo_data.D.rows() > 0 && fmo_data.D.cols() > 0) {
      std::cout << "Matrix density: "
                << (100.0 * fmo_data.D.nonZeros() /
                    (fmo_data.D.rows() * fmo_data.D.cols()))
                << "%" << std::endl;
    }
    std::cout << "PTV indices: " << fmo_data.ptv_indices.size() << std::endl;
    std::cout << "Rectum indices: " << fmo_data.rectum_indices.size()
              << std::endl;
    std::cout << "Bladder indices: " << fmo_data.bladder_indices.size()
              << std::endl;
  }
};

#endif // FMO_UTILS_HPP