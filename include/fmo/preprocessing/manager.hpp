/**
 * @file utils.hpp
 * @brief Classe FMODataManager per il caricamento e la manipolazione dei dati
 * FMO
 * @details La classe gestisce internamente un'istanza di FMOData e fornisce
 * metodi per caricare e manipolare le matrici sparse D (matrici di influenza) e
 * i vettori di indici VOILIST (Volume of Interest List) da file .mat. Gli
 * indici VOILIST sono utilizzati per identificare le regioni di interesse (ROI)
 * come PTV, retto e vescica.
 *
 */

#pragma once

#include "fmo/preprocessing/fmo.hpp"
#include <Eigen/Sparse>
#include <iostream>
#include <matio.h>
#include <string>
#include <string_view>
#include <vector>

namespace ROIFiles {
constexpr std::string_view default_directory = "../data/prostate";
constexpr std::string_view ptv_file =
    "PTV_68_VOILIST.mat"; // prendiamo di riferimento solo la regione tumorale
                          // senza considerare le sue aree di contorno
constexpr std::string_view rectum_file = "Rectum_VOILIST.mat";
constexpr std::string_view bladder_file = "Bladder_VOILIST.mat";
}; // namespace ROIFiles

/**
 * @class FMODataManagers
 * @brief Manager per il caricamento e la manipolazione dei dati FMO
 * @details Mantiene internamente un'istanza di FMOData e fornisce metodi per:
 *   - Caricamento di matrici D singole e multiple da file .mat
 *   - Concatenazione orizzontale di matrici D da multiple configurazioni di
 * fasci
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
  static Eigen::SparseMatrix<double>
  loadDMatrixHelper(const std::string &filepath, int &out_rows, int &out_cols) {
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
   * @brief Costruttore da FMOData esistente
   */
  explicit FMODataManager(const FMOData &data) : fmo_data(data) {}

  /**
   * @brief Imposta lo stato interno da una struttura FMOData
   */
  void setData(const FMOData &data) { fmo_data = data; }

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

    // std::cout << "✓ Matrice D caricata: " << rows << " x " << cols
    //           << " (nnz: " << D.nonZeros() << ")" << std::endl;
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
      // std::cout << "✓ PTV_VOILIST caricato: " << indices.size() << " voxel"
      //           << std::endl;
    } else if (roi_type == "rectum") {
      fmo_data.rectum_indices = indices;
      // std::cout << "✓ RECTUM_VOILIST caricato: " << indices.size() << "
      // voxel"
      //           << std::endl;
    } else if (roi_type == "bladder") {
      fmo_data.bladder_indices = indices;
      // std::cout << "✓ BLADDER_VOILIST caricato: " << indices.size() << "
      // voxel"
      //           << std::endl;
    } else {
      throw std::runtime_error("Tipo di ROI non riconosciuto: " + roi_type);
    }
  }

  /**
   * @brief Carica le ROI standard dalla directory predefinita.
   * @param directory Directory contenente i file VOILIST.
   */
  void loadDefaultROIs(
      const std::string &directory = std::string(ROIFiles::default_directory)) {
    loadVoiList(directory + "/" + std::string(ROIFiles::ptv_file), "v", "ptv");
    loadVoiList(directory + "/" + std::string(ROIFiles::rectum_file), "v",
                "rectum");
    loadVoiList(directory + "/" + std::string(ROIFiles::bladder_file), "v",
                "bladder");
  }

  /**
   * @brief Estrae le sottomatrici dedicate a ciascuna ROI (PTV, Retto, Vescica)
   *        e precalcola i vettori dei pesi per la dose media.
   * @details Esegue una singola scansione delle colonne della matrice sparsa D globale,
   *          ripartendo gli elementi non nulli nelle rispettive sottomatrici.
   */
  void extractROISubmatrices() {
    if (fmo_data.total_voxels <= 0 || fmo_data.total_beamlets <= 0) {
      throw std::runtime_error(
          "Impossibile estrarre sottomatrici: matrice D non inizializzata.");
    }

    const int total_voxels = fmo_data.total_voxels;
    const int cols = fmo_data.total_beamlets;

    // 1. Tabelle di lookup inverso (O(1)) per mappare gli indici globali dei voxel agli indici locali
    std::vector<int> map_ptv(total_voxels, -1);
    for (size_t i = 0; i < fmo_data.ptv_indices.size(); ++i) {
      map_ptv[fmo_data.ptv_indices[i]] = static_cast<int>(i);
    }

    std::vector<int> map_rectum(total_voxels, -1);
    for (size_t i = 0; i < fmo_data.rectum_indices.size(); ++i) {
      map_rectum[fmo_data.rectum_indices[i]] = static_cast<int>(i);
    }

    std::vector<int> map_bladder(total_voxels, -1);
    for (size_t i = 0; i < fmo_data.bladder_indices.size(); ++i) {
      map_bladder[fmo_data.bladder_indices[i]] = static_cast<int>(i);
    }

    // 2. Vettori di triplette per costruire le sottomatrici
    std::vector<Eigen::Triplet<double>> triplets_ptv;
    std::vector<Eigen::Triplet<double>> triplets_rectum;
    std::vector<Eigen::Triplet<double>> triplets_bladder;

    // Inizializza i pesi lineari medi
    fmo_data.rectum_weights.assign(cols, 0.0);
    fmo_data.bladder_weights.assign(cols, 0.0);

    const double rectum_norm =
        fmo_data.rectum_indices.empty()
            ? 1.0
            : static_cast<double>(fmo_data.rectum_indices.size());
    const double bladder_norm =
        fmo_data.bladder_indices.empty()
            ? 1.0
            : static_cast<double>(fmo_data.bladder_indices.size());

    // 3. Singola scansione per colonne della matrice CSC D
    for (int col = 0; col < cols; ++col) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(fmo_data.D, col); it;
           ++it) {
        int r = it.row();
        double val = it.value();

        if (int sub_r = map_ptv[r]; sub_r >= 0) {
          triplets_ptv.emplace_back(sub_r, col, val);
        }
        if (int sub_r = map_rectum[r]; sub_r >= 0) {
          triplets_rectum.emplace_back(sub_r, col, val);
          fmo_data.rectum_weights[col] += val / rectum_norm;
        }
        if (int sub_r = map_bladder[r]; sub_r >= 0) {
          triplets_bladder.emplace_back(sub_r, col, val);
          fmo_data.bladder_weights[col] += val / bladder_norm;
        }
      }
    }

    fmo_data.global_nonzeros = fmo_data.D.nonZeros();

    // 4. Costruzione e compressione delle sottomatrici sparse
    fmo_data.D_ptv.resize(fmo_data.ptv_indices.size(), cols);
    fmo_data.D_ptv.setFromTriplets(triplets_ptv.begin(), triplets_ptv.end());
    fmo_data.D_ptv.makeCompressed();

    fmo_data.D_rectum.resize(fmo_data.rectum_indices.size(), cols);
    fmo_data.D_rectum.setFromTriplets(triplets_rectum.begin(),
                                      triplets_rectum.end());
    fmo_data.D_rectum.makeCompressed();

    fmo_data.D_bladder.resize(fmo_data.bladder_indices.size(), cols);
    fmo_data.D_bladder.setFromTriplets(triplets_bladder.begin(),
                                       triplets_bladder.end());
    fmo_data.D_bladder.makeCompressed();

    fmo_data.is_preprocessed = true;

    // 5. Deallocazione della matrice globale D (~145 MB di RAM liberati)
    fmo_data.D = Eigen::SparseMatrix<double>();
  }

  /**
   * @brief Carica matrice di influenza e ROI necessarie al problema FMO.
   * @param data_directory Directory contenente le matrici D.
   * @param gantry_angles Angoli dei gantry da concatenare.
   * @param roi_directory Directory contenente i file VOILIST.
   */
  void loadAllData(const std::string &data_directory,
                   const std::vector<int> &gantry_angles,
                   const std::string &roi_directory =
                       std::string(ROIFiles::default_directory)) {
    loadGlobalDMatrixFromAngles(data_directory, gantry_angles);
    loadDefaultROIs(roi_directory);
    extractROISubmatrices();
  }

  /**
   * @brief Concatena orizzontalmente matrici D di più fasci nello stato interno
   * @details Crea la matrice globale D = [D^(1) | D^(2) | ... | D^(K)]
   * @param beam_filepaths Vettore di percorsi ai file .mat contenenti le
   * matrici D
   * @throw std::runtime_error Se la lista è vuota o le matrici hanno righe
   * incoerenti
   */
  void
  loadGlobalDMatrixFromFiles(const std::vector<std::string> &beam_filepaths) {
    if (beam_filepaths.empty()) {
      throw std::runtime_error("Lista di fasci vuota");
    }

    std::vector<Eigen::SparseMatrix<double>> D_matrices;
    int M = -1;
    int total_N = 0;

    // std::cout << "Caricamento matrici D per " << beam_filepaths.size()
    // << " fascio/i..." << std::endl;

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
      // std::cout << "  Fascio " << (k + 1) << ": " << rows << " x " << cols
      // << " (nnz: " << D.nonZeros() << ")" << std::endl;
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

    // std::cout << "Matrice globale concatenata: " << M << " x " << total_N
    // << " (nnz: " << D_global.nonZeros() << ")" << std::endl;
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
  const std::vector<int> &getPTVIndices() const { return fmo_data.ptv_indices; }

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
   * @brief Ritorna la sottomatrice PTV
   */
  const Eigen::SparseMatrix<double> &getDPTV() const { return fmo_data.D_ptv; }

  /**
   * @brief Ritorna la sottomatrice Retto
   */
  const Eigen::SparseMatrix<double> &getDRectum() const {
    return fmo_data.D_rectum;
  }

  /**
   * @brief Ritorna la sottomatrice Vescica
   */
  const Eigen::SparseMatrix<double> &getDBladder() const {
    return fmo_data.D_bladder;
  }

  /**
   * @brief Ritorna i pesi lineari medi precalcolati per il Retto
   */
  const std::vector<double> &getRectumWeights() const {
    return fmo_data.rectum_weights;
  }

  /**
   * @brief Ritorna i pesi lineari medi precalcolati per la Vescica
   */
  const std::vector<double> &getBladderWeights() const {
    return fmo_data.bladder_weights;
  }

  /**
   * @brief Verifica se le sottomatrici ROI sono state estratte
   */
  bool isPreprocessed() const { return fmo_data.is_preprocessed; }

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
    int nnz = fmo_data.is_preprocessed ? fmo_data.global_nonzeros
                                       : fmo_data.D.nonZeros();
    std::cout << "Influence matrix size: " << fmo_data.total_voxels << " x "
              << fmo_data.total_beamlets
              << (fmo_data.is_preprocessed
                      ? " (deallocated after ROI extraction)"
                      : "")
              << std::endl;
    std::cout << "Non-zero elements: " << nnz << std::endl;
    if (fmo_data.total_voxels > 0 && fmo_data.total_beamlets > 0 && nnz > 0) {
      std::cout << "Matrix density: "
                << (100.0 * nnz /
                    (static_cast<double>(fmo_data.total_voxels) *
                     fmo_data.total_beamlets))
                << "%" << std::endl;
    }
    std::cout << "PTV indices: " << fmo_data.ptv_indices.size() << std::endl;
    std::cout << "Rectum indices: " << fmo_data.rectum_indices.size()
              << std::endl;
    std::cout << "Bladder indices: " << fmo_data.bladder_indices.size()
              << std::endl;

    if (fmo_data.is_preprocessed) {
      std::cout << "--- Preprocessed ROI Submatrices ---" << std::endl;
      std::cout << "D_ptv:     " << fmo_data.D_ptv.rows() << " x "
                << fmo_data.D_ptv.cols()
                << " (nnz: " << fmo_data.D_ptv.nonZeros() << ")" << std::endl;
      std::cout << "D_rectum:  " << fmo_data.D_rectum.rows() << " x "
                << fmo_data.D_rectum.cols()
                << " (nnz: " << fmo_data.D_rectum.nonZeros() << ")" << std::endl;
      std::cout << "D_bladder: " << fmo_data.D_bladder.rows() << " x "
                << fmo_data.D_bladder.cols()
                << " (nnz: " << fmo_data.D_bladder.nonZeros() << ")"
                << std::endl;
    }
  }
};