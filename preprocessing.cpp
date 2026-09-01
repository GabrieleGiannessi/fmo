#include <iostream>
#include <vector>
#include <string>
#include <matio.h>
#include <Eigen/Sparse>

// Struttura per memorizzare il problema FMO pre-elaborato
struct FMOData {
    int total_voxels = 0;       // M (es. 690373 per Prostate)
    int total_beamlets = 0;     // N (somma dei bixel dei fasci scelti)
    Eigen::SparseMatrix<double> D; // Matrice globale concatenata [M x N]
    std::vector<int> ptv_indices;
    std::vector<int> rectum_indices;
    std::vector<int> bladder_indices;
};

// Funzione per caricare una singola matrice sparsa D da un file .mat
Eigen::SparseMatrix<double> loadDMatrix(const std::string& filepath, int& out_rows, int& out_cols) {
    mat_t *matfp = Mat_Open(filepath.c_str(), MAT_ACC_RDONLY);
    if (!matfp) {
        throw std::runtime_error("Errore nell'apertura del file: " + filepath);
    }

    matvar_t *matvar = Mat_VarRead(matfp, "D");
    if (!matvar || matvar->class_type != MAT_C_SPARSE) {
        Mat_Close(matfp);
        throw std::runtime_error("Variabile D sparsa non trovata in: " + filepath);
    }

    mat_sparse_t *sparseData = static_cast<mat_sparse_t*>(matvar->data);
    out_rows = matvar->dims[0];
    out_cols = matvar->dims[1];

    // Creazione della matrice sparsa Eigen a partire dal formato CSC di matio
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(sparseData->ndata);

    double *data = static_cast<double*>(sparseData->data);
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
std::vector<int> loadVoiList(const std::string& filepath, const std::string& var_name) {
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
        double *data = static_cast<double*>(matvar->data);
        for (int i = 0; i < num_elements; ++i) {
            // Conversione da 1-based (MATLAB) a 0-based (C++)
            indices[i] = static_cast<int>(data[i]) - 1;
        }
    } else if (matvar->data_type == MAT_T_INT32) {
        int32_t *data = static_cast<int32_t*>(matvar->data);
        for (int i = 0; i < num_elements; ++i) {
            indices[i] = data[i] - 1;
        }
    }

    Mat_VarFree(matvar);
    Mat_Close(matfp);
    return indices;
}