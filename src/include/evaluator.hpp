/**
 * @file evaluator.hpp
 * @brief Classe per la valutazione della bontà degli individui di una
 * popolazione del modello FMO.
 * @details La classe Evaluator fornisce metodi per calcolare il punteggio di
 * fitness di un individuo in base agli obiettivi di ottimizzazione del problema
 * FMO. Gli obiettivi includono la minimizzazione della dose agli organi a
 * rischio (OAR) e la massimizzazione della dose al target (PTV). La classe
 * utilizza la matrice di influenza D e i vettori di indici VOILIST per
 * calcolare le dosi agli organi a rischio e al target, e restituisce un
 * punteggio di fitness complessivo per l'individuo.
 */

#include <fmo.hpp>
#include "individual.hpp"
class Evaluator {
public:
  /**
   * @brief Costruttore della classe Evaluator
   * @param fmo_data Riferimento alla struttura FMOData contenente la matrice di influenza D e i vettori di indici VOILIST
   */
  FMOData fmo_data;
  Evaluator(const FMOData &fmo_data) : fmo_data(fmo_data) {}
private: 

//metodo di calcolo delle dosi per ogni voxel utilizzando la matrice di influenza D
std::vector<double> computeDoses(const std::vector<double> &beamlet_intensities) {
    std::vector<double> doses(fmo_data.total_voxels, 0.0);
    for (int i = 0; i < fmo_data.total_voxels; ++i) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(fmo_data.D, i); it; ++it) {
        int beamlet_index = it.col();
        doses[i] += it.value() * beamlet_intensities[beamlet_index];
      }
    }
    return doses;
  }

double computeFitnessPVT(Individual<double> &individual) {
    std::vector<double> doses = computeDoses(individual.genes);
    double total_dose_ptv = 0.0;
    for (int index : fmo_data.ptv_indices) {
      total_dose_ptv += doses[index];
    }
    return total_dose_ptv / fmo_data.ptv_indices.size();
  }

  double computeFitnessRectum(Individual<double> &individual, const std::vector<int> &rectum_indices) {
    std::vector<double> doses = computeDoses(individual.genes);
    double total_dose_oar = 0.0;
    for (int index : rectum_indices) {
      total_dose_oar += doses[index];
    }
    return total_dose_oar / rectum_indices.size();
  }

    double computeFitnessBladder(Individual<double> &individual, const std::vector<int> &bladder_indices) {
        std::vector<double> doses = computeDoses(individual.genes);
        double total_dose_oar = 0.0;
        for (int index : bladder_indices) {
        total_dose_oar += doses[index];
        }
        return total_dose_oar / bladder_indices.size();
    }


};