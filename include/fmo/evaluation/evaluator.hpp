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
#pragma once
#include "fmo/core/fitness.hpp"
#include "fmo/core/individual.hpp"
#include "fmo/preprocessing/fmo.hpp"
#include <algorithm>
#include <stdexcept>

#define PTV_DOSE 68 // dose prescritta al target PTV

class Evaluator {
private:
  const FMOData &fmo_data;

  // Calcolo legacy con la matrice globale D (usato come fallback se non pre-elaborato)
  std::vector<double> computeDosesLegacy(const std::vector<double> &beamlet_intensities) const {
    if (beamlet_intensities.size() != static_cast<size_t>(fmo_data.D.cols())) {
      throw std::invalid_argument("Numero di intensita' beamlet non coerente "
                                  "con le colonne della matrice D");
    }

    std::vector<double> doses(fmo_data.total_voxels, 0.0);
    for (int beamlet_index = 0; beamlet_index < fmo_data.D.cols();
         ++beamlet_index) {
      double intensity = beamlet_intensities[beamlet_index];
      if (intensity == 0.0) continue;
      for (Eigen::SparseMatrix<double>::InnerIterator it(fmo_data.D,
                                                         beamlet_index);
           it; ++it) {
        doses[it.row()] += it.value() * intensity;
      }
    }
    return doses;
  }

public:
  /**
   * @brief Costruttore della classe Evaluator
   * @param fmo_data Riferimento alla struttura FMOData contenente matrici e ROI
   */
  explicit Evaluator(const FMOData &fmo_data) : fmo_data(fmo_data) {}

  /**
   * @brief Calcola le dosi per i soli voxel appartenenti al target PTV.
   * @details Utilizza la sottomatrice D_ptv (solo righe PTV), allocando un buffer
   *          di dimensioni ridotte (es. ~54 KB invece di ~24 MB).
   */
  std::vector<double>
  computeDosesPTV(const std::vector<double> &beamlet_intensities) const {
    if (beamlet_intensities.size() != static_cast<size_t>(fmo_data.D_ptv.cols())) {
      throw std::invalid_argument("Numero di intensita' beamlet non coerente "
                                  "con le colonne della matrice D_ptv");
    }

    std::vector<double> doses(fmo_data.D_ptv.rows(), 0.0);
    for (int col = 0; col < fmo_data.D_ptv.cols(); ++col) {
      double intensity = beamlet_intensities[col];
      if (intensity == 0.0) continue;
      for (Eigen::SparseMatrix<double>::InnerIterator it(fmo_data.D_ptv, col);
           it; ++it) {
        doses[it.row()] += it.value() * intensity;
      }
    }
    return doses;
  }

  /**
   * @brief Calcola le dosi per i soli voxel del Retto (OAR).
   */
  std::vector<double>
  computeDosesRectum(const std::vector<double> &beamlet_intensities) const {
    if (beamlet_intensities.size() != static_cast<size_t>(fmo_data.D_rectum.cols())) {
      throw std::invalid_argument("Numero di intensita' beamlet non coerente "
                                  "con le colonne della matrice D_rectum");
    }

    std::vector<double> doses(fmo_data.D_rectum.rows(), 0.0);
    for (int col = 0; col < fmo_data.D_rectum.cols(); ++col) {
      double intensity = beamlet_intensities[col];
      if (intensity == 0.0) continue;
      for (Eigen::SparseMatrix<double>::InnerIterator it(fmo_data.D_rectum, col);
           it; ++it) {
        doses[it.row()] += it.value() * intensity;
      }
    }
    return doses;
  }

  /**
   * @brief Calcola le dosi per i soli voxel della Vescica (OAR).
   */
  std::vector<double>
  computeDosesBladder(const std::vector<double> &beamlet_intensities) const {
    if (beamlet_intensities.size() != static_cast<size_t>(fmo_data.D_bladder.cols())) {
      throw std::invalid_argument("Numero di intensita' beamlet non coerente "
                                  "con le colonne della matrice D_bladder");
    }

    std::vector<double> doses(fmo_data.D_bladder.rows(), 0.0);
    for (int col = 0; col < fmo_data.D_bladder.cols(); ++col) {
      double intensity = beamlet_intensities[col];
      if (intensity == 0.0) continue;
      for (Eigen::SparseMatrix<double>::InnerIterator it(fmo_data.D_bladder, col);
           it; ++it) {
        doses[it.row()] += it.value() * intensity;
      }
    }
    return doses;
  }

  /**
   * @brief Metodo legacy di calcolo globale delle dosi su tutti i voxel.
   */
  std::vector<double>
  computeDoses(const std::vector<double> &beamlet_intensities) const {
    return computeDosesLegacy(beamlet_intensities);
  }

  /**
   * @brief Calcola la fitness del target PTV (penalità quadratica underdose).
   * @param doses_ptv Vettore delle dosi dei soli voxel del PTV.
   */
  double computeFitnessPTV(const std::vector<double> &doses_ptv) const {
    if (doses_ptv.empty()) return 0.0;
    double total_dose_ptv = 0.0;
    for (double d : doses_ptv) {
      double underdose = PTV_DOSE - d;
      if (underdose > 0.0) {
        total_dose_ptv += underdose * underdose;
      }
    }
    return total_dose_ptv / doses_ptv.size();
  }

  // Alias per retrocompatibilità
  double computeFitnessPVT(const std::vector<double> &doses,
                           const Individual &individual) const {
    if (fmo_data.is_preprocessed && doses.size() == fmo_data.D_ptv.rows()) {
      return computeFitnessPTV(doses);
    }
    double total_dose_ptv = 0.0;
    for (int index : fmo_data.ptv_indices) {
      total_dose_ptv += std::pow(std::max(0.0, PTV_DOSE - doses[index]), 2);
    }
    return total_dose_ptv / fmo_data.ptv_indices.size();
  }

  /**
   * @brief Calcola la fitness del Retto (dose media).
   * @details Se disponibili, sfrutta i pesi precalcolati O(N). Altrimenti calcola
   *          tramite scansione della sottomatrice D_rectum.
   */
  double computeFitnessRectum(const std::vector<double> &beamlet_intensities) const {
    if (fmo_data.D_rectum.rows() == 0 && fmo_data.rectum_indices.empty()) {
      return 0.0;
    }

    // Modalità veloce con pesi precalcolati O(N)
    if (!fmo_data.rectum_weights.empty()) {
      double fitness = 0.0;
      for (size_t col = 0; col < beamlet_intensities.size(); ++col) {
        fitness += fmo_data.rectum_weights[col] * beamlet_intensities[col];
      }
      return fitness;
    }

    // Calcolo diretto dalla sottomatrice D_rectum
    if (fmo_data.D_rectum.rows() > 0) {
      double total_dose = 0.0;
      for (int col = 0; col < fmo_data.D_rectum.cols(); ++col) {
        double intensity = beamlet_intensities[col];
        if (intensity == 0.0) continue;
        for (Eigen::SparseMatrix<double>::InnerIterator it(fmo_data.D_rectum, col);
             it; ++it) {
          total_dose += it.value() * intensity;
        }
      }
      return total_dose / fmo_data.D_rectum.rows();
    }

    return 0.0;
  }

  // Overload legacy con signature (doses, individual)
  double computeFitnessRectum(const std::vector<double> &doses,
                              const Individual &individual) const {
    if (fmo_data.is_preprocessed && !fmo_data.rectum_weights.empty()) {
      return computeFitnessRectum(individual.genes);
    }
    double total_dose_rectum = 0.0;
    for (int index : fmo_data.rectum_indices) {
      total_dose_rectum += doses[index];
    }
    return total_dose_rectum / fmo_data.rectum_indices.size();
  }

  /**
   * @brief Calcola la fitness della Vescica (dose media).
   * @details Se disponibili, sfrutta i pesi precalcolati O(N). Altrimenti calcola
   *          tramite scansione della sottomatrice D_bladder.
   */
  double computeFitnessBladder(const std::vector<double> &beamlet_intensities) const {
    if (fmo_data.D_bladder.rows() == 0 && fmo_data.bladder_indices.empty()) {
      return 0.0;
    }

    // Modalità veloce con pesi precalcolati O(N)
    if (!fmo_data.bladder_weights.empty()) {
      double fitness = 0.0;
      for (size_t col = 0; col < beamlet_intensities.size(); ++col) {
        fitness += fmo_data.bladder_weights[col] * beamlet_intensities[col];
      }
      return fitness;
    }

    // Calcolo diretto dalla sottomatrice D_bladder
    if (fmo_data.D_bladder.rows() > 0) {
      double total_dose = 0.0;
      for (int col = 0; col < fmo_data.D_bladder.cols(); ++col) {
        double intensity = beamlet_intensities[col];
        if (intensity == 0.0) continue;
        for (Eigen::SparseMatrix<double>::InnerIterator it(fmo_data.D_bladder, col);
             it; ++it) {
          total_dose += it.value() * intensity;
        }
      }
      return total_dose / fmo_data.D_bladder.rows();
    }

    return 0.0;
  }

  // Overload legacy con signature (doses, individual)
  double computeFitnessBladder(const std::vector<double> &doses,
                               const Individual &individual) const {
    if (fmo_data.is_preprocessed && !fmo_data.bladder_weights.empty()) {
      return computeFitnessBladder(individual.genes);
    }
    double total_dose_bladder = 0.0;
    for (int index : fmo_data.bladder_indices) {
      total_dose_bladder += doses[index];
    }
    return total_dose_bladder / fmo_data.bladder_indices.size();
  }

  /**
   * @brief Valuta la fitness dell'individuo.
   * @details Se fmo_data e' pre-elaborato, calcola unicamente le dosi per il PTV
   *          (~54 KB) e calcola analiticamente la dose media di Retto e Vescica,
   *          evitando l'allocazione di 24 MB e la scansione di 3 milioni di voxel.
   */
  Fitness evaluate(const Individual &individual) const {
    if (fmo_data.is_preprocessed) {
      std::vector<double> doses_ptv = computeDosesPTV(individual.genes);
      return Fitness(computeFitnessPTV(doses_ptv),
                     computeFitnessRectum(individual.genes),
                     computeFitnessBladder(individual.genes));
    }

    // Modalità legacy su matrice globale D
    std::vector<double> doses = computeDosesLegacy(individual.genes);
    return Fitness(computeFitnessPVT(doses, individual),
                   computeFitnessRectum(doses, individual),
                   computeFitnessBladder(doses, individual));
  }
};