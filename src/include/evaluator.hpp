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
#include "fitness.hpp"
#include "individual.hpp"
#include "src/preprocessing/fmo.hpp"
#include <algorithm>

#define PTV_DOSE 68 //dose prescritta al target PTV

class Evaluator {
public:
  /**
   * @brief Costruttore della classe Evaluator
   * @param fmo_data Riferimento alla struttura FMOData contenente la matrice di
   * influenza D e i vettori di indici VOILIST
   */
  FMOData fmo_data;
  Evaluator(const FMOData &fmo_data) : fmo_data(fmo_data) {}

  // metodo di calcolo delle dosi per ogni voxel utilizzando la matrice di
  // influenza D
  std::vector<double>
  computeDoses(const std::vector<double> &beamlet_intensities) {
    std::vector<double> doses(fmo_data.total_voxels, 0.0);
    for (int i = 0; i < fmo_data.total_voxels; ++i) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(fmo_data.D, i); it;
           ++it) {
        int beamlet_index = it.col();
        doses[i] += it.value() * beamlet_intensities[beamlet_index];
      }
    }
    return doses;
  }

  //metodi di calcolo delle singole fitness function per PTV, retto e vescica
  double computeFitnessPVT(const std::vector<double> &doses, const Individual &individual) {
    double total_dose_ptv = 0.0;
    for (int index : fmo_data.ptv_indices) {
      total_dose_ptv += std::pow(std::max(0.0, PTV_DOSE - doses[index]), 2);
    }
    return total_dose_ptv / fmo_data.ptv_indices.size();
  }

  double computeFitnessRectum(const std::vector<double> &doses, const Individual &individual) {
    double total_dose_rectum = 0.0;
    for (int index : fmo_data.rectum_indices) {
      total_dose_rectum += doses[index];
    }
    return total_dose_rectum / fmo_data.rectum_indices.size();
  }

  double computeFitnessBladder(const std::vector<double> &doses, const Individual &individual) {
    double total_dose_bladder = 0.0;
    for (int index : fmo_data.bladder_indices) {
      total_dose_bladder += doses[index];
    }
    return total_dose_bladder / fmo_data.bladder_indices.size();
  }

  //metodo di calcolo delle fitness function, ritorna un Fitness object contenente i punteggi di fitness per PTV, retto e vescica
  Fitness evaluate(const Individual &individual) {
    std::vector<double> doses = computeDoses(individual.genes);
    return Fitness(computeFitnessPVT(doses, individual),
                   computeFitnessRectum(doses, individual),
                   computeFitnessBladder(doses, individual));
  }
};