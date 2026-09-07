#pragma once

/**
 * @file metrics.hpp
 * @brief file che racchiude le metriche di prestazione usate per misurare la
 * stabilità dell'algoritmo tra le sue versioni (sequenziale e parallelo tramite
 * ff o OpenMP)
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <fmo/core/population.hpp>

/**
 * @brief Salva una popolazione su file in formato CSV.
 * @param pop Popolazione da salvare.
 * @param filepath Percorso del file di destinazione.
 */
inline void savePopulation(const Population &pop, const std::string &filepath) {
  std::ofstream out(filepath);
  if (!out.is_open()) {
    throw std::runtime_error("Impossibile aprire il file per la scrittura: " + filepath);
  }
  const size_t N = pop.size();
  const size_t G = (N > 0) ? pop.getIndividual(0).genes.size() : 0;
  out << "# population_size=" << N << " num_genes=" << G << "\n";
  out << "# rank,crowding_distance,ptv_fitness,rectal_fitness,bladder_fitness,genes...\n";
  out << std::setprecision(14);
  for (size_t i = 0; i < N; ++i) {
    const auto &ind = pop.getIndividual(i);
    out << ind.getRank() << ","
        << ind.getCrowdingDistance() << ","
        << ind.getFitness().getPTVFitness() << ","
        << ind.getFitness().getRectalFitness() << ","
        << ind.getFitness().getBladderFitness();
    for (double g : ind.genes) {
      out << "," << g;
    }
    out << "\n";
  }
}

/**
 * @brief Carica una popolazione da file CSV generato tramite savePopulation.
 * @param filepath Percorso del file da leggere.
 * @return Popolazione ricostruita.
 */
inline Population loadPopulation(const std::string &filepath) {
  std::ifstream in(filepath);
  if (!in.is_open()) {
    throw std::runtime_error("Impossibile aprire il file per la lettura: " + filepath);
  }
  Population pop;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::stringstream ss(line);
    std::string token;

    if (!std::getline(ss, token, ',')) continue;
    int rank = std::stoi(token);

    if (!std::getline(ss, token, ',')) continue;
    double cd = std::stod(token);

    if (!std::getline(ss, token, ',')) continue;
    double ptv = std::stod(token);

    if (!std::getline(ss, token, ',')) continue;
    double rectum = std::stod(token);

    if (!std::getline(ss, token, ',')) continue;
    double bladder = std::stod(token);

    std::vector<double> genes;
    while (std::getline(ss, token, ',')) {
      genes.push_back(std::stod(token));
    }

    Individual ind(genes, Fitness(ptv, rectum, bladder));
    ind.setRank(rank);
    ind.setCrowdingDistance(cd);
    pop.addIndividual(ind);
  }
  return pop;
}

/**
 * Struttura usata per il calcolo della Root Mean Square Error (RMSE) dei valori
 * di fitness delle popolazioni generate tramite NSGA2 sequenziale e parallelo
 */
struct PopulationDiscrepancy {
  double rmse_ptv;
  double rmse_rectum;
  double rmse_bladder;
  double rmse_total_fitness;
  double rmse_genes; // Scarto quadratico sui dosaggi dei singoli bixel
};

/**
 * @brief Calcola lo Spread Metric (Deb et al., 2002) sul primo fronte di Pareto
 * (F_1).
 * @details Misura l'uniformità di distribuzione spaziale delle soluzioni non
 * dominate. In assenza di una frontiera ottima teorica a priori, gli estremi
 * della frontiera vengono identificati come i minimi/massimi osservati
 * all'interno dello stesso fronte.
 * @param population Popolazione finale valutata e classificata (con rank
 * assegnati).
 * @return Valore di spread Delta in [0, 1] (0 indica distanza tra le soluzioni
 * perfettamente uniforme).
 */
inline double computeSpreadMetric(const Population &population) {
  // 1. Estrai gli indici degli individui appartenenti al
  // primo fronte di Pareto (rank == 1)
  std::vector<size_t> front_1_indices;
  for (size_t i = 0; i < population.size(); ++i) {
    if (population.getIndividual(i).rank == 1) {
      front_1_indices.push_back(i);
    }
  }

  size_t Q = front_1_indices.size();
  // Se il fronte contiene meno di 3 soluzioni, la dispersione interna non è
  // calcolabile
  if (Q <= 2) {
    return 0.0;
  }

  // Ordino gli indici per PTV fitness
  std::sort(front_1_indices.begin(), front_1_indices.end(),
            [&population](size_t a, size_t b) {
              return population.getIndividual(a).getFitness().getPTVFitness() <
                     population.getIndividual(b).getFitness().getPTVFitness();
            });

  // Calcola le distanze euclidee consecutive nello spazio degli obiettivi
  std::vector<double> d(Q - 1);
  double sum_d = 0.0;
  for (size_t i = 0; i < Q - 1; ++i) {
    const auto &fit_curr =
        population.getIndividual(front_1_indices[i]).getFitness();
    const auto &fit_next =
        population.getIndividual(front_1_indices[i + 1]).getFitness();

    double dp = fit_next.getPTVFitness() - fit_curr.getPTVFitness();
    double dr = fit_next.getRectalFitness() - fit_curr.getRectalFitness();
    double db = fit_next.getBladderFitness() - fit_curr.getBladderFitness();

    d[i] = std::sqrt(dp * dp + dr * dr + db * db);
    sum_d += d[i];
  }

  double d_bar = sum_d / static_cast<double>(Q - 1);
  if (d_bar <= 1e-9) {
    return 0.0; // Tutte le soluzioni collassano nello stesso punto
  }

  // 4. Somma delle deviazioni assolute rispetto alla media
  double sum_deviations = 0.0;
  for (double dist : d) {
    sum_deviations += std::abs(dist - d_bar);
  }

  // 5. Calcolo dello spread relativo normalizzato (Delta)
  return sum_deviations / (static_cast<double>(Q - 1) * d_bar);
}

/**
 * @brief Calcola del Root Mean Square Error (RMSE) tra la popolazione finale
 *        ottenuta in sequenziale e quella ottenuta in parallelo.
 * @param pop_seq Popolazione finale prodotta dalla baseline sequenziale.
 * @param pop_par Popolazione finale prodotta dal motore parallelo (OMP /
 * FastFlow).
 * @return Struttura contenente gli RMSE suddivisi per obiettivo e sui geni.
 */
inline PopulationDiscrepancy computePopulationRMSE(Population pop_seq,
                                                   Population pop_par) {
  if (pop_seq.size() != pop_par.size() || pop_seq.size() == 0) {
    throw std::invalid_argument(
        "Le popolazioni devono avere la medesima dimensione non nulla.");
  }

  const size_t N = pop_seq.size();

  // 1. Ordina entrambe le popolazioni per PTV fitness (e tie-break su altri obiettivi)
  // per allineare gli individui corrispondenti
  auto sortByObjectives = [](const Individual &a, const Individual &b) {
    if (a.getFitness().getPTVFitness() != b.getFitness().getPTVFitness()) {
      return a.getFitness().getPTVFitness() < b.getFitness().getPTVFitness();
    }
    if (a.getFitness().getRectalFitness() != b.getFitness().getRectalFitness()) {
      return a.getFitness().getRectalFitness() < b.getFitness().getRectalFitness();
    }
    return a.getFitness().getBladderFitness() < b.getFitness().getBladderFitness();
  };

  std::sort(pop_seq.individuals.begin(), pop_seq.individuals.end(), sortByObjectives);
  std::sort(pop_par.individuals.begin(), pop_par.individuals.end(), sortByObjectives);

  double sum_sq_ptv = 0.0;
  double sum_sq_rectum = 0.0;
  double sum_sq_bladder = 0.0;
  double sum_sq_genes = 0.0;
  size_t total_gene_count = 0;

  // 2. Accumula gli scarti quadratici tra individui allineati
  for (size_t i = 0; i < N; ++i) {
    const auto &f_seq = pop_seq.getIndividual(i).getFitness();
    const auto &f_par = pop_par.getIndividual(i).getFitness();

    double diff_ptv = f_seq.getPTVFitness() - f_par.getPTVFitness();
    double diff_rec = f_seq.getRectalFitness() - f_par.getRectalFitness();
    double diff_bla = f_seq.getBladderFitness() - f_par.getBladderFitness();

    sum_sq_ptv += diff_ptv * diff_ptv;
    sum_sq_rectum += diff_rec * diff_rec;
    sum_sq_bladder += diff_bla * diff_bla;

    // Scarto quadratico sui geni (intensità dei beamlet)
    const auto &genes_seq = pop_seq.getIndividual(i).genes;
    const auto &genes_par = pop_par.getIndividual(i).genes;
    for (size_t g = 0; g < genes_seq.size(); ++g) {
      double diff_g = genes_seq[g] - genes_par[g];
      sum_sq_genes += diff_g * diff_g;
    }
    total_gene_count += genes_seq.size();
  }

  PopulationDiscrepancy res;
  res.rmse_ptv = std::sqrt(sum_sq_ptv / N);
  res.rmse_rectum = std::sqrt(sum_sq_rectum / N);
  res.rmse_bladder = std::sqrt(sum_sq_bladder / N);
  res.rmse_total_fitness =
      std::sqrt((sum_sq_ptv + sum_sq_rectum + sum_sq_bladder) / (3.0 * N));
  res.rmse_genes =
      (total_gene_count > 0) ? std::sqrt(sum_sq_genes / total_gene_count) : 0.0;

  return res;
}