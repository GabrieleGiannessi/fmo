#include <ff/ff.hpp>
#include <ff/node.hpp>
#include <ff/parallel_for.hpp>

#include "fmo/evaluation/evaluator.hpp"
#include "fmo/nsga2/nsga-utils.hpp"
#include "fmo/nsga2/offspring.hpp"
#include "fmo/utilities/hpc_helpers.hpp"

using namespace ff;

/**
 * @brief valutazione della fitness usando il costrutto ParallelFor di FastFlow
 * @param population popolazione di individui
 * @param evaluator oggetto evaluator per il calcolo delle fitness
 * @param nw number of workers
 */
void evaluatePopulationFFParFor(Population &population, Evaluator &evaluator,
                                ParallelFor &pf) {
  // partizionamento statico di default
  // l'intervallo [0, size) viene diviso in blocchi contigui equi tra i worker.
  pf.parallel_for(0, population.size(), [&population, &evaluator](const size_t i) {
    auto &ind = population.getIndividual(i);
    ind.setFitness(evaluator.evaluate(ind));
  });
}

struct EvalTask {
  Individual *ptr_i; // puntatore all'individuo della popolazione
};

class SolutionEmitter : public ff_node_t<EvalTask> {
private:
  Population &population;

public:
  SolutionEmitter(Population &pop) : population(pop) {}

  EvalTask *svc(EvalTask *) {
    for (size_t i = 0; i < population.size(); ++i) {
      // Alloca un puntatore al task e lo spara nel canale lock-free
      ff_send_out(new EvalTask{&population.getIndividual(i)});
    }
    // Segnala la fine dello stream ai worker
    return EOS;
  }
};

class EvalWorker : public ff_node_t<EvalTask> {
private:
  Evaluator &evaluator;

public:
  EvalWorker(Evaluator &e) : evaluator(e) {}

  EvalTask *svc(EvalTask *ind) {
    if (ind == nullptr) {
      return EOS;
    }
    ind->ptr_i->setFitness(evaluator.evaluate(*(ind->ptr_i)));
    return ind;
  }
};

class SolutionCollector : public ff_node_t<EvalTask> {
public:
  EvalTask *svc(EvalTask *task) override {
    // Libera il descrittore del task
    delete task;
    return GO_ON; // Rimane in ascolto del prossimo task
  }
};

void evaluatePopulationFFFarm(Population &population, Evaluator &evaluator,
                              int nw) {

  SolutionEmitter emitter(population);
  SolutionCollector collector;

  std::vector<ff_node*> workers;
  for (int i = 0; i < nw; i++) {
    workers.push_back(new EvalWorker(evaluator));
  }

  ff_farm f;
  f.add_emitter(&emitter);
  f.add_workers(workers);
  f.add_collector(&collector);

  f.run_and_wait_end();
}

/**
 * @brief generazione della nuova prole distribuendo il calcolo dei nuovi
 * individui su più cores (nw). In particolare, ogni thread riceve il seme base
 * e calcola il suo PRNG per effettuare poi le chiamate ai metodi di selezione a
 * torneo, crossover e mutazione polinomiale.
 */
Population generatePopulationOffspringFF(const Population &population,
                                         double eta_c, double eta_m,
                                         uint64_t base_seed, int nw);
