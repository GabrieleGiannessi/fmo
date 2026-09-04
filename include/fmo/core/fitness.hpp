/**
 * @file fitness.hpp
 * @brief Definizione della classe Fitness per rappresentare il punteggio di
 * fitness di un individuo. Si implementano più fitness function in relazione
 * agli obiettivi di ottimizzazione.
 * @details La classe Fitness rappresenta il punteggio di fitness di un
 * individuo nella popolazione genetica. Contiene metodi per accedere e
 * modificare il punteggio di fitness dei vari individui. La classe implementa
 * diverse funzioni di fitness in relazione agli obiettivi di ottimizzazione,
 * ovvero la minimizzazione della dose agli organi a rischio (OAR) e la
 * massimizzazione della dose al target (PTV). La classe Fitness è progettata
 * per essere utilizzata in combinazione con la classe Individual, che
 * rappresenta un individuo nella popolazione genetica.
 */
#pragma once
class Fitness {
public:
  double value_target_ptv;  // Punteggio di fitness relativo al target PTV
  double value_oar_rectal;  // Punteggio di fitness relativo al retto (OAR)
  double value_oar_bladder; // Punteggio di fitness relativo alla vescica (OAR)

  Fitness(double value_target_ptv, double value_oar_rectal,
          double value_oar_bladder)
      : value_target_ptv(value_target_ptv), value_oar_rectal(value_oar_rectal),
        value_oar_bladder(value_oar_bladder) {}

  Fitness() : value_target_ptv(0.0), value_oar_rectal(0.0),
              value_oar_bladder(0.0) {}
private:
  // getters
  //  metodo per ottenere il punteggio di fitness relativo al target PTV
  double getValueTargetPTV() const { return value_target_ptv; }

  // metodo per ottenere il punteggio di fitness relativo al retto (OAR)
  double getValueOARRectal() const { return value_oar_rectal; }

  // metodo per ottenere il punteggio di fitness relativo alla vescica (OAR)
  double getValueOARBladder() const { return value_oar_bladder; }

  // setters
  //  metodo di modifica del punteggio di fitness relativo al target PTV
  void setValueTargetPTV(double value) { this->value_target_ptv = value; }

  // metodo di modifica del punteggio di fitness relativo al retto (OAR)
  void setValueOARRectal(double value) { this->value_oar_rectal = value; }

  // metodo di modifica del punteggio di fitness relativo alla vescica (OAR)
  void setValueOARBladder(double value) { this->value_oar_bladder = value; }

  
};