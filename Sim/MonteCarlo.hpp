/*
 *
 * MCSim - Monte Carlo simulation of lattice polymers and proteins
 *
 * Copyright (C) 2006 - 2019 Thomas Wuest <monte.carlo.coder@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#ifndef MONTECARLO_H
#define MONTECARLO_H

#include "Model.hpp"


/*
  MonteCarlo class:

  This class provides the important control parameters used for
  Wang-Landau sampling (e.g. initial/final modification factors,
  flatness criterion etc.) and for Metropolis sampling
  (e.g. temperature) and contains the routines for these two Monte
  Carlo algorithms (i.e. Wang-Landau sampling and Metropolis
  sampling). In this code, the Monte Carlo algorithms and the specific
  physical models are strictly separated form each other. This makes
  it possible that the various algorithms can be used - principally -
  for any model that conforms to the "Model" class.
*/

class MonteCarlo {

public:

  // constructor
  // initializes the GSL random number generator and the Monte Carlo
  // algorithms either anew or from a Monte Carlo restart file
  // 1. input file name
  MonteCarlo(const char*);

  // destructor
  // "deletes" the GSL random number generator
  ~MonteCarlo();

  // the following routines perform a particular type of Monte Carlo
  // sampling and can be applied to any "physical model" class object
  // that is derived from the "Model" class:

  // Wang-Landau sampling (WLS)
  // 1. pointer to a "Histogram" object
  // 2. pointer to a "Model" object
  void WangLandauSampling(Histogram*, Model*);

  // Metropolis sampling (MS)
  // 1. pointer to a "Model" object
  void MetropolisSampling(Model*);

  // the following two parameters are "global" and initialized in the
  // "MonteCarlo" class constructor
  static int PhysicalModel;  // type of physical model
  static int MCAlgorithm;    // type of MC algorithm

private:

  // random number generator type
  int RngType;
  // random number generator seed
  unsigned long int RngSeed;

  // Wang-Landau sampling parameters:
  double ModFactorInit;          // initial modification factor (ln)
  double ModFactorDivider;       // modification factor divider
  double ModFactorThreshold;     // final modification factor (ln)
  double FlatnessMeasure;        // histogram flatness criterion
  unsigned long int MinHitsBin;  // min. hits per histogram bin
  // histogram check type, see "Histogram" class
  int HistogramCheckType;
  // # MC steps between successive flatness checks
  unsigned long int HistogramCheckInterval;
  // DOS update type, see "Histogram" class
  int DOSUpdateType;
  // # MC steps between successive storage of "Histogram" data
  unsigned long int DOSDumpInterval;
  // flag: 1 = mask all histogram entries as visited
  int MaskAllFlag;

  // Metropolis sampling parameters:
  unsigned long int EquilibrationMCSteps;  // # MC steps for equilibration
  unsigned long int SamplingMCSteps;       // # MC steps for sampling
  unsigned long int MeasuringInterval;     // # MC steps between successive measurements
  double SamplingTemperature;              // sampling temperature

  unsigned long int MCMoves;          // current # MC steps
  unsigned long int MCMovesAccepted;  // current # accepted MC steps
  unsigned long int MCMovesMem;
  unsigned long int MCMovesIteration;
  double ModFactor;                   // current modification factor
  int WLIteration;                    // current Wang-Landau iteration
  ObservableType Emin;                // current energy minimum

  // saves the Monte Carlo simulation "state" in order to restart the
  // simulation later
  void SaveMCState();

  // prints error messages on standard error stream
  // 1. error code
  // 2. string argument (optional)
  void ErrorMsg(int, const char* = NULL);

};


#endif
