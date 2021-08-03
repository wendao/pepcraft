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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <limits.h>

#include "Globals.hpp"
#include "MonteCarlo.hpp"


int MonteCarlo::PhysicalModel = 0;
int MonteCarlo::MCAlgorithm = 0;


MonteCarlo::MonteCarlo(const char* filename)
{

  char line[200];
  char pname[51];
  FILE* f;

  // initialize 'MonteCarlo' class parameters

  RngType = 0;
  RngSeed = 0;
  // Wang-Landau sampling
  ModFactorInit = 1.0;
  ModFactorDivider = 2.0;
  ModFactorThreshold = 1.0e-8;
  FlatnessMeasure = 0.8;
  MinHitsBin = 0;
  HistogramCheckType = 0;
  HistogramCheckInterval = 1000000;
  DOSUpdateType = 1;
  DOSDumpInterval = 0;
  MaskAllFlag = 0;
  // Metropolis sampling
  EquilibrationMCSteps = 0;
  SamplingMCSteps = 0;
  MeasuringInterval = 1;
  SamplingTemperature = 1.0;

  // read 'MonteCarlo' class parameters from input file
  f = fopen(filename, "r");
  if (f == NULL) ErrorMsg(1, filename);

  while (fgets(line, sizeof(line), f) != NULL) {

    if ((line[0] != '#') && (line[0] != '\n')) {

      if (sscanf(line, "%50s", pname) == 1) {

	if (!strcmp(pname, "PhysicalModel")) {
	  if (sscanf(line, "%*50s %d", &PhysicalModel) != 1) ErrorMsg(2, pname);
	  if ((PhysicalModel < 0) || (PhysicalModel > 1)) ErrorMsg(3, pname);
	  continue;
	}

	if (!strcmp(pname, "MCAlgorithm")) {
	  if (sscanf(line, "%*50s %d", &MCAlgorithm) != 1) ErrorMsg(2, pname);
	  if ((MCAlgorithm < 0) || (MCAlgorithm > 1)) ErrorMsg(3, pname);
	  continue;
	}

	if (!strcmp(pname, "RngType")) {
	  if (sscanf(line, "%*50s %d", &RngType) != 1) ErrorMsg(2, pname);
	  if ((RngType < 0) || (RngType > 2)) ErrorMsg(3, pname);
	  continue;
	}

	if (!strcmp(pname, "RngSeed")) {
	  if (sscanf(line, "%*50s %lu", &RngSeed) != 1) ErrorMsg(2, pname);
	  continue;
	}

	if (!strcmp(pname, "ModFactorInit")) {
	  if (sscanf(line, "%*50s %lf", &ModFactorInit) != 1) ErrorMsg(2, pname);
	  if (ModFactorInit < 0.0) ErrorMsg(3, pname);
	  continue;
	}

	if (!strcmp(pname, "ModFactorDivider")) {
	  if (sscanf(line, "%*50s %lf", &ModFactorDivider) != 1) ErrorMsg(2, pname);
	  if (ModFactorDivider <= 0.0) ErrorMsg(3, pname);
	  continue;
	}

	if (!strcmp(pname, "ModFactorThreshold")) {
	  if (sscanf(line, "%*50s %lf", &ModFactorThreshold) != 1) ErrorMsg(2, pname);
	  if (ModFactorThreshold < 0.0) ErrorMsg(3, pname);
	  continue;
	}

	if (!strcmp(pname, "FlatnessMeasure")) {
	  if (sscanf(line, "%*50s %lf", &FlatnessMeasure) != 1) ErrorMsg(2, pname);
	  if (FlatnessMeasure < 0.0) ErrorMsg(3, pname);
	  continue;
	}

	if (!strcmp(pname, "MinHitsBin")) {
	  if (sscanf(line, "%*50s %lu", &MinHitsBin) != 1) ErrorMsg(2, pname);
	  continue;
	}

	if (!strcmp(pname, "HistogramCheckType")) {
	  if (sscanf(line, "%*50s %d", &HistogramCheckType) != 1) ErrorMsg(2, pname);
	  if ((HistogramCheckType < 0) || (HistogramCheckType > 1)) ErrorMsg(3, pname);
	  continue;
	}

	if (!strcmp(pname, "DOSUpdateType")) {
	  if (sscanf(line, "%*50s %d", &DOSUpdateType) != 1) ErrorMsg(2, pname);
	  if ((DOSUpdateType < 0) || (DOSUpdateType > 1)) ErrorMsg(3, pname);
	  continue;
	}

	if (!strcmp(pname, "EquilibrationMCSteps")) {
	  if (sscanf(line, "%*50s %lu", &EquilibrationMCSteps) != 1) ErrorMsg(2, pname);
	  continue;
	}

	if (!strcmp(pname, "SamplingMCSteps")) {
	  if (sscanf(line, "%*50s %lu", &SamplingMCSteps) != 1) ErrorMsg(2, pname);
	  continue;
	}

	if (!strcmp(pname, "SamplingTemperature")) {
	  if (sscanf(line, "%*50s %lf", &SamplingTemperature) != 1) ErrorMsg(2, pname);
	  if (SamplingTemperature <= 0.0) ErrorMsg(3, pname);
	  continue;
	}

	if (!strcmp(pname, "MeasuringInterval")) {
	  if (sscanf(line, "%*50s %lu", &MeasuringInterval) != 1) ErrorMsg(2, pname);
	  if (MeasuringInterval < 1) ErrorMsg(3, pname);
	  continue;
	}

	if (!strcmp(pname, "HistogramCheckInterval")) {
	  if (sscanf(line, "%*50s %lu", &HistogramCheckInterval) != 1) ErrorMsg(2, pname);
	  if (HistogramCheckInterval < 1) ErrorMsg(3, pname);
	  continue;
	}

	if (!strcmp(pname, "DOSDumpInterval")) {
	  if (sscanf(line, "%*50s %lu", &DOSDumpInterval) != 1) ErrorMsg(2, pname);
	  continue;
	}

	if (!strcmp(pname, "MaskAllFlag")) {
	  if (sscanf(line, "%*50s %d", &MaskAllFlag) != 1) ErrorMsg(2, pname);
	  if ((MaskAllFlag < 0) || (MaskAllFlag > 1)) ErrorMsg(3, pname);
	  continue;
	}

      }

      else ErrorMsg(4);

    }

  }

  fclose(f);

  // now, construct class instance...

  // initialize random number generator
  switch (RngType) {
  case  1 : {rng = gsl_rng_alloc(gsl_rng_ranlxd1); break;}
  case  2 : {rng = gsl_rng_alloc(gsl_rng_mt19937); break;}
  default : {rng = gsl_rng_alloc(gsl_rng_ranlxd2);}
  }
  if (rng == NULL) ErrorMsg(5);

  f = fopen("mcs_restart.dat", "r");

  if (f == NULL) {

    RestartFlag = 0;

    MCMoves = 0;
    MCMovesAccepted = 0;
    MCMovesMem = 0;
    MCMovesIteration = 0;
    ModFactor = ModFactorInit;
    WLIteration = 1;
    Emin = 0;
    //    Emin = -1000000;

    if (RngSeed == 0) RngSeed = time(NULL) % ULONG_MAX;
    gsl_rng_set(rng, RngSeed);

  }

  else {

    RestartFlag = 1;

    fscanf(f, "%*[^\n]"); fscanf(f, "%*c");
    fscanf(f, "%*[^\n]"); fscanf(f, "%*c");
    fscanf(f, "%*[^\n]"); fscanf(f, "%*c");

    if (fscanf(f, "%*s %lu", &MCMoves) != 1) ErrorMsg(6);
    if (fscanf(f, "%*s %lu", &MCMovesAccepted) != 1) ErrorMsg(6);
    if (fscanf(f, "%*s %lu", &MCMovesMem) != 1) ErrorMsg(6);
    if (fscanf(f, "%*s %lu", &MCMovesIteration) != 1) ErrorMsg(6);
    if (fscanf(f, "%*s %lf", &ModFactor) != 1) ErrorMsg(6);
    if (fscanf(f, "%*s %d", &WLIteration) != 1) ErrorMsg(6);
    if (fscanf(f, "%*s %ld", &Emin) != 1) ErrorMsg(6);

    fclose(f);

    f = fopen("rng_state", "r");
    if (f == NULL) ErrorMsg(5);
    if (gsl_rng_fread(f, rng) != 0) ErrorMsg(5);
    fclose(f);

  }

  printf("\n### Monte Carlo Simulation Specifications ##############\n\n");

  if (RestartFlag)
    printf("Simulation: RESTART\n");
  else
    printf("Simulation: NEW\n");

  printf("\nRandom Number Generator:\n\n");
  printf(" Type  %s\n", gsl_rng_name(rng));
  if (RestartFlag) printf(" Seed  restart...\n");
  else printf(" Seed  %lu\n", RngSeed);
  printf(" Min/Max  %lu/%lu\n", gsl_rng_min(rng), gsl_rng_max(rng));

  switch (MCAlgorithm) {

  case  1 : {
    printf("\nMetropolis Sampling:\n\n");
    printf(" Number of Equilibration MC Steps  %lu\n", EquilibrationMCSteps);
    printf(" Number of Sampling MC Steps  %lu\n", SamplingMCSteps);
    printf(" Measuring Interval  %lu\n", MeasuringInterval);
    printf(" Sampling Temperature  %15.8e\n", SamplingTemperature);
    break;
  }

  default : {
    printf("\nWang-Landau (WL) Sampling:\n\n");
    printf(" Initial Modification Factor  %15.8e\n", ModFactorInit);
    printf(" Modification Factor Divider  %15.8e\n", ModFactorDivider);
    printf(" Modification Factor Threshold  %15.8e\n", ModFactorThreshold);
    printf(" Histogram Flatness Criterion  %15.8e\n", FlatnessMeasure);
    printf(" Min. Hits per Histogram Bin  %lu\n", MinHitsBin);
    printf(" Histogram Check Type  %d\n", HistogramCheckType);
    printf(" Histogram Check Interval  %lu\n", HistogramCheckInterval);
    printf(" DOS Update Type  %d\n", DOSUpdateType);
    printf(" DOS Dump Interval  %lu\n", DOSDumpInterval);
    printf(" MaskAllFlag  %d\n", MaskAllFlag);
  }

  }

  printf("\nAlgorithm Initialization:\n\n");
  printf(" MCMoves  %lu\n", MCMoves);
  printf(" MCMovesAccepted  %lu\n", MCMovesAccepted);
  printf(" MCMovesMem  %lu\n", MCMovesMem);
  printf(" MCMovesIteration  %lu\n", MCMovesIteration);
  printf(" ModFactor  %18.10e\n", ModFactor);
  printf(" WLIteration  %d\n", WLIteration);
  printf(" Emin  %ld\n", Emin);

  printf("\n########################################################\n");
}


MonteCarlo::~MonteCarlo()
{

  gsl_rng_free(rng);

}


void MonteCarlo::SaveMCState()
{

  FILE* f;

  f = fopen("mcs_restart.dat", "w");

  fprintf(f, "# Monte Carlo restart state\n");
  fprintf(f, "# Do not edit this file!\n");
  fprintf(f, "# Note: Random number generator state stored in 'rng_state'.\n");

  fprintf(f, "MCMoves  %lu\n", MCMoves);
  fprintf(f, "MCMovesAccepted  %lu\n", MCMovesAccepted);
  fprintf(f, "MCMovesMem  %lu\n", MCMovesMem);
  fprintf(f, "MCMovesIteration  %lu\n", MCMovesIteration);
  fprintf(f, "ModFactor  %18.10e\n", ModFactor);
  fprintf(f, "WLIteration  %d\n", WLIteration);
  fprintf(f, "Emin  %ld\n", Emin);

  fclose(f);

  f = fopen("rng_state", "w");
  if (gsl_rng_fwrite(f, rng) != 0) ErrorMsg(7);
  fclose(f);

}


void MonteCarlo::ErrorMsg(int message, const char* arg)
{

  fprintf(stderr, "\n");
  fprintf(stderr, "Error in 'MonteCarlo': ");

  switch (message) {

  case  1 : {
    fprintf(stderr, "Input file '%s' not found.\n", arg);
    break;
  }

  case  2 : {
    fprintf(stderr, "Parameter '%s' has invalid format.\n", arg);
    break;
  }

  case  3 : {
    fprintf(stderr, "Parameter '%s' is out of range.\n", arg);
    break;
  }

  case  4 : {
    fprintf(stderr, "Unreadable parameter.\n");
    break;
  }

  case  5 : {
    fprintf(stderr, "Random number generator initialization failed.\n");
    break;
  }

  case  6 : {
    fprintf(stderr, "Monte Carlo restart file unreadable.\n");
    break;
  }

  case  7 : {
    fprintf(stderr, "Saving random number generator state failed.\n");
    break;
  }

  case  8 : {
    fprintf(stderr, "Invalid model start for simulation.\n");
    break;
  }

  default : fprintf(stderr, "%s.\n", arg);

  }

  fprintf(stderr, "Program aborted.\n\n");
  exit(1);

}


void MonteCarlo::WangLandauSampling(Histogram* h, Model* m)
{

  unsigned long int n;
  long int prev, cur;

  if ((!RestartFlag) && (MaskAllFlag)) h -> MaskAll();

  m -> PrintObservables();
  prev = h -> GetIndex(m -> observable);
  if (prev < 0) ErrorMsg(8);

  while (ModFactor >= ModFactorThreshold) {

    for (n = 0; n < HistogramCheckInterval; n++) {

      MCMoves++;
      m -> DoMCMove();

      if (m -> MoveProposal) {

        cur = h -> GetIndex(m -> observable);

        if ((cur >= 0) && (gsl_rng_uniform(rng) < (exp(h -> GetDOS(prev) - h -> GetDOS(cur))))) {

          prev = cur;
          MCMovesAccepted++;

          // in the "HPModel" class the measured observable is "number
          // of non-bonded HH contacts"; the energy is minus this
          // quantity, therefore the minus sign
          // --> needs to be adapted for each physical model
          if (-m -> observable[0] < Emin) {
            Emin = -m -> observable[0];
            printf("New Emin = %ld ( MC moves = %lu )\n", Emin, MCMoves);
            fflush(stdout);
            //m -> WriteState(1, "confsEmin.mol2");
            //m -> WriteState(3, "confsEmin.xyz");
            m -> WriteState(2, "confsEmin.pdb");
            m->get_current_conf()->PrintInt();
          }
          else if (-m -> observable[0] == Emin) {
            // check if it is the same as the saved one
            m->get_current_conf()->PrintInt();
          }
          h -> CheckItinerancy(prev, MCMoves, MCMovesMem);
        }
        else m -> UnDoMCMove();

      }

      h -> Update(prev, n, DOSUpdateType, ModFactor);

      if ((DOSDumpInterval) && ((MCMoves % DOSDumpInterval) == 0)) {
        h -> SaveState(MCMoves / DOSDumpInterval, NULL, false);
      }

    }

    if (h -> CheckHistogram(HistogramCheckType, FlatnessMeasure, MinHitsBin)) {
      printf("DOS: iteration = %d ( mod. factor = %15.8e , MC moves : %lu %lu )\n",
	     WLIteration, ModFactor, MCMoves, MCMoves - MCMovesIteration);
      h -> SaveState(WLIteration, "hdata_iteration");
      h -> PrintNormDOS();
      WLIteration++;
      MCMovesIteration = MCMoves;
      ModFactor /= ModFactorDivider;
      h -> ResetHist();
    }
    //else printf("Histogram not flat ( MC moves = %lu )\n", MCMoves);

    h -> SaveState();
    //m -> WriteState(4, "coords_current.xyz");
    SaveMCState();
    fflush(stdout);

  }

  printf("Wang-Landau (WL) sampling done.\n");
  printf("MC moves = %lu\n", MCMoves);
  printf("Accepted MC moves = %lu\n", MCMovesAccepted);
  printf("Acceptance ratio = %15.8e\n", double(MCMovesAccepted) / double(MCMoves));

  h -> SaveState(0, "hdata_final.dat");
  h -> PrintNormDOS("dos.dat");

}


void MonteCarlo::MetropolisSampling(Model* m)
{

  unsigned long int n;
  ObservableType prev, cur;

  m -> PrintObservables();

  prev = -(m -> observable[0]);

  if (MCMoves == 0) {

    printf("Equilibrating...\n");
    fflush(stdout);

    for (n = 1; n <= EquilibrationMCSteps; n++) {

      m -> DoMCMove();

      if (m -> MoveProposal) {

        cur = -(m -> observable[0]);

        if (gsl_rng_uniform(rng) < (exp(double(prev - cur) / SamplingTemperature))) {

          prev = cur;

        }

        else m -> UnDoMCMove();

      }

      if ((n % MeasuringInterval) == 0) {
        printf("%lu %ld\n", n, m -> observable[0]);
        fflush(stdout);
      }

    }

    printf("Equilibration done.\n");
    fflush(stdout);
    //m -> WriteState(3, "conf_equi.xyz");
    m -> WriteState(2, "conf_equi.pdb");
  }

  while (MCMoves < SamplingMCSteps) {

    MCMoves++;
    m -> DoMCMove();

    if (m -> MoveProposal) {

      cur = -(m -> observable[0]);

      if (gsl_rng_uniform(rng) < (exp(double(prev - cur) / SamplingTemperature))) {

	prev = cur;
	MCMovesAccepted++;

	if (-m -> observable[0] < Emin) {
	  Emin = -m -> observable[0];
	  printf("New Emin = %ld ( MC moves = %lu )\n", Emin, MCMoves);
	  fflush(stdout);
	  //m -> WriteState(1, "confsEmin.mol2");
	  //m -> WriteState(3, "confsEmin.xyz");
    m -> WriteState(2, "confsEmin.pdb");
	}

      }

      else m -> UnDoMCMove();

    }

    if ((MCMoves % MeasuringInterval) == 0) {

      printf("\nMetropolis sampling: %lu MC moves\n", MCMoves);
      m -> PrintObservables();
      printf("\n");
      //m -> WriteState(4, "coords_current.xyz");
      m -> get_current_conf() -> PrintInt();
      SaveMCState();
      fflush(stdout);

    }

  }

  printf("Metropolis sampling done.\n");
  printf("MC moves = %lu\n", MCMoves);
  printf("Accepted MC moves = %lu\n", MCMovesAccepted);
  printf("Acceptance ratio = %15.8e\n", double(MCMovesAccepted) / double(MCMoves));

}
