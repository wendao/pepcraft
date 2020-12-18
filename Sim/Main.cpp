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

#include "Globals.hpp"
#include "MonteCarlo.hpp"
#include "Histogram.hpp"
#include "Model.hpp"
#include "HPModel.hpp"
#include "HIPModel.hpp"


// this is the main program
int main(int argc, char* argv[])
{

  MonteCarlo* mcs;
  Histogram* histogram = NULL;
  Model* model;

  // check whether there is an input file argument in the program call
  if (argc != 2) {
    printf("\n");
    printf("Usage: MCSim <input file>\n");
    printf("\n");
    exit(0);
  };

  // instantiate a "MonteCarlo" class object
  printf("db00");
  mcs = new MonteCarlo(argv[1]);

  // select and instantiate a "physical model" class object
  // (currently, only "HPModel")
 
  switch (MonteCarlo::PhysicalModel) {

    // in principal, any "physical model" class derived from the basis
    // class "Model" can be instantiated here

    case  1 : {
      model = new HIPModel(argv[1]);
      break;
    }

  default : {
    model = new HPModel(argv[1]);
  }

  }

  // run Monte Carlo simulation, either Wang-Landau sampling (default)
  // or Metropolis sampling
  switch (MonteCarlo::MCAlgorithm) {

  case  1 : {
    mcs -> MetropolisSampling(model);
    break;
  }

  default : {
    // instantiate a "Histogram" class object for Wang-Landau sampling
    histogram = new Histogram(argv[1]);
    mcs -> WangLandauSampling(histogram, model);
  }

  }

  delete model;
  if (histogram != NULL) delete histogram;
  delete mcs;

  printf("\nSimulation successfully completed.\n");

}
