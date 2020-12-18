# MCSim - Monte Carlo simulation of lattice polymers and proteins

Copyright (C) 2006 - 2019 Thomas Wüst

Contact & feedback: Thomas Wüst, <monte.carlo.coder@gmail.com>

Please send questions, comments, etc. to the above e-mail address.

---

## Table of Contents

This README file briefly covers the following topics:

1. Introduction

2. Credits / Acknowledgments

3. Structure of the program **MCSim**
    - List and short description of source files

4. Compilation

5. Running an example simulation (HP model 64mer in 2D)
    - Testing correct execution of the code
    - Short description of input parameters
    - Explanation of output files
    - Check-pointing / restarting mechanism

---

## Introduction

**MCSim** is a program for simulating coarse-grained lattice models of
polymers and proteins. Simplified lattice models are employed to study
fundamental questions in the statistical physics of polymers (e.g.,
phase transitions and critical phenomena) and bio-physics (e.g.,
protein folding). They are often the only means of investigation for
research questions in cases where experiments are not feasible or
detailed, atomistic models are computationally too demanding or
impractical to focus on the analysis of a particular phenomena.

At the core of **MCSim** are highly optimized data structures and
implementations for sophisticated Monte Carlo algorithms (e.g.,
Wang-Landau sampling) and efficient Monte Carlo trial move sets (e.g.,
"pull" and "bond-rebridging" moves) for polymers and proteins on the
simple (hyper-)cubic lattice in any dimension d >= 2.

In its present form, the code runs very efficient Wang-Landau
simulations of the hydrophobic-polar (HP) protein model and the
interacting self-avoiding walk (ISAW) polymer model in free space. In
particular, **MCSim** is able to compute the density of states (DOS) for
HP sequences with more than 100 monomers and ISAWs up to several
thousand monomers over their entire energy ranges with very high
numerical accuracy. The code has also proven to be very powerful for
the combinatorial optimization problem of finding new ground state
configurations for long HP sequences. **MCSim** is the first "blind
search" code that found the lowest ground state configuration for a
benchmark 103mer HP sequence.

Thanks to its efficient and reusable building blocks, the code can be
easily extended for Monte Carlo simulations of other physical systems,
such as protein interaction with or polymer grafted to surfaces,
protein membrane systems, multi-chain models, and proteins in crowded
environments.

---

## Credits / Acknowledgments

The first version of this code was developed by Thomas Wüst in 2006,
at the Center for Simulational Physics (CSP), University of
Georgia. Since then, Thomas Wüst has been continuously extending,
improving and optimizing the code and the current version of the
program is the result of these developments.

The code has been the basis for the following publications where its
underlying concepts are described in detail:

1. Thomas Wüst and David P. Landau, "Versatile approach to access
    the low temperature thermodynamics of lattice polymers and
    proteins," _Phys. Rev. Lett._ **102**, 178101 (2009),
    [ref](https://doi.org/10.1103/PhysRevLett.102.178101 "https://doi.org/10.1103/PhysRevLett.102.178101")
    [arXiv](https://arxiv.org/abs/1503.04433 "https://arxiv.org/abs/1503.04433").

2. Thomas Wüst and David P. Landau, "Optimized Wang-Landau sampling
    of lattice polymers: Ground state search and folding
    thermodynamics of HP model proteins," _J. Chem. Phys._ **137**, 064903
    (2012),
    [ref](https://doi.org/10.1063/1.4742969 "https://doi.org/10.1063/1.4742969")
    [arXiv](https://arxiv.org/abs/1207.3974 "https://arxiv.org/abs/1207.3974").

**_Please cite these two references whenever you use this code or parts
of it in your own research / work._**

Thanks to its efficiency, ease of use and extensibility, this code has
laid the foundation for the study of numerous research questions where
computer simulations of coarse-grained lattice models can serve as an
efficient means of understanding protein folding. Under the direction
of Prof. David P. Landau several CSP members have used the code or
added new features for the needs of their research. Thus, over the
years, many more publications resulted from these efforts. However,
the present code solely refers to the above two references.

---

## Structure of **MCSim**

Conceptual design: In this code the "physical model" and the Monte
Carlo algorithm are strictly separated from each other. This principle
enables the easy integration of different physical models and reuse of
existing Monte Carlo algorithms.

The program is written in C/C++.

The program consists of three major objects (C++ classes):

- `MonteCarlo`:

  > Monte Carlo algorithms, i.e., Wang-Landau sampling and Metropolis
    sampling.

- `Histogram`:

  > data structures and routines for the storing and handling of
    histograms.

- `Model`:

  > virtual base class that defines the necessary routines and data
    structures for the interaction of the physical model and the Monte
    Carlo algorithm.

Derived class from the virtual base class "Model":

- `HPModel`:

  > essential physical model class that defines the necessary data
    structures (e.g., monomer coordinates representations) and
    routines (e.g., Monte Carlo trial move sets, calculation of
    observables) for the hydrophobic-polar (HP) protein model and the
    interacting self-avoiding walk (ISAW) polymer model in free space
    on the simple (hyper-)cubic lattice in d dimensions.

### Short description of source files:

The header files (`*.H`) in the following list of source files usually
contain short descriptions to their corresponding classes and comments
on each parameter, data structure or routine.

- `Vector.C/.H`:

  > classes **Vector** and **Matrix** that comprise data structures
    and routines for vector / vector and vector / matrix operations,
    respectively.

- `Field.C/.H`:

  > virtual base class **Field** and its derived classes **BitTable**
    and **HashTable** that comprise data structures and routines
    (e.g., query, insert, delete operations) for the handling of bit
    and hash tables, respectively.

- `Globals.C/.H`:

  > generic parameters and routines, e.g., Boltzmann constant,
    factorial function.

- `Histogram.C/.H`:

  > class **Histogram** that comprises data structures and routines
    for the handling of multi-dimensional arrays such as histograms.

- `MonteCarlo.C/.H`:

  > class **MonteCarlo** that comprises the Monte Carlo algorithms,
    i.e., Wang-Landau sampling and Metropolis sampling.

- `Model.H`:

  > virtual base class **Model** that comprises the necessary routines
    and data structures for the interaction of the physical model and
    Monte Carlo algorithms.

- `HPModel.C/.H`:

  > derived physical model class **HPModel** (from the base class
    `Model`) that comprises data structures (e.g., monomer coordinates
    representations) and routines (e.g., Monte Carlo trial move sets,
    calculation of observables) for the hydrophobic-polar (HP) protein
    model and interacting self-avoiding walk (ISAW) polymer model in
    free space on the simple (hyper-)cubic lattice in d dimensions.

- `Main.C`:

  > `main()` routine which simply selects and instantiates a physical
    model class and runs the Monte Carlo simulation using either
    Wang-Landau or Metropolis sampling.

---

## Compilation

A `Makefile` is included with the source code. To compile the program,
run `make` in the directory of the `Makefile` and source files.

The program depends only on one external library, namely the GNU
Scientific Library [GSL](https://www.gnu.org/software/gsl/
"https://www.gnu.org/software/gsl/"). If GSL is installed correctly on
your system, the Makefile should be able to locate the relevant
directories of the GSL library automatically. Otherwise, you need to
specify the locations of the GSL header and library files manually.

Experimentation with different compilers has shown that the GNU
compiler [gcc/g++](https://gcc.gnu.org/ "https://gcc.gnu.org/") gives
very satisfactory performance.

`-O3 -mtune=generic` (GNU compiler) gives an optimized binary that
runs on both Intel and AMD CPUs (other options such as
`-march=amdfam10` do not change the performance significantly).

---

## Running an example simulation (HP model 64mer in 2D)

In the directory `./example` there is an input file (`seq.input`) to
run a Wang-Landau simulation to obtain the density of states (DOS) of
an 64mer HP protein in two dimensions (2D).

### Testing correct execution of the code:

Once you have successfully built the binary, change to the `./example`
directory and run

> `../MCSim seq.input > seq.output`

For an executable compiled with all optimizations turned on, this
simulation takes approximately 10 minutes on a present-day workstation
or notebook computer. Upon completion of the simulation, various files
are written to this directory (`./example`). The directory
`./example_complete` contains the correct files / results of the same
simulation. You should find exactly the same output (because the HP
model is a discrete model). You can check whether this is the case
with the Linux `diff` command

> `diff -r ./example ./example_complete`

If there are no differences (no output from the command), your code
runs correctly. This is also the case if you observe only rounding
errors in the very last digits of floating point numbers stored to
files (these discrepancies can arise due to slightly different
hardware architectures).

### Short description of input parameters:

Each class reads its needed input parameters / data from the same
input file (in the above example from the file `seq.input`). Lines in
the input file that start with a `#` or new line (`\n`) character are
ignored. Thus, comments need to start with a `#` character. The name
of an input parameter and its value(s) must be separated by at least
one space character. If not otherwise stated, the order of occurrences
of input parameters in the input file is irrelevant.

Some input parameters are initialized with reasonable default
values. These are overwritten in case the corresponding input
parameter receives an explicit value assigned in the input file.

**Important note:** Generally, the program checks the input for
erroneous or conflicting input parameters and aborts with a short
error message. Thus, it should not be possible to make the code crash
by specifying wrong parameters (feedback welcome if this should not be
the case). However, this does not mean that the program always
produces _correct_ results. It is easily possible to make the
simulation run endlessly if unreasonable values for certain input
parameters are chosen. Good knowledge of the underlying concepts, the
types of simulations performed and the influence of the input
parameters on the outcome is essential to make good use of this code.

In the following the most important input parameters are briefly
explained and for each parameter its type (e.g., integer value) and
default value(s) (in parantheses) are specified. There are several
other input parameters and many combinations thereof that are not
described here. The presented selection represents a reasonable set to
start with.

#### General input parameters:

- `RngType`

  > type of random number generator, see [GSL's random number generator algorithms](https://www.gnu.org/software/gsl/doc/html/rng.html#random-number-generator-algorithms "https://www.gnu.org/software/gsl/doc/html/rng.html#random-number-generator-algorithms")
  >
  > integer value in the interval [0-2] (0):  
  > 0 = ranlxd2 (very good randomness but relatively slow)  
  > 1 = ranlxd1 (good randomness, slightly faster than ranlxd2)  
  > 2 = mt19937 (reasonable randomness and very fast)

- `RngSeed`

  > random number seed
  >
  > integer value >= 0 (0)  
  > if equal to zero (= 0), the current date/time is used

- `MCAlgorithm`

  > type of Monte Carlo algorithm
  >
  > integer value in the interval [0-1] (0)  
  > 0 = Wang-Landau sampling  
  > 1 = Metropolis sampling

- `PhysicalModel`

  > type of physical model
  >
  > integer value (0), only one model currently defined  
  > 0 = HP / ISAW model in free space

#### Input parameters for the HP / ISAW model:

- `HPModelDimension`

  > spatial dimension
  >
  > integer value >= 2 (2)

- `OccupancyFieldType`

  > To perform self-avoidance and nearest neighbor checks efficiently,
    the coordinates of the monomers are mapped to an array data
    structure; this input parameter determines the type of this data
    structure.
  >
  > integer value in [0,1] (0)
  >
  > 0 = "bit table":

  >> fast but consumes a lot of memory; memory usage scales with N to
     the power of d (N = number of monomers, d = spatial dimension);
     this is the preferred choice for N <= ca. 500;

  > 1 = "hash table":

  >> almost as fast but consumes much less memory. See e.g.,
     T. H. Cormen et al., "Introduction to Algorithms," MIT Press
     (2009).

- `MoveFractions`

  > relative frequency of pull, bond-rebridging and pivot moves
  >
  > two floating point values in the interval [0,1] (1.0 1.0, i.e.,
    only pull moves)
  >
  > For instance, `MoveFractions 0.75 0.90` means 75% pull moves, 15%
    bond-rebridging moves and 10% pivot moves; at each Monte Carlo
    step, the type of trial move is chosen randomly; i.e., the
    specified values indicate probabilities.

- `HPSequence`

  > specification of the HP sequence (H and P) for the HP model
  >
  > sequence of chars ["H","P"] (no default values)
  >
  > this specification can go over several lines of the input file;
    each line must start with the keyword `HPSequence`; the length of
    the HP sequence per line must not exceed 100 characters. For this
    input parameter, the order of occurrences is relevant because the
    final HP sequence is the concatenation of the pieces from each
    occurrence.

- `NumberOfMonomers`

  > polymer length (i.e., number of "H" monomers) for the ISAW model
  >
  > integer value >= 3 (no default value)
  >
  > **Note:** It is mandatory to specify either the input parameter
    `HPSequence` or the input parameter `NumberOfMonomers`; however,
    both of these input parameters must not occur in the input file.

#### Input parameters for Wang-Landau sampling:

- `HistogramDims`

  > dimensions of the histogram data structures
  >
  > two integer values >= 1 (1 1)
  >
  > Here we consider only Wang-Landau sampling for one physical
    observable; thus, the arrays for the density of states (DOS) and
    histogram are one dimensional, i.e., `HistogramDims 1 1`.
  >
  > **Note:** `HistogramDims` must be specified before the input
    parameter `HistogramSpecs` (see below).

- `HistogramSpecs`

  > lower boundary (inclusive), upper boundary (exclusive), and bin
    width of the density of states and histogram arrays, respectively
  >
  > three integer values (no default values)
  >
  > **Note:** If the chosen Monte Carlo algorithm is Wang-Landau
    sampling, it is mandatory to specify values for this input
    parameter.

  > **Note:** For the HP model, the measured observable is "number of
    non-bonded HH contacts" and **_not_** energy; therefore, the lower
    boundary in the example input file (`seq.input`) is zero.

- `ModFactorInit`

  > natural logarithm of the initial modification factor
  >
  > floating point value >= 0 (1.0)

- `ModFactorDivider`

  > `log(f_new) = log(f_old) / ModFactorDivider`
  >
  > floating point value > 0 (2.0)

- `ModFactorThreshold`

  > natural logarithm of the final modification factor
  >
  > floating point value >= 0 (1.0e-8)

- `FlatnessMeasure`

  > flatness criterion
  >
  > floating point value >= 0 (0.8)

- `HistogramCheckInterval`

  > number of Monte Carlo steps between subsequent checks of the
    histogram flatness
  >
  > integer value >= 1 (1000000)

#### Input parameters for Metropolis sampling:

- `EquilibrationMCSteps`

  > number of Monte Carlo steps for equilibration
  >
  > integer value >= 0 (0)

- `SamplingMCSteps`

  > number of Monte Carlo steps for sampling
  >
  > integer value >= 0 (0)

- `MeasuringInterval`

  > number of Monte Carlo steps between successive measurements
  >
  > integer value >= 1 (1)

- `SamplingTemperature`

  > sampling temperature in reduced units
  >
  > floating point value > 0 (1.0)

---

## Explanation of output files

In addition to the standard output (in the above example redirected to
the file `seq.output`), additional files are created / updated during
the simulation. In the following these output files are briefly
explained for the example of a Wang-Landau simulation of an HP / ISAW
model:

- `seq.input`

  > general input file (can have any filename)

- `seq.output`

  > standard output (can be redirected to any file); first, all input
    parameters with their corresponding values, important
    characteristics and the initial state of the physical model, and
    the histogram specifications are returned. During the simulation,
    useful information is displayed; for instance, the current
    Wang-Landau iteration and DOS, newly found observable values, or
    the number of Monte Carlo steps needed to make a round trip
    through the histogram (up-/down itinerary); at the end of the
    simulation, the overall acceptance ratio and Monte Carlo trial
    move statistics is returned.
  >
  > **Note:** The exact output varies depending on the selection of
    the Monte Carlo algorithm and physical model.

- `confsEmin.xyz`

  > stores the coordinates of each newly found minimal energy HP model
    conformation

- `coords_current.xyz`
- `hdata_current.dat`
- `mcs_restart.dat`
- `rng_state`

  > These four files are required for the check-pointing / restarting
    mechanism (see below). They store the current configuration of the
    HP protein / ISAW polymer (`coords_current.xyz`), the current
    histogram (`hdata_current.dat`), some key values of the Monte
    Carlo algorithm (`mcs_restart.dat`) and the current state of the
    random number generator (`rng_state`), respectively.

- `hdata_iterationXXXXXXXX.dat`

  > histogram data after every completed Wang-Landau iteration
    (`XXXXXXXX` indicating the iteration)

- `dos.dat`

  > final, normalized density of states (DOS); the first column lists
    all the values of the measured observable sampled during the
    simulation within the boundaries specified in the input file (see
    input parameter `HistogramSpecs` above); the second column lists
    the corresponding values of the DOS
  >
  > **Note:** For the HP / ISAW model, the observable is "number of
    non-bonded HH contacts" and not energy (which is just the negative
    of the observable).

- `hdata_final.dat`

  > final histogram data (same as `hdata_current.dat` upon completion
    of the Wang-Landau simulation)

---

## Check-pointing / restarting mechanism

Monte Carlo simulations such as Wang-Landau or Metropolis sampling may
often last for hours or even days. To avoid losing the intermediate
results of a long-running Monte Carlo simulation due to an
(un-)predictable interruption of the program, a _check-pointing /
restarting mechanism_ has been implemented in the code. Thus, during a
Wang-Landau or Metropolis simulation, it is possible to stop the
program at any time (except at the very beginning) and restart it at a
later time. To restart the simulation, just call

> `../MCSim seq.input >> seq.output`

from the directory where the simulation has been started originally
(note the append `>>` here). No other user interaction is required.

A restarted simulation gives exactly the same results as if the
simulation were not interrupted; this is true in particular for the
calculated DOS which is the most important result from the
simulation. Otherwise, after a restart only small deviations in the
last digits of some floating point numbers in the `hdata_*` files may
arise and the final Monte Carlo trial move statistics differs because
it is (currently) only kept in memory and reset for each new run.

The check-pointing / restarting mechanism may fail if the program was
interrupted during an output operation to any of the files
`coords_current.xyz`, `hdata_current.dat`, `mcs_restart.dat` or
`rng_state`. However, because the computer time for these operations
is usually marginal compared to the overall run time, this situation
is very unlikely.

This check-pointing / restarting mechanism is very useful for running
simulations on a high performance computing environment where usually
a batch queuing system is in place with a limited wall clock time per
job.
