#
#
# MCSim - Monte Carlo simulation of lattice polymers and proteins
#
# Copyright (C) 2006 - 2019 Thomas Wuest <monte.carlo.coder@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#

.SUFFIXES:
.SUFFIXES: .cpp .o


#CFLAGS = -U_FORTIFY_SOURCE -g -pg -Wall -Wextra -pedantic
CFLAGS = -U_FORTIFY_SOURCE -O3 -mtune=generic


# GNU Scientific Library (GSL) support:
GSL_C := $(shell gsl-config --cflags)
GSL_L := $(shell gsl-config --libs)


# Object files:
OBJS = Vector.o      \
       Field.o       \
       Globals.o     \
       Histogram.o   \
       MonteCarlo.o  \
       HPModel.o     \
       HIPModel.o     \
       Main.o


# Target file:
TARGET = MCSim

.cpp.o:
	$(CXX) $(CFLAGS) $(GSL_C) -c $*.cpp

$(TARGET): $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $(GSL_L)


Vector.o: Vector.cpp Vector.hpp
Field.o: Field.cpp Field.hpp
Globals.o: Globals.cpp Globals.hpp
Histogram.o: Histogram.cpp Histogram.hpp
MonteCarlo.o: MonteCarlo.cpp Globals.hpp Histogram.hpp Model.hpp MonteCarlo.hpp
HPModel.o: HPModel.cpp Vector.hpp Field.hpp Globals.hpp Histogram.hpp Model.hpp HPModel.hpp
HIPModel.o: HIPModel.cpp Vector.hpp Field.hpp Globals.hpp Histogram.hpp Model.hpp HIPModel.hpp
Main.o: Main.cpp Globals.hpp MonteCarlo.hpp Histogram.hpp Model.hpp HPModel.hpp HIPModel.hpp


cleanup:
	rm -f *.o $(TARGET)
