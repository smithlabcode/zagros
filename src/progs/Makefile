#    Copyright (C) 2011 University of Southern California,
#                       Philip J. Uren, Andrew D. Smith
#
#    Authors: Philip J. Uren, Andrew D. Smith
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

ifndef SMITHLAB_CPP
$(error Must define SMITHLAB_CPP variable)
endif

PROGS = zagros

CXX = g++
CFLAGS = -g -Wall -fPIC -fmessage-length=50
OPTFLAGS = -O3
DEBUGFLAGS =
COMMON_DIR = ../rbp_common/
LIBS = -lgsl -lgslcblas

INCLUDEDIRS = $(SMITHLAB_CPP) $(COMMON_DIR)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
endif

ifdef DEBUG
CFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CFLAGS += $(OPTFLAGS)
endif

all:    $(PROGS)

%.o: %.cpp %.hpp
	$(CXX) $(CFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

zagros : $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o RNG.o FileIterator.o MappedRead.o)

zagros : $(addprefix $(COMMON_DIR), Model.o IO.o)

install: all
	@mkdir -p ../../bin
	@install -m 755 $(PROGS) ../../bin

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~
	@-rm -rf *.dSYM
	
test:
	make -C ../test SMITHLAB_CPP=$(SMITHLAB_CPP) test; \

.PHONY: clean
