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

PROGS = zagros thermo extractDEs scan

CXX = g++
CFLAGS = -g -Wall -fPIC -fmessage-length=50
OPTFLAGS = -O3
DEBUGFLAGS =
COMMON_DIR = ../rbp_common/
LIBS =

INCLUDEDIRS = $(SMITHLAB_CPP) $(COMMON_DIR)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
ifeq "$(shell if [ `sysctl -n kern.osrelease | cut -d . -f 1` -ge 13 ];\
              then echo 'true'; fi)" "true"
CFLAGS += -stdlib=libstdc++
endif
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

scan : $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
  smithlab_utils.o OptionParser.o RNG.o MappedRead.o)
scan : $(addprefix $(COMMON_DIR), Model.o IO.o RNA_Utils.o Part_Func.o)

zagros : $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o RNG.o MappedRead.o)

zagros : $(addprefix $(COMMON_DIR), Model.o IO.o RNA_Utils.o Part_Func.o)

thermo : $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o RNG.o MappedRead.o)

thermo : $(addprefix $(COMMON_DIR), RNA_Utils.o Part_Func.o IO.o)

simulate : $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o RNG.o MappedRead.o)

simulate : $(addprefix $(COMMON_DIR), RNA_Utils.o Part_Func.o IO.o)

extractDEs : $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o MappedRead.o)

extractDEs : $(addprefix $(COMMON_DIR), IO.o)

install: all
	@mkdir -p ../../bin
	@install -m 755 $(PROGS) ../../bin

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~
	@-rm -rf *.dSYM
	
test:
	make -C ../test SMITHLAB_CPP=$(SMITHLAB_CPP) test; \

.PHONY: clean
