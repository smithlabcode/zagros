#    Copyright (C) 2008 Cold Spring Harbor Laboratory, 
#                       University of Southern California and
#                       Andrew D. Smith
#
#    Authors: Emad Bahrami Samani, Philip J. Uren, Andrew D. Smith
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

package = zagros
version = 1.0.0
tarname = $(package)
distdir = $(tarname)-$(version)

RBP = $(shell pwd)
export PATH := $(shell pwd):$(PATH)

BINDIR = $(RBP)/bin
export SMITHLAB_CPP := $(shell pwd)/src/smithlab_cpp

all:
	@make -C src RBP=$(RBP) OPT=1 install

install:
	@export PATH=$(PATH)
	@make -C src RBP=$(RBP) OPT=1 install

developmentDocs:
	@doxygen $(RBP)/src/doxygen.config
.PHONY: developmentDocs

test:
	@make -C src OPT=1 test
.PHONY: test 

clean:
	@rm -rf $(RBP)/bin
	@make -C src RBP=$(RBP) clean
.PHONY: clean

distclean: clean
	@make -C src OPT=1 clean
	@rm -rf $(RBP)/bin
	@rm -rf $(distdir) $(distdir).tar.gz
.PHONY: distclean

dist: $(distdir).tar.gz

distcheck : $(distdir).tar.gz
	gzip -cd $(distdir).tar.gz | tar xvf -
	cd $(distdir) && ./configure && $(MAKE) all && $(MAKE) test
	cd $(distdir) && $(MAKE) clean
	rm -rf $(distdir)
	@echo "*** Package $(distdir).tar.gz is ready for distribution"

$(distdir).tar.gz : $(distdir)
	tar chof - $(distdir) | gzip -9 -c > $@
	rm -rf $(distdir)

$(distdir) : FORCE
	# make the directory structure 
	mkdir -p $(distdir)/src
	mkdir -p $(distdir)/src/rbp_common
	mkdir -p $(distdir)/src/progs
	mkdir -p $(distdir)/src/smithlab_cpp
	# copy top level files
	cp Makefile $(distdir)
	cp README.TXT $(distdir)
	# copy smithlab_cpp src files
	cp src/Makefile $(distdir)/src
	cp src/smithlab_cpp/Makefile $(distdir)/src/smithlab_cpp
	cp src/smithlab_cpp/*.cpp $(distdir)/src/smithlab_cpp
	cp src/smithlab_cpp/*.hpp $(distdir)/src/smithlab_cpp
	# copy the programs
	cp src/progs/extractDEs.cpp $(distdir)/src/progs
	cp src/progs/thermo.cpp $(distdir)/src/progs
	cp src/progs/zagros.cpp $(distdir)/src/progs
	cp src/progs/Makefile $(distdir)/src/progs
	# copy the common files
	cp src/rbp_common/*.cpp $(distdir)/src/rbp_common
	cp src/rbp_common/*.hpp $(distdir)/src/rbp_common
	cp src/rbp_common/Makefile $(distdir)/src/rbp_common
.PHONY: dist

FORCE:
	-rm $(distdir).tar.gz > /dev/null 2>&1
	-rm -rf $(distdir) > /dev/null 2>&1
.PHONY: FORCE

