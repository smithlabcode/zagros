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
	@rm -rf $(RBP)/bin
	@rm -rf $(RBP)/lib
	@rm -rf $(RBP)/include
	@rm -rf $(RBP)/developmentDocs
.PHONY: distclean
