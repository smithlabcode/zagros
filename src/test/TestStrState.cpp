/**
  \file TestZTP.cpp
  \brief This source file defines a set of unit tests for the ZTP distribution
         class.

  \authors Philip J. Uren, Andrew D. Smith

  \section copyright Copyright Details
  Copyright (C) 2011
  University of Southern California,
  Philip J. Uren, Andrew D. Smith

  \section license License Details
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  \section bugs Known Bugs

  \section history Revision History
****/

#include "gtest/gtest.h"
#include "Model.hpp"
#include "part_func.hpp"
#include "RNA_Utils.hpp"

/**
 * \brief test the zero-truncated poisson PDF with a few sample values with
 *        known answers
 */
TEST(TestStrState, testCalculateUnpairedState) {
  std::string input = "CACCAAAAAAGGAG";
  const double TOL = 0.0001;
  double rightAnswer[] = {0.967633,1,0.159809,0.1831,1,1,1,1,1,1,0.182688,0.159622,1,0.968231};

  std::vector< std::vector<double> > tmp =
    RNAUtils::calculateWeightedNumberOfStructures(input);

  std::vector<double> alpha =
    RNAUtils::calculateUnpairedState(input, tmp[0][input.length()-1]);
  for (int i=0; i<14; i++)
    EXPECT_NEAR(rightAnswer[i], alpha[i], TOL);
}

/**
 * \brief test the zero-truncated poisson PDF with a few sample values with
 *        known answers
 */
TEST(TestStrState, testCalculatePairedState) {
  std::string input = "CACCAAAAAAGGAG";
  const double TOL = 0.0001;
  double rightAnswer[] = {0.0385151,0,0.999778,0.972062,0,0,0,0,0,0,0.972552,1,0,0.0378029};

  std::vector< std::vector<double> > tmp =
    RNAUtils::calculateWeightedNumberOfStructures(input);

  std::vector<double> gamma =
    RNAUtils::calculatePairedState(input, tmp[0][input.length()-1]);
  for (int i=0; i<14; i++)
    EXPECT_NEAR(rightAnswer[i], gamma[i], TOL);
}

