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
TEST(TestGPV, testGetProbabilityVector) {
  std::string input = "CACCAAAAAAGGAG";
  const double TOL = 0.0001;
  double rightAnswer[] = {0.0323673,0,0.840191,0.8169,0,0,0,0,0,0,0.817312,0.840378,0,0.0317688};
  std::vector<double> pv =
    RNAUtils::getBasePairProbabilityVector(input);
  for (int i=0; i<14; i++)
    EXPECT_NEAR(rightAnswer[i], pv[i],TOL);
}

