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
TEST(TestWNSS, testWeightedNumberofSecondaryStructres) {
  std::string input = "CACCAAAAAAGGAG";
  const double TOL = 0.0001;
  double rightAnswer[14][14] =  {{1,1,1,1,1,1,1,1,1,1,1.81731,2.52496,2.52496,2.52469},
                                 {0,1,1,1,1,1,1,1,1,1,1.81686,2.52424,2.52424,2.52398},
                                 {0,0,1,1,1,1,1,1,1,1,1.81686,2.83973,2.83973,2.83924},
                                 {0,0,0,1,1,1,1,1,1,1,1.8157,1.81561,1.81561,1.81507},
                                 {0,0,0,0,1,1,1,1,1,1,1,1,1,1},
                                 {0,0,0,0,0,1,1,1,1,1,1,1,1,1},
                                 {0,0,0,0,0,0,1,1,1,1,1,1,1,1},
                                 {0,0,0,0,0,0,0,1,1,1,1,1,1,1},
                                 {0,0,0,0,0,0,0,0,1,1,1,1,1,1},
                                 {0,0,0,0,0,0,0,0,0,1,1,1,1,1},
                                 {0,0,0,0,0,0,0,0,0,0,1,1,1,1},
                                 {0,0,0,0,0,0,0,0,0,0,0,1,1,1},
                                 {0,0,0,0,0,0,0,0,0,0,0,0,1,1},
                                 {0,0,0,0,0,0,0,0,0,0,0,0,0,1}};

  std::vector<std::vector<double> > ns =
    RNAUtils::calculateWeightedNumberOfStructures(input);
  for (int i=0; i<14; i++)
    for (int j=0; j<14; j++)
      EXPECT_NEAR(rightAnswer[i][j], ns[i][j], TOL);
}

