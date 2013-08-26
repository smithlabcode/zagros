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
TEST(TestBOIM, testBaseOccurrenceIndicatorMatrix) {
  std::string input = "CACCAAAAAAGGAG";
  double rightAnswer[4][14] = { {0,1,0,0,1,1,1,1,1,1,0,0,1,0},
				{1,0,1,1,0,0,0,0,0,0,0,0,0,0},
				{0,0,0,0,0,0,0,0,0,0,1,1,0,1},
				{0,0,0,0,0,0,0,0,0,0,0,0,0,0} };

  std::vector<std::vector<size_t> > nu =
    RNAUtils::baseOccurrenceIndicatorMatrix(input);
  for (int i=0; i<4; i++)
    for (int j=0; j<14; j++)
      EXPECT_EQ(rightAnswer[i][j], nu[i][j]);
}

