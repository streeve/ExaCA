// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAfunctions.hpp"

#include "mpi.h"

#include <cmath>
#include <iostream>

/*************************** FUNCTIONS CALLED THROUGH MAIN SUBROUTINES ***************************/

//*****************************************************************************/
int FindItBounds(int RankX, int RankY, int MyXSlices, int MyYSlices) {
    int ItBounds;
    // If X and Y coordinates are not on edges, Case 0: iteratation over neighbors 0-25 possible
    // If Y coordinate is on lower edge, Case 1: iteration over only neighbors 9-25 possible
    // If Y coordinate is on upper edge, Case 2: iteration over only neighbors 0-16 possible
    // If X coordinate is on lower edge, Case 3: iteration over only neighbors
    // 0,1,3,4,6,8,9,10,11,13,15,17,18,20,21,22,24 If X coordinate is on upper edge, Case 4: iteration over only
    // neighbors 0,2,3,4,5,7,9,10,12,14,16,17,19,20,21,23,25 If X/Y coordinates are on lower edge, Case 5: iteration
    // over only neighbors 9,10,11,13,15,17,18,20,21,22,24 If X coordinate is on upper edge/Y on lower edge, Case 6: If
    // X coordinate is on lower edge/Y on upper edge, Case 7: If X/Y coordinates are on upper edge, Case 8:
    if (RankY == 0) {
        if (RankX == 0) {
            ItBounds = 5;
        }
        else if (RankX == MyXSlices - 1) {
            ItBounds = 6;
        }
        else {
            ItBounds = 1;
        }
    }
    else if (RankY == MyYSlices - 1) {
        if (RankX == 0) {
            ItBounds = 7;
        }
        else if (RankX == MyXSlices - 1) {
            ItBounds = 8;
        }
        else {
            ItBounds = 2;
        }
    }
    else {
        if (RankX == 0) {
            ItBounds = 3;
        }
        else if (RankX == MyXSlices - 1) {
            ItBounds = 4;
        }
        else {
            ItBounds = 0;
        }
    }
    return ItBounds;
}

// Create a view of size "NumberOfOrientation" of the misorientation of each possible grain orientation with the X, Y,
// or Z directions (dir = 0, 1, or 2, respectively)
ViewF_H MisorientationCalc(int NumberOfOrientations, ViewF_H GrainUnitVector, int dir) {

    ViewF_H GrainMisorientation(Kokkos::ViewAllocateWithoutInitializing("GrainMisorientation"), NumberOfOrientations);
    // Find the smallest possible misorientation between the specified direction, and this grain orientations' 6
    // possible 001 directions (where 62.7 degrees is the largest possible misorientation between two 001 directions
    // for a cubic crystal system)
    for (int n = 0; n < NumberOfOrientations; n++) {
        float MisorientationAngleMin = 62.7;
        for (int ll = 0; ll < 3; ll++) {
            float Misorientation = std::abs((180 / M_PI) * acos(GrainUnitVector(9 * n + 3 * ll + dir)));
            if (Misorientation < MisorientationAngleMin) {
                MisorientationAngleMin = Misorientation;
            }
        }
        GrainMisorientation(n) = MisorientationAngleMin;
    }
    return GrainMisorientation;
}
