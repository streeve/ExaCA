// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_GHOST_HPP
#define EXACA_GHOST_HPP

#include "CAtypes.hpp"
#include "Halo.hpp"

#include <Kokkos_Core.hpp>

// Load data (GrainID, DOCenter, DiagonalLength) into ghost nodes if the given RankY is associated with a 1D halo region
KOKKOS_INLINE_FUNCTION void loadghostnodes(const int GhostGID, const float GhostDOCX, const float GhostDOCY,
                                           const float GhostDOCZ, const float GhostDL, const int BufSizeX,
                                           const int MyYSlices, const int RankX, const int RankY, const int RankZ,
                                           const bool AtNorthBoundary, const bool AtSouthBoundary,
                                           ViewBuffer BufferSouthSend, ViewBuffer BufferNorthSend) {

    if ((RankY == 1) && (!(AtSouthBoundary))) {
        int GNPosition = RankZ * BufSizeX + RankX;
        BufferSouthSend(GNPosition).GrainID = GhostGID;
        BufferSouthSend(GNPosition).DOCenterX = GhostDOCX;
        BufferSouthSend(GNPosition).DOCenterY = GhostDOCY;
        BufferSouthSend(GNPosition).DOCenterZ = GhostDOCZ;
        BufferSouthSend(GNPosition).DiagonalLength = GhostDL;
    }
    else if ((RankY == MyYSlices - 2) && (!(AtNorthBoundary))) {
        int GNPosition = RankZ * BufSizeX + RankX;
        BufferNorthSend(GNPosition).GrainID = GhostGID;
        BufferNorthSend(GNPosition).DOCenterX = GhostDOCX;
        BufferNorthSend(GNPosition).DOCenterY = GhostDOCY;
        BufferNorthSend(GNPosition).DOCenterZ = GhostDOCZ;
        BufferNorthSend(GNPosition).DiagonalLength = GhostDL;
    }
}
void GhostNodes1D(int, int, int NeighborRank_North, int NeighborRank_South, int nx, int MyYSlices, int MyYOffset,
                  NList NeighborX, NList NeighborY, NList NeighborZ, ViewI CellType, ViewF DOCenter, ViewI GrainID,
                  ViewF GrainUnitVector, ViewF DiagonalLength, ViewF CritDiagonalLength, int NGrainOrientations,
                  Halo &halo, int BufSizeX, int BufSizeZ, int ZBound_Low);

#endif
