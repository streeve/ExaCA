// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_GHOST_HPP
#define EXACA_GHOST_HPP

#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

class Halo {

  public:
    int mpi_rank;
    int mpi_size;
    int global_x;
    int global_y;
    int global_z;
    // 1D decomposition in Y: Each MPI rank has a subset consisting of of local_y cells, out of ny cells in Y
    int local_y;
    // Each MPI rank's subdomain is offset by offset_y cells from the lower bound of the domain in Y
    int offset_y;
    // Total cells in the local domain
    long int local_xyz;
    // Process IDs of neighboring MPI ranks on the grid: positive/negative Y directions are North/South
    int neighbor_north;
    int neighbor_south;
    // Whether or not each MPI rank is at a global domain boundary
    bool boundary_north;
    bool boundary_south;
    int buffer_size_x;
    int buffer_size_z;

    Buffer2D BufferSouthSend;
    Buffer2D BufferNorthSend;
    Buffer2D BufferSouthRecv;
    Buffer2D BufferNorthRecv;

    Halo(int id, int np, int nx, int ny, int nz);

    // Load data (GrainID, DOCenter, DiagonalLength) into ghost nodes if the given RankY is associated with a
    // 1D halo region
    KOKKOS_INLINE_FUNCTION void fill(const double GhostGID, const double GhostDOCX, const double GhostDOCY,
                                     const double GhostDOCZ, const double GhostDL, const int RankX, const int RankY,
                                     const int RankZ) const {

        if ((RankY == 1) && (!(boundary_south))) {
            int GNPosition = RankZ * buffer_size_x + RankX;
            BufferSouthSend(GNPosition, 0) = GhostGID;
            BufferSouthSend(GNPosition, 1) = GhostDOCX;
            BufferSouthSend(GNPosition, 2) = GhostDOCY;
            BufferSouthSend(GNPosition, 3) = GhostDOCZ;
            BufferSouthSend(GNPosition, 4) = GhostDL;
        }
        else if ((RankY == local_y - 2) && (!(boundary_north))) {
            int GNPosition = RankZ * buffer_size_x + RankX;
            BufferNorthSend(GNPosition, 0) = GhostGID;
            BufferNorthSend(GNPosition, 1) = GhostDOCX;
            BufferNorthSend(GNPosition, 2) = GhostDOCY;
            BufferNorthSend(GNPosition, 3) = GhostDOCZ;
            BufferNorthSend(GNPosition, 4) = GhostDL;
        }
    }

    void resizeBuffers(const int size_x, const int size_z);
    void gather(NList NeighborX, NList NeighborY, NList NeighborZ, ViewI CellType, ViewF DOCenter, ViewI GrainID,
                ViewF GrainUnitVector, ViewF DiagonalLength, ViewF CritDiagonalLength, int NGrainOrientations,
                int ZBound_Low);

    // Only used externally in one place.
    int calcLocalCellsY();
    int calcOffsetY();

  protected:
    void domainDecomposition();
    void initialDecomposition();
    void addGhosts();
};

#endif
