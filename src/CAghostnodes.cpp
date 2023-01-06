// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAghostnodes.hpp"
#include "CAfunctions.hpp"
#include "CAupdate.hpp"
#include "mpi.h"

#include <cmath>
#include <vector>

//*****************************************************************************/
// Constructor
Halo::Halo(int id, int np, int nx, int ny, int nz)
    : mpi_rank(id)
    , mpi_size(np)
    , global_x(nx)
    , global_y(ny)
    , global_z(nz)
    , BufferSouthSend(Buffer2D("BufferSouthSend", 0, 0))
    , BufferNorthSend(Buffer2D("BufferNorthSend", 0, 0))
    , BufferSouthRecv(Buffer2D("BufferSouthRecv", 0, 0))
    , BufferNorthRecv(Buffer2D("BufferNorthRecv", 0, 0)) {
    buffer_size_x = 0;
    buffer_size_z = 0;
    domainDecomposition();
}

// Decompose the domain into subdomains on each MPI rank: Calculate local_y and offset_y for each rank, where
// each subdomain contains "local_y" in Y, offset from the full domain origin by "offset_y" cells in Y
void Halo::domainDecomposition() {

    // Compare total MPI ranks to total Y cells.
    if (mpi_size > global_y)
        throw std::runtime_error("Error: Cannot run with more MPI ranks than cells in Y (decomposition direction).");

    // Determine which subdomains are at which locations on the grid relative to the others
    initialDecomposition();
    // Determine, for each MPI process id, the local grid size in x and y (and the offsets in x and y relative to
    // the overall simulation domain)
    offset_y = calcOffsetY();
    local_y = calcLocalCellsY();

    // Add ghost nodes at subdomain overlaps
    addGhosts();

    local_xyz = global_x * local_y * global_z; // Number of cells on this MPI rank
}

//*****************************************************************************/
// Determine the mapping of processors to grid data
void Halo::initialDecomposition() {

    if (mpi_size > 1) {
        neighbor_south = mpi_rank - 1;
        neighbor_north = mpi_rank + 1;
        if (mpi_rank == 0)
            neighbor_south = MPI_PROC_NULL;
        if (mpi_rank == mpi_size - 1)
            neighbor_north = MPI_PROC_NULL;
    }
    else {
        // No MPI communication
        neighbor_north = MPI_PROC_NULL;
        neighbor_south = MPI_PROC_NULL;
    }

    // Store whether each MPI rank is on each boundary or not
    if (neighbor_north == MPI_PROC_NULL)
        boundary_north = true;
    else
        boundary_north = false;
    if (neighbor_south == MPI_PROC_NULL)
        boundary_south = true;
    else
        boundary_south = false;
}

//*****************************************************************************/
int Halo::calcLocalCellsY() {
    int YRemoteMPSlices = 0;
    int YSlicesPerP = global_y / mpi_size;
    int YRemainder = global_y % mpi_size;
    if (YRemainder == 0) {
        YRemoteMPSlices = YSlicesPerP;
    }
    else {
        if (YRemainder > mpi_rank) {
            YRemoteMPSlices = YSlicesPerP + 1;
        }
        else {
            YRemoteMPSlices = YSlicesPerP;
        }
    }
    return YRemoteMPSlices;
}

//*****************************************************************************/
int Halo::calcOffsetY() {
    int RemoteYOffset = 0;
    int YSlicesPerP = global_y / mpi_size;
    int YRemainder = global_y % mpi_size;
    if (YRemainder == 0) {
        RemoteYOffset = mpi_rank * YSlicesPerP;
    }
    else {
        if (YRemainder > mpi_rank) {
            RemoteYOffset = mpi_rank * (YSlicesPerP + 1);
        }
        else {
            RemoteYOffset =
                (YSlicesPerP + 1) * (YRemainder - 1) + YSlicesPerP + 1 + (mpi_rank - YRemainder) * YSlicesPerP;
        }
    }
    return RemoteYOffset;
}

// Add ghost nodes to the appropriate subdomains (added where the subdomains overlap, but not at edges of physical
// domain)
void Halo::addGhosts() {

    // Add halo regions in Y direction if this subdomain borders subdomains on other processors
    // If only 1 rank in the y direction, no halo regions - subdomain is coincident with overall simulation domain
    // If multiple ranks in the y direction, either 1 halo region (borders another rank's subdomain in either the +y
    // or -y direction) or 2 halo regions (if it borders other rank's subdomains in both the +y and -y directions)
    if (neighbor_north != MPI_PROC_NULL)
        local_y++;
    if (neighbor_south != MPI_PROC_NULL) {
        local_y++;
        // Also adjust subdomain offset, as these ghost nodes were added on the -y side of the subdomain
        offset_y--;
    }
}

void Halo::resizeBuffers(const int size_x, const int size_z) {
    buffer_size_x = size_x;
    buffer_size_z = size_z;
    auto buffer_size = buffer_size_x * buffer_size_z;

    // Send/recv buffers for ghost node data should be initialized with zeros
    Kokkos::resize(BufferSouthSend, buffer_size, 5);
    Kokkos::resize(BufferNorthSend, buffer_size, 5);
    Kokkos::resize(BufferSouthRecv, buffer_size, 5);
    Kokkos::resize(BufferNorthRecv, buffer_size, 5);
}

// 1D domain decomposition: update ghost nodes with new cell data from Nucleation and CellCapture routines
void Halo::gather(NList NeighborX, NList NeighborY, NList NeighborZ, ViewI CellType, ViewF DOCenter, ViewI GrainID,
                  ViewF GrainUnitVector, ViewF DiagonalLength, ViewF CritDiagonalLength, int NGrainOrientations,
                  int ZBound_Low) {

    std::vector<MPI_Request> SendRequests(2, MPI_REQUEST_NULL);
    std::vector<MPI_Request> RecvRequests(2, MPI_REQUEST_NULL);

    // Send data to each other rank (MPI_Isend)
    auto buffer_size = 5 * buffer_size_x * buffer_size_z;
    MPI_Isend(BufferSouthSend.data(), buffer_size, MPI_DOUBLE, neighbor_south, 0, MPI_COMM_WORLD, &SendRequests[0]);
    MPI_Isend(BufferNorthSend.data(), buffer_size, MPI_DOUBLE, neighbor_north, 0, MPI_COMM_WORLD, &SendRequests[1]);

    // Receive buffers for all neighbors (MPI_Irecv)
    MPI_Irecv(BufferSouthRecv.data(), buffer_size, MPI_DOUBLE, neighbor_south, 0, MPI_COMM_WORLD, &RecvRequests[0]);
    MPI_Irecv(BufferNorthRecv.data(), buffer_size, MPI_DOUBLE, neighbor_north, 0, MPI_COMM_WORLD, &RecvRequests[1]);

    // unpack in any order
    bool unpack_complete = false;
    while (!unpack_complete) {
        // Get the next buffer to unpack from rank "unpack_index"
        int unpack_index = MPI_UNDEFINED;
        MPI_Waitany(2, RecvRequests.data(), &unpack_index, MPI_STATUS_IGNORE);
        // If there are no more buffers to unpack, leave the while loop
        if (MPI_UNDEFINED == unpack_index) {
            unpack_complete = true;
        }
        // Otherwise unpack the next buffer.
        else {
            int RecvBufSize = buffer_size_x * buffer_size_z;
            Kokkos::parallel_for(
                "BufferUnpack", RecvBufSize, KOKKOS_LAMBDA(const int &BufPosition) {
                    int RankX, RankY, RankZ, NewGrainID;
                    long int CellLocation;
                    double DOCenterX, DOCenterY, DOCenterZ, NewDiagonalLength;
                    bool Place = false;
                    RankZ = BufPosition / buffer_size_x;
                    RankX = BufPosition % buffer_size_x;
                    // Which rank was the data received from?
                    if ((unpack_index == 0) && (neighbor_south != MPI_PROC_NULL)) {
                        // Data receieved from South
                        RankY = 0;
                        CellLocation = RankZ * global_x * local_y + local_y * RankX + RankY;
                        int GlobalCellLocation = CellLocation + ZBound_Low * global_x * local_y;
                        if ((BufferSouthRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferSouthRecv(BufPosition, 0));
                            DOCenterX = BufferSouthRecv(BufPosition, 1);
                            DOCenterY = BufferSouthRecv(BufPosition, 2);
                            DOCenterZ = BufferSouthRecv(BufPosition, 3);
                            NewDiagonalLength = BufferSouthRecv(BufPosition, 4);
                        }
                    }
                    else if ((unpack_index == 1) && (neighbor_north != MPI_PROC_NULL)) {
                        // Data received from North
                        RankY = local_y - 1;
                        CellLocation = RankZ * global_x * local_y + local_y * RankX + RankY;
                        int GlobalCellLocation = CellLocation + ZBound_Low * global_x * local_y;
                        if ((BufferNorthRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferNorthRecv(BufPosition, 0));
                            DOCenterX = BufferNorthRecv(BufPosition, 1);
                            DOCenterY = BufferNorthRecv(BufPosition, 2);
                            DOCenterZ = BufferNorthRecv(BufPosition, 3);
                            NewDiagonalLength = BufferNorthRecv(BufPosition, 4);
                        }
                    }
                    if (Place) {
                        int GlobalZ = RankZ + ZBound_Low;
                        int GlobalCellLocation = GlobalZ * global_x * local_y + RankX * local_y + RankY;

                        // Update this ghost node cell's information with data from other rank
                        GrainID(GlobalCellLocation) = NewGrainID;
                        DOCenter((long int)(3) * CellLocation) = static_cast<float>(DOCenterX);
                        DOCenter((long int)(3) * CellLocation + (long int)(1)) = static_cast<float>(DOCenterY);
                        DOCenter((long int)(3) * CellLocation + (long int)(2)) = static_cast<float>(DOCenterZ);
                        int MyOrientation = getGrainOrientation(GrainID(GlobalCellLocation), NGrainOrientations);
                        DiagonalLength(CellLocation) = static_cast<float>(NewDiagonalLength);
                        // Global coordinates of cell center
                        double xp = RankX + 0.5;
                        double yp = RankY + offset_y + 0.5;
                        double zp = GlobalZ + 0.5;
                        // Calculate critical values at which this active cell leads to the activation of a neighboring
                        // liquid cell
                        calcCritDiagonalLength(CellLocation, xp, yp, zp, DOCenterX, DOCenterY, DOCenterZ, NeighborX,
                                               NeighborY, NeighborZ, MyOrientation, GrainUnitVector,
                                               CritDiagonalLength);
                        CellType(GlobalCellLocation) = Active;
                    }
                });
        }
    }

    // Wait on send requests
    MPI_Waitall(2, SendRequests.data(), MPI_STATUSES_IGNORE);
    Kokkos::fence();
}
