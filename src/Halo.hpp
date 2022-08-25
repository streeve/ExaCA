// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_HALO_HPP
#define EXACA_HALO_HPP

#include "CAtypes.hpp"
#include <Kokkos_Core.hpp>

#include "mpi.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Class specifying ExaCA buffers for halo exchange operations

class Halo {

  public:

    ViewBuffer BufferSouthSend, BufferSouthRecv, BufferNorthSend, BufferNorthRecv;
    //-Constructor
    Halo(const int BufSizeX, const int BufSizeZ);
    // Reset to zeros
    void reset(const int Buffersize);
    // Destructor
    ~Halo();
};

#endif
