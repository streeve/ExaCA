// Copyright Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_PRINT_ADIOS2_HPP
#define EXACA_PRINT_ADIOS2_HPP

#include "CAprint.hpp"

#include "mpi.h"

#include <Kokkos_Core.hpp>
#include <adios2.h>

struct PrintADIOS2 : public Print {
    using base = Print;
    using base::base;

    void printIntralayer(const int id, const int np, const int layernumber, const double deltat, const int cycle,
                         const Grid &grid, CellData<MemorySpace> &celldata, Temperature<MemorySpace> &temperature,
                         Interface<MemorySpace> &interface, Orientation<MemorySpace> &orientation) {
        adios2::fstream oStream(time_series_filename, adios2::fstream::out, MPI_COMM_WORLD);
    }

    template <typename Print3DViewType>
    void printViewData(const int id, std::ofstream &output_fstream, const Grid &grid, const bool current_layer_only,
                       std::string data_label, std::string var_name_label, Print3DViewType view_data_whole_domain) {}
};

#endif
