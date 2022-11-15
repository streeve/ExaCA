// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_TYPES_HPP
#define EXACA_TYPES_HPP

#include <Kokkos_Core.hpp>

enum TypeNames {
    Wall = 0,
    Solid = 1,
    Active = 2,
    TemporaryUpdate = 3,
    TemporaryInit = 4,
    Liquid = 5,
    TempSolid = 6,
    FutureActive = 7
};

// Use Kokkos::DefaultExecutionSpace
typedef Kokkos::View<float *> ViewF;
typedef Kokkos::View<int *> ViewI;
typedef Kokkos::View<int **> ViewI2D;
typedef Kokkos::View<int *, Kokkos::MemoryTraits<Kokkos::Atomic>> View_a;
typedef Kokkos::View<float ***> ViewF3D;

using exe_space = Kokkos::DefaultExecutionSpace::execution_space;
using device_memory_space = Kokkos::DefaultExecutionSpace::memory_space;
typedef typename exe_space::array_layout layout;
typedef Kokkos::View<float *, layout, Kokkos::HostSpace> ViewF_H;
typedef Kokkos::View<float **, layout, Kokkos::HostSpace> ViewF2D_H;
typedef Kokkos::View<float ***, layout, Kokkos::HostSpace> ViewF3D_H;
typedef Kokkos::View<int *, layout, Kokkos::HostSpace> ViewI_H;
typedef Kokkos::View<int **, layout, Kokkos::HostSpace> ViewI2D_H;
typedef Kokkos::View<int ***, layout, Kokkos::HostSpace> ViewI3D_H;

typedef Kokkos::Array<int, 26> NList;

namespace Impl {

template <typename... Types>
struct MemberTypes {
    static constexpr std::size_t size = sizeof...(Types);
};

template <std::size_t M, typename T, typename... Types>
struct MemberTypeAtIndexImpl;

template <typename T, typename... Types>
struct MemberTypeAtIndexImpl<0, T, Types...> {
    using type = T;
};

template <std::size_t M, typename T, typename... Types>
struct MemberTypeAtIndexImpl {
    using type = typename MemberTypeAtIndexImpl<M - 1, Types...>::type;
};

//! Get the type of the member at a given index.
template <std::size_t M, typename... Types>
struct MemberTypeAtIndex;

//! Get the type of the member at a given index.
template <std::size_t M, typename... Types>
struct MemberTypeAtIndex<M, MemberTypes<Types...>> {
    //! Member type.
    using type = typename MemberTypeAtIndexImpl<M, Types...>::type;
};

} // namespace Impl

// Forward declaration.
template <typename Types>
struct CommData;

template <typename... Types>
struct CommData<Impl::MemberTypes<Types...>> {};

typedef Impl::MemberTypes<int, float, float, float, float> CommMemberTypes;
typedef Kokkos::View<CommData<CommMemberTypes> *> ViewBuffer;

#endif
