/*
 * NOTICE
 *
 * This work was produced at the Lawrence Livermore National Laboratory (LLNL)
 * under contract no. DE-AC-52-07NA27344 (Contract 44) between the U.S.
 * Department of Energy (DOE) and Lawrence Livermore National Security, LLC
 * (LLNS) for the operation of LLNL. The rights of the Federal Government are
 * reserved under Contract 44.
 *
 * DISCLAIMER
 *
 * This work was prepared as an account of work sponsored by an agency of the
 * United States Government. Neither the United States Government nor Lawrence
 * Livermore National Security, LLC nor any of their employees, makes any
 * warranty, express or implied, or assumes any liability or responsibility
 * for the accuracy, completeness, or usefulness of any information, apparatus,
 * product, or process disclosed, or represents that its use would not infringe
 * privately-owned rights. Reference herein to any specific commercial products,
 * process, or service by trade name, trademark, manufacturer or otherwise does
 * not necessarily constitute or imply its endorsement, recommendation, or
 * favoring by the United States Government or Lawrence Livermore National
 * Security, LLC. The views and opinions of authors expressed herein do not
 * necessarily state or reflect those of the United States Government or
 * Lawrence Livermore National Security, LLC, and shall not be used for
 * advertising or product endorsement purposes.
 *
 * NOTIFICATION OF COMMERCIAL USE
 *
 * Commercialization of this product is prohibited without notifying the
 * Department of Energy (DOE) or Lawrence Livermore National Security.
 */

#ifndef KRIPKE_ARCH_SWEEPSUBDOMAINS
#define KRIPKE_ARCH_SWEEPSUBDOMAINS

#include <Kripke.h>
#include <Kripke/VarTypes.h>

namespace Kripke {
namespace Arch {

template<typename AL>
struct Policy_SweepSubdomains;

template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_Sequential, LayoutT_DGZ>> {
  using ExecPolicy =
    KernelPolicy<
      For<0, loop_exec,  // direction
        For<1, loop_exec, // group
          For<2, loop_exec, // k
            For<3, loop_exec, // j
              For<4, loop_exec, // i
                Lambda<0>
              >
            >
          >
        >
      >
    >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_Sequential, LayoutT_DZG>> {
  using ExecPolicy =
    KernelPolicy<
      For<0, loop_exec,  // direction
        For<2, loop_exec, // k
          For<3, loop_exec, // j
            For<4, loop_exec, // i
              For<1, loop_exec, // group
                Lambda<0>
              >
            >
          >
        >
      >
    >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_Sequential, LayoutT_GDZ>> {
  using ExecPolicy =
    KernelPolicy<
      For<1, loop_exec, // group
        For<0, loop_exec,  // direction
          For<2, loop_exec, // k
            For<3, loop_exec, // j
              For<4, loop_exec, // i
                Lambda<0>
              >
            >
          >
        >
      >
    >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_Sequential, LayoutT_GZD>> {
  using ExecPolicy =
    KernelPolicy<
      For<1, loop_exec, // group
        For<2, loop_exec, // k
          For<3, loop_exec, // j
            For<4, loop_exec, // i
              For<0, loop_exec,  // direction
                Lambda<0>
              >
            >
          >
        >
      >
    >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_Sequential, LayoutT_ZDG>> {
  using ExecPolicy =
    KernelPolicy<
      For<2, loop_exec, // k
        For<3, loop_exec, // j
          For<4, loop_exec, // i
            For<0, loop_exec,  // direction
              For<1, loop_exec, // group
                Lambda<0>
              >
            >
          >
        >
      >
    >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_Sequential, LayoutT_ZGD>> {
  using ExecPolicy =
    KernelPolicy<
      For<2, loop_exec, // k
        For<3, loop_exec, // j
          For<4, loop_exec, // i
            For<1, loop_exec, // group
              For<0, loop_exec,  // direction
                Lambda<0>
              >
            >
          >
        >
      >
    >;
};




#ifdef KRIPKE_USE_OPENMP


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_OpenMP, LayoutT_DGZ>> {


  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<0,1>, // direction, group
        For<2, loop_exec, // k
          For<3, loop_exec, // j
            For<4, loop_exec, // i
              Lambda<0>
            >
          >
        >
      >
    >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_OpenMP, LayoutT_DZG>> {
  using ExecPolicy =
    KernelPolicy<
      For<0, omp_parallel_for_exec,  // direction
        For<2, loop_exec, // k
          For<3, loop_exec, // j
            For<4, loop_exec, // i
              For<1, loop_exec, // group
                Lambda<0>
              >
            >
          >
        >
      >
    >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_OpenMP, LayoutT_GDZ>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<1,0>, // group, direction
        For<2, loop_exec, // k
          For<3, loop_exec, // j
            For<4, loop_exec, // i
              Lambda<0>
            >
          >
        >
      >
    >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_OpenMP, LayoutT_GZD>> {
  using ExecPolicy =
    KernelPolicy<
      For<1, omp_parallel_for_exec, // group
        For<2, loop_exec, // k
          For<3, loop_exec, // j
            For<4, loop_exec, // i
              For<0, loop_exec,  // direction
                Lambda<0>
              >
            >
          >
        >
      >
    >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_OpenMP, LayoutT_ZDG>> {
  using ExecPolicy =
    KernelPolicy<
      Hyperplane<2, seq_exec, ArgList<3,4>, omp_parallel_collapse_exec,
        For<0, loop_exec,  // direction
          For<1, loop_exec, // group
            Lambda<0>
          >
        >
      >
    >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_OpenMP, LayoutT_ZGD>> {
  using ExecPolicy =
    KernelPolicy<
      Hyperplane<2, seq_exec, ArgList<3,4>, omp_parallel_collapse_exec,
        For<1, loop_exec, // group
          For<0, loop_exec,  // direction
            Lambda<0>
          >
        >
      >
    >;
};

#endif // KRIPKE_USE_OPENMP


#ifdef KRIPKE_USE_CUDA
template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_CUDA, LayoutT_DGZ>> {
  using ExecPolicy =
          KernelPolicy<
            CudaKernel<
              For<0, cuda_block_exec,
                For<1, cuda_block_exec,

                      Hyperplane<
                        2, cuda_seq_syncthreads_exec,
                        ArgList<3, 4>, cuda_thread_exec,

                        Lambda<0>
                      >
                >
              >
            >
          >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_CUDA, LayoutT_DZG>> {
    using ExecPolicy =
            KernelPolicy<
              CudaKernel<
                For<0, cuda_block_exec,
                  For<1, cuda_block_exec,

                        Hyperplane<
                          2, cuda_seq_syncthreads_exec,
                          ArgList<3, 4>, cuda_thread_exec,

                          Lambda<0>
                        >
                  >
                >
              >
            >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_CUDA, LayoutT_GDZ>> {
    using ExecPolicy =
            KernelPolicy<
              CudaKernel<
                For<0, cuda_block_exec,
                  For<1, cuda_block_exec,

                        Hyperplane<
                          2, cuda_seq_syncthreads_exec,
                          ArgList<3, 4>, cuda_thread_exec,

                          Lambda<0>
                        >
                  >
                >
              >
            >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_CUDA, LayoutT_GZD>> {
    using ExecPolicy =
            KernelPolicy<
              CudaKernel<
                For<0, cuda_block_exec,
                  For<1, cuda_block_exec,

                        Hyperplane<
                          2, cuda_seq_syncthreads_exec,
                          ArgList<3, 4>, cuda_thread_exec,

                          Lambda<0>
                        >
                  >
                >
              >
            >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_CUDA, LayoutT_ZDG>> {
    using ExecPolicy =
            KernelPolicy<
              CudaKernel<
                For<0, cuda_block_exec,
                  For<1, cuda_block_exec,

                        Hyperplane<
                          2, cuda_seq_syncthreads_exec,
                          ArgList<3, 4>, cuda_thread_exec,

                          Lambda<0>
                        >
                  >
                >
              >
            >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_CUDA, LayoutT_ZGD>> {
    using ExecPolicy =
            KernelPolicy<
              CudaKernel<
                For<0, cuda_block_exec,
                  For<1, cuda_block_exec,

                        Hyperplane<
                          2, cuda_seq_syncthreads_exec,
                          ArgList<3, 4>, cuda_thread_exec,

                          Lambda<0>
                        >
                  >
                >
              >
            >;
};
#endif // KRIPKE_USE_CUDA

}
}

#endif
