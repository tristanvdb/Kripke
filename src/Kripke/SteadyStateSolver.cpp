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

#include <Kripke/SteadyStateSolver.h>
#include <Kripke.h>
#include <Kripke/Core/Comm.h>
#include <Kripke/Kernel.h>
#include <Kripke/ParallelComm.h>
#include <Kripke/Core/PartitionSpace.h>
#include <Kripke/Timing.h>
#include <Kripke/SweepSolver.h>
#include <vector>
#include <stdio.h>

using namespace Kripke::Core;

void fields_summary(Kripke::Core::DataStore &data_store, size_t rank, size_t iter, size_t step) {
  data_store.getVariable<Kripke::Field_Flux_psi       >("psi"    ).summary(rank, iter, step);
  data_store.getVariable<Kripke::Field_Flux_rhs       >("rhs"    ).summary(rank, iter, step);
  data_store.getVariable<Kripke::Field_Moments_phi    >("phi"    ).summary(rank, iter, step);
  data_store.getVariable<Kripke::Field_Moments_phi_out>("phi_out").summary(rank, iter, step);
}

/**
  Run solver iterations.
*/
int Kripke::SteadyStateSolver (Kripke::Core::DataStore &data_store, size_t max_iter, bool block_jacobi, bool compute_errors)
{
  KRIPKE_TIMER(data_store, Solve);

  PartitionSpace &pspace = data_store.getVariable<PartitionSpace>("pspace");

  Kripke::Core::Comm const &comm = data_store.getVariable<Kripke::Core::Comm>("comm");
  if(comm.rank() == 0){
    printf("\n");
    printf("Steady State Solve\n");
    printf("==================\n\n");
  }

  // Intialize unknowns
  Kripke::Kernel::kConst(data_store.getVariable<Kripke::Field_Flux_psi>("psi"), 0.0);

  fields_summary(data_store, comm.rank(), 0, 0);

  // Loop over iterations
  double part_last = 0.0;
  for(size_t iter = 0;iter < max_iter;++ iter){

    /*
     * Compute the RHS:  rhs = LPlus*S*L*psi + Q
     */

    // Discrete to Moments transformation (phi = L*psi)
    Kripke::Kernel::kConst(data_store.getVariable<Field_Moments_phi>("phi"), 0.0);
    Kripke::Kernel::LTimes(data_store);

    fields_summary(data_store, comm.rank(), iter, 1);

    // Compute Scattering Source Term (psi_out = S*phi)
    Kripke::Kernel::kConst(data_store.getVariable<Kripke::Field_Moments_phi_out>("phi_out"), 0.0);
    Kripke::Kernel::scattering(data_store);

    fields_summary(data_store, comm.rank(), iter, 2);

    // Compute External Source Term (psi_out = psi_out + Q)
    Kripke::Kernel::source(data_store);

    fields_summary(data_store, comm.rank(), iter, 3);

    // Moments to Discrete transformation (rhs = LPlus*psi_out)
    Kripke::Kernel::kConst(data_store.getVariable<Kripke::Field_Flux_rhs>("rhs"), 0.0);
    Kripke::Kernel::LPlusTimes(data_store);

    fields_summary(data_store, comm.rank(), iter, 4);

    /*
     * Sweep (psi = Hinv*rhs)
     */
    {
      // Create a list of all groups
      int num_subdomains = pspace.getNumSubdomains(SPACE_PQR);
      std::vector<SdomId> sdom_list(num_subdomains);
      for(SdomId i{0};i < num_subdomains;++ i){
        sdom_list[*i] = i;
      }

      // Sweep everything
      Kripke::SweepSolver(data_store, sdom_list, block_jacobi);
    }

    fields_summary(data_store, comm.rank(), iter, 5);

    /*
     * Population edit and convergence test
     */
    double part = Kripke::Kernel::population(data_store);
    if(comm.rank() == 0){
      printf("  iter %zd: particle count=%e, change=%e\n", iter, part, (part-part_last)/part);
      fflush(stdout);
    }
    part_last = part;

  }

  fields_summary(data_store, comm.rank(), max_iter, 6);

  if (compute_errors)
    Kripke::Kernel::error_norms(data_store);

  if(comm.rank() == 0){
    printf("  Solver terminated\n");
  }

  return(0);
}




