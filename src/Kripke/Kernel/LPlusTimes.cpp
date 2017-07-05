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

#include<Kripke/Kernel.h>
#include<Kripke/Grid.h>
#include<Kripke/SubTVec.h>
#include<Kripke/Timing.h>

void Kripke::Kernel::LPlusTimes(Kripke::DataStore &data_store)
{
  KRIPKE_TIMER(data_store, LPlusTimes);

  Grid_Data *grid_data = &data_store.getVariable<Grid_Data>("grid_data");

  // Outer parameters
  int num_moments = grid_data->total_num_moments;

  // Loop over Subdomains
  int num_subdomains = grid_data->subdomains.size();
  for (int sdom_id = 0; sdom_id < num_subdomains; ++ sdom_id){
    Subdomain &sdom = grid_data->subdomains[sdom_id];

    // Get dimensioning
    int num_zones = sdom.num_zones;
    int num_local_groups = sdom.num_groups;
    int num_groups = sdom.phi_out->groups;
    int group0 = sdom.group0;
    int num_local_directions = sdom.num_directions;
    int num_groups_zones = num_local_groups*num_zones;
    
    // Zero RHS
    sdom.rhs->clear(0.0);

    // Get pointers
    double const * KRESTRICT phi_out = sdom.phi_out->ptr() + group0*num_zones;
    double const * KRESTRICT ell_plus = sdom.ell_plus->ptr();
    double       * KRESTRICT rhs = sdom.rhs->ptr();

    for (int d = 0; d < num_local_directions; d++) {      
      double       * KRESTRICT rhs_d = rhs + d*num_groups_zones;
      double const * KRESTRICT ell_plus_d = ell_plus + d*num_moments;
      
      for(int nm = 0;nm < num_moments;++nm){
        double const ell_plus_d_nm = ell_plus_d[nm];
        double const * KRESTRICT phi_out_nm = phi_out + nm*num_groups*num_zones;

        for(int gz = 0;gz < num_groups_zones; ++ gz){
          rhs_d[gz] += ell_plus_d_nm * phi_out_nm[gz];
        }
      }
    }
  }
}