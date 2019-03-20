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

#include <Kripke/Generate.h>

#include <Kripke/Core/Comm.h>
#include <Kripke/Core/Field.h>
#include <Kripke/Kernel.h>
#include <Kripke/Core/PartitionSpace.h>
#include <Kripke/Core/Set.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>

using namespace Kripke;
using namespace Kripke::Core;

#ifdef KRIPKE_USE_ZFP
double double_zfp_psi::zfp_rate{64};
size_t double_zfp_psi::cached_zfp_blocks{1024};
double double_zfp_phi::zfp_rate{64};
size_t double_zfp_phi::cached_zfp_blocks{1024};
double double_zfp_exclude_1_fast::zfp_rate{64};
size_t double_zfp_exclude_1_fast::cached_zfp_blocks{1024};
#endif

void Kripke::Generate::generateData(Kripke::Core::DataStore &data_store,
    InputVariables const &input_vars)
{

  PartitionSpace &pspace = data_store.getVariable<PartitionSpace>("pspace");


  // Create a set to span angular the flux
  Set const &dir_set   = data_store.getVariable<Set>("Set/Direction");
  Set const &group_set = data_store.getVariable<Set>("Set/Group");
  Set const &zone_set  = data_store.getVariable<Set>("Set/Zone");
  ProductSet<3> *flux_set = new ProductSet<3>(pspace, SPACE_PQR,
      dir_set, group_set, zone_set);

  data_store.addVariable("Set/Flux", flux_set);

  ArchLayoutV al_v = data_store.getVariable<ArchLayout>("al").al_v;
  
  // Create Solution and RHS fields
  createField<Field_Flux>(data_store, "psi", al_v, *flux_set);
  createField<Field_Flux>(data_store, "rhs", al_v, *flux_set);


  // Create a set to span moments of the angular flux
  Set const &moment_set   = data_store.getVariable<Set>("Set/Moment");
  ProductSet<3> *fluxmoment_set = new ProductSet<3>(pspace, SPACE_PR,
        moment_set, group_set, zone_set);

  data_store.addVariable("Set/FluxMoment", fluxmoment_set);


  // Create flux moment and source moment fields
  createField<Field_Moments>(data_store, "phi", al_v, *fluxmoment_set);
  createField<Field_Moments>(data_store, "phi_out", al_v, *fluxmoment_set);


  // Create "plane data" to hold face-centered values while sweeping
  Set const &zonei_set = data_store.getVariable<Set>("Set/ZoneI");
  Set const &zonej_set = data_store.getVariable<Set>("Set/ZoneJ");
  Set const &zonek_set = data_store.getVariable<Set>("Set/ZoneK");
  Set const &iplane_set = data_store.newVariable<ProductSet<4>>("Set/IPlane", pspace, SPACE_PQR, dir_set, group_set, zonej_set, zonek_set);
  Set const &jplane_set = data_store.newVariable<ProductSet<4>>("Set/JPlane", pspace, SPACE_PQR, dir_set, group_set, zonei_set, zonek_set);
  Set const &kplane_set = data_store.newVariable<ProductSet<4>>("Set/KPlane", pspace, SPACE_PQR, dir_set, group_set, zonei_set, zonej_set);
  createField<Field_IPlane>(data_store, "i_plane", al_v, iplane_set);
  createField<Field_JPlane>(data_store, "j_plane", al_v, jplane_set);
  createField<Field_KPlane>(data_store, "k_plane", al_v, kplane_set);

  // Create a set to span scattering transfer matrix
  Set const &material_set   = data_store.getVariable<Set>("Set/Material");
  Set const &legendre_set   = data_store.getVariable<Set>("Set/Legendre");
  Set const &global_group_set = data_store.getVariable<Set>("Set/GlobalGroup");
  ProductSet<4> *sigs_set = new ProductSet<4>(pspace, SPACE_NULL,
      material_set, legendre_set, global_group_set, global_group_set);

  data_store.addVariable("Set/SigmaS", sigs_set);


  // Create storage for the scattering transfer matrix
  createField<Field_SigmaS>(data_store, "data/sigs", al_v, *sigs_set);
  auto &field_sigs = data_store.getVariable<Field_SigmaS>("data/sigs");

  // Zero out entire matrix
  Kripke::Kernel::kConst(field_sigs, 0.0);

  // Assign basic diagonal data to matrix
  for(auto sdom_id : field_sigs.getWorkList()){

    // Assign diagonal to the user input for each material
    // Assume each group has same behavior
    auto sigs = field_sigs.getView(sdom_id);
    int global_num_groups = global_group_set.size(sdom_id);
    Legendre n{0};
    for(Material mat{0};mat < 3;++ mat){
      RAJA::forall<RAJA::seq_exec>(
        RAJA::TypedRangeSegment<GlobalGroup>(0, global_num_groups),
        [=](GlobalGroup g){
          sigs(mat, n, g, g) = input_vars.sigs[*mat];
      });
    }
  }



}




