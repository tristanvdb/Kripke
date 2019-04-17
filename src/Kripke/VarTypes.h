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

#ifndef KRIPKE_VARTYPES_H__
#define KRIPKE_VARTYPES_H__

#include <Kripke.h>
#include <Kripke/ArchLayout.h>
#include <Kripke/Core/VarLayout.h>
#include <Kripke/Core/Set.h>
#include <Kripke/Core/Field.h>

namespace Kripke {

  RAJA_INDEX_VALUE(Dimension, "Dimension");
  RAJA_INDEX_VALUE(Direction, "Direction");
  RAJA_INDEX_VALUE(GlobalGroup, "GlobalGroup");
  RAJA_INDEX_VALUE(Group, "Group");
  RAJA_INDEX_VALUE(Legendre, "Legendre");
  RAJA_INDEX_VALUE(Material, "Material");
  RAJA_INDEX_VALUE(MixElem, "MixElem");
  RAJA_INDEX_VALUE(Moment, "Moment");
  RAJA_INDEX_VALUE(Zone, "Zone");
  RAJA_INDEX_VALUE(ZoneI, "ZoneI");
  RAJA_INDEX_VALUE(ZoneJ, "ZoneJ");
  RAJA_INDEX_VALUE(ZoneK, "ZoneK");
  using ZoneX = std::tuple<ZoneI,ZoneJ,ZoneK>;

#if defined(KRIPKE_USE_ZFP)

#if 0
  // There is no case where replacing a whole field by a ZFP array make sense...
  struct double_zfp : public Kripke::Core::field_storage_config {
    using type = double;
    constexpr static bool uses_zfp;
  };
#endif

  // Descriptors for array of zfp array
  //    - `exclude` refers to the number of dimension being excluded from the user defined array
  //         -> we exclude dimension *after* the layout is applied
  //    - `fast_dims` refers to whether we put the faster or upper dimensions form the user defined array
  //         -> `fast_dims == true` => the faster dimension are used
  struct double_zfp_exclude_2_fast : public Kripke::Core::field_storage_config {
    using type = double;
    constexpr static bool uses_zfp = true; // true/false does not matter... FIXME use enum of available user defined array
    constexpr static size_t exclude = 2;
    constexpr static bool fast_dims = true;
  };

  // Exclude: Direction & Group
  using Field_Flux_psi = Kripke::Core::Field<double_zfp_exclude_2_fast, Direction, Group, ZoneI, ZoneJ, ZoneK>; // used in: SweepSubdomain+, LTimes, Population
  using Field_Flux_rhs = Kripke::Core::Field<double_zfp_exclude_2_fast, Direction, Group, ZoneI, ZoneJ, ZoneK>; // used in: LPlusTime+, SweepSubdomain

  // Exclude: Moment & Group
  using Field_Moments_phi     = Kripke::Core::Field<double_zfp_exclude_2_fast, Direction, Group, ZoneI, ZoneJ, ZoneK>; // used in: LTimes+, Scattering
  using Field_Moments_phi_out = Kripke::Core::Field<double_zfp_exclude_2_fast, Direction, Group, ZoneI, ZoneJ, ZoneK>; // used in: Scattering+, Source, LPlusTime
#else
  using Field_Flux_psi = Kripke::Core::Field<double, Direction, Group, ZoneI, ZoneJ, ZoneK>; // used in: SweepSubdomain+, LTimes
  using Field_Flux_rhs = Kripke::Core::Field<double, Direction, Group, ZoneI, ZoneJ, ZoneK>; // used in: LPlusTime+, SweepSubdomain

  using Field_Moments_phi     = Kripke::Core::Field<double, Moment, Group, ZoneI, ZoneJ, ZoneK>; // used in: LTimes+, Scattering
  using Field_Moments_phi_out = Kripke::Core::Field<double, Moment, Group, ZoneI, ZoneJ, ZoneK>; // used in: Scattering+, Source, LPlusTime
#endif

  using Field_IPlane     = Kripke::Core::Field<double, Direction, Group, ZoneJ, ZoneK>;
  using Field_JPlane     = Kripke::Core::Field<double, Direction, Group, ZoneI, ZoneK>;
  using Field_KPlane     = Kripke::Core::Field<double, Direction, Group, ZoneI, ZoneJ>;

  using Field_Ell     = Kripke::Core::Field<double, Moment, Direction>;
  using Field_EllPlus = Kripke::Core::Field<double, Direction, Moment>;

  using Field_Speed  = Kripke::Core::Field<double, Material, GlobalGroup>;
  using Field_SigmaT = Kripke::Core::Field<double, Material, GlobalGroup>;
  using Field_SigmaS = Kripke::Core::Field<double, Material, Legendre, GlobalGroup, GlobalGroup>;

  using Field_Direction2Double = Kripke::Core::Field<double, Direction>; // used three times quadrature/xcos, quadrature/ycos, and quadrature/zcos
  using Field_Direction2Int    = Kripke::Core::Field<int, Direction>;

  using Field_Adjacency        = Kripke::Core::Field<GlobalSdomId, Dimension>;

  using Field_Moment2Legendre  = Kripke::Core::Field<Legendre, Moment>;

  using Field_ZoneI2Double  = Kripke::Core::Field<double, ZoneI>;
  using Field_ZoneJ2Double  = Kripke::Core::Field<double, ZoneJ>;
  using Field_ZoneK2Double  = Kripke::Core::Field<double, ZoneK>;
  using Field_Zone2Double   = Kripke::Core::Field<double, ZoneI, ZoneJ, ZoneK>;
  using Field_Zone2Int      = Kripke::Core::Field<int, ZoneI, ZoneJ, ZoneK>;
  using Field_Zone2MixElem  = Kripke::Core::Field<MixElem, ZoneI, ZoneJ, ZoneK>;

  using Field_MixElem2Double   = Kripke::Core::Field<double, MixElem>;
  using Field_MixElem2Material = Kripke::Core::Field<Material, MixElem>;
  using Field_MixElem2Zone     = Kripke::Core::Field<ZoneX, MixElem>;

  using Field_SigmaTZonal = Kripke::Core::Field<double, Group, ZoneI, ZoneJ, ZoneK>;

  template<typename T>
  struct DefaultOrder{};

  template<typename A, typename L>
  struct DefaultOrder<ArchLayoutT<A, L>> : DefaultOrder<L> {};

  template<>
  struct DefaultOrder<LayoutT_DGZ>{
    using type = camp::list<long, Dimension, Material, Direction, Legendre, Moment, GlobalGroup, Group, Zone, ZoneK, ZoneJ, ZoneI, MixElem>;
  };

  template<>
  struct DefaultOrder<LayoutT_DZG>{
    using type = camp::list<long, Dimension, Material, Direction, Legendre, Moment, Zone, ZoneK, ZoneJ, ZoneI, GlobalGroup, Group, MixElem>;
  };

  template<>
  struct DefaultOrder<LayoutT_GDZ>{
    using type = camp::list<long, Dimension, Material, GlobalGroup, Group, Direction, Legendre, Moment, Zone, ZoneK, ZoneJ, ZoneI, MixElem>;
  };
  
  template<>
  struct DefaultOrder<LayoutT_GZD>{
    using type = camp::list<long, Dimension, Material, GlobalGroup, Group, Zone, ZoneK, ZoneJ, ZoneI, MixElem, Direction, Legendre, Moment>;
  };

  template<>
  struct DefaultOrder<LayoutT_ZDG>{
    using type = camp::list<long, Dimension, Material, Zone, ZoneK, ZoneJ, ZoneI, MixElem, Direction, Legendre, Moment, GlobalGroup, Group>;
  };

  template<>
  struct DefaultOrder<LayoutT_ZGD>{
    using type = camp::list<long, Dimension, Material, Zone, ZoneK, ZoneJ, ZoneI, MixElem, GlobalGroup, Group, Direction, Legendre, Moment>;
  };


  template<typename AL>
  struct SdomAL;

  template<typename A, typename L>
  struct SdomAL<ArchLayoutT<A, L>>
  {
    using al_t = ArchLayoutT<A, L>;
    using arch_t = A;
    using layout_t = L;

    using order_t = typename DefaultOrder<L>::type;

    Kripke::SdomId sdom_id;

    template<typename FieldType>
    auto getView(FieldType &field) const ->
      decltype(field.template getViewOrder<order_t>(sdom_id))
    {
      return field.template getViewOrder<order_t>(sdom_id);
    }
    
    template<typename FieldType>
    auto getView(FieldType &field, Kripke::SdomId sdom) const ->
      decltype(field.template getViewOrder<order_t>(sdom))
    {
      return field.template getViewOrder<order_t>(sdom);
    }
  };

  template<typename AL>
  SdomAL<AL> getSdomAL(AL, Kripke::SdomId sdom_id){ 
    return SdomAL<AL>{sdom_id};
  }


  template<typename FieldType, typename SetType>
  RAJA_INLINE
  FieldType &createField(Core::DataStore &data_store, std::string const &name, ArchLayoutV al_v, SetType const &set)
  {
    FieldType *field = nullptr;
    dispatchLayout(al_v.layout_v, [&](auto layout_t){
      using order_t = typename DefaultOrder<decltype(layout_t)>::type; 
     
      field = new FieldType(set, order_t{});
      data_store.addVariable(name, field);
    });

    return *field;
  };

} // namespace Kripke


#endif
