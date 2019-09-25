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

#ifndef KRIPKE_CORE_FIELD_H__
#define KRIPKE_CORE_FIELD_H__

#include <Kripke.h>
#include <Kripke/Core/VarLayout.h>
#include <Kripke/Core/DataStore.h>
#include <Kripke/Core/DomainVar.h>
#include <Kripke/Core/Set.h>
#include <vector>

#ifdef KRIPKE_USE_CHAI
#define DEBUG
#include <chai/ManagedArray.hpp>
#undef DEBUG
#endif

#if defined(KRIPKE_USE_ZFP)
#  define KRIPKE_DECOUPLE_STORAGE_FROM_ELEMENTS
#  include "zfparray1.h"
#  include "zfparray2.h"
#  include "zfparray3.h"
#  include "zfpvarray2.h"
#  include "zfpvarray3.h"
#endif

namespace Kripke {

namespace Config {

namespace Field {

  enum class StorageEngine {
    raw,
#if defined(KRIPKE_USE_ZFP)
    zfp,
    arc_a,
    arc_p
#endif
  };

  struct StaticBase {};

  template <
    typename ElementT,
    StorageEngine engine__ = StorageEngine::raw,
    size_t exclude__ = 0,
    bool fast_dims__ = true
  >
  struct Static : public StaticBase {
    using type = double;
    constexpr static StorageEngine engine = engine__;
    constexpr static size_t exclude = exclude__;
    constexpr static bool fast_dims = fast_dims__;
  };

  struct Runtime {
#ifdef KRIPKE_USE_ZFP
    size_t block_size{512};
    size_t cached_blocks{1024};
    double compression_rate{0.};
    double precision{0.};
    size_t accuracy{0};
#endif
  };

  extern std::map<std::string, ::Kripke::Config::Field::Runtime> field_config_map;
}
}

namespace Core {

#if defined(KRIPKE_DECOUPLE_STORAGE_FROM_ELEMENTS)

  namespace StorageHelper {

    namespace ViewConfig {
      template <typename Atype>
      struct array : public RAJA::view_config_user_defined_array {
        using array_type = Atype;
      };

      template <typename Atype, size_t N, size_t exclude, bool fast_dims>
      struct array_of_array;

      template <typename Atype, size_t N, size_t exclude>
      struct array_of_array<Atype, N, exclude, true> : public RAJA::view_config_array_of_user_defined_array {
        using array_type = Atype;
        constexpr static bool user_array_on_fast_dims = true;
        constexpr static size_t num_slow_dims = exclude;
      };

      template <typename Atype, size_t N, size_t exclude>
      struct array_of_array<Atype, N, exclude, false> : public RAJA::view_config_array_of_user_defined_array {
        using array_type = Atype;
        constexpr static bool user_array_on_fast_dims = false;
        constexpr static size_t num_slow_dims = N - exclude;
      };
    }

    namespace ArraySelector {

      template <typename ELEMENT, size_t N>
      struct ZFP;

      template <typename ELEMENT>
      struct ZFP<ELEMENT, 1>{
        using type = ::zfp::array1<ELEMENT>;

        static inline type * alloc(size_t size, ::Kripke::Config::Field::Runtime const & fldcfg) {
          KRIPKE_ASSERT(fldcfg.compression_rate > 0., "Compression rate must be set...");
          return new type(size, fldcfg.compression_rate, NULL, fldcfg.cached_blocks * fldcfg.block_size);
        }

        template <typename LayoutT>
        static inline void layout(type * & store, typename type::pointer & baseptr, LayoutT const & layout) {
          store->resize(layout.sizes[0]);
          baseptr = &((*store)[0]);
        }
      };

      template <typename ELEMENT>
      struct ZFP<ELEMENT, 2>{
        using type = ::zfp::array2<ELEMENT>;

        static inline type * alloc(size_t size, ::Kripke::Config::Field::Runtime const & fldcfg) {
          KRIPKE_ASSERT(fldcfg.compression_rate > 0., "Compression rate must be set...");
          return new type(size, 1, fldcfg.compression_rate, NULL, fldcfg.cached_blocks * fldcfg.block_size);
        }

        template <typename LayoutT>
        static inline void layout(type * & store, typename type::pointer & baseptr, LayoutT const & layout) {
          store->resize(layout.sizes[0], layout.sizes[1]);
        }
      };

      template <typename ELEMENT>
      struct ZFP<ELEMENT, 3>{
        using type = ::zfp::array3<ELEMENT>;

        static inline type * alloc(size_t size, ::Kripke::Config::Field::Runtime const & fldcfg) {
          KRIPKE_ASSERT(fldcfg.compression_rate > 0., "Compression rate must be set...");
          return new type(size, 1, 1, fldcfg.compression_rate, NULL, fldcfg.cached_blocks * fldcfg.block_size);
        }

        template <typename LayoutT>
        static inline void layout(type * & store, typename type::pointer & baseptr, LayoutT const & layout) {
          store->resize(layout.sizes[0], layout.sizes[1], layout.sizes[2]);
          baseptr = &((*store)[0]);
        }
      };

      template <typename ELEMENT, size_t N>
      struct ARC_P;

      template <typename ELEMENT>
      struct ARC_P<ELEMENT, 2>{
        using type = ::zfp::varray2<ELEMENT>;

        static inline type * alloc(size_t size, ::Kripke::Config::Field::Runtime const & fldcfg) {
          KRIPKE_ASSERT(fldcfg.precision > 0., "Precision must be set...");
          return new type(size, 1, fldcfg.precision, NULL, fldcfg.cached_blocks * fldcfg.block_size);
        }

        template <typename LayoutT>
        static inline void layout(type * & store, typename type::pointer & baseptr, LayoutT const & layout) {
          store->resize(layout.sizes[0], layout.sizes[1]);
          baseptr = &((*store)[0]);
        }
      };

      template <typename ELEMENT>
      struct ARC_P<ELEMENT, 3>{
        using type = ::zfp::varray3<ELEMENT>;

        static inline type * alloc(size_t size, ::Kripke::Config::Field::Runtime const & fldcfg) {
          KRIPKE_ASSERT(fldcfg.precision > 0., "Precision must be set...");
          return new type(size, 1, 1, fldcfg.precision, NULL, fldcfg.cached_blocks * fldcfg.block_size);
        }

        template <typename LayoutT>
        static inline void layout(type * & store, typename type::pointer & baseptr, LayoutT const & layout) {
          store->resize(layout.sizes[0], layout.sizes[1], layout.sizes[2]);
          baseptr = &((*store)[0]);
        }
      };

      template <typename ELEMENT, size_t N>
      struct ARC_A;

      template <typename ELEMENT>
      struct ARC_A<ELEMENT, 2>{
        using type = ::zfp::varray2<ELEMENT>;

        static inline type * alloc(size_t size, ::Kripke::Config::Field::Runtime const & fldcfg) {
          KRIPKE_ASSERT(fldcfg.accuracy > 0., "Accuracy must be set...");
          return new type(size, 1, (uint)fldcfg.accuracy, NULL, fldcfg.cached_blocks * fldcfg.block_size);
        }

        template <typename LayoutT>
        static inline void layout(type * & store, typename type::pointer & baseptr, LayoutT const & layout) {
          store->resize(layout.sizes[0], layout.sizes[1]);
          baseptr = &((*store)[0]);
        }
      };

      template <typename ELEMENT>
      struct ARC_A<ELEMENT, 3>{
        using type = ::zfp::varray3<ELEMENT>;

        static inline type * alloc(size_t size, ::Kripke::Config::Field::Runtime const & fldcfg) {
          KRIPKE_ASSERT(fldcfg.accuracy > 0., "Accuracy must be set...");
          return new type(size, 1, 1, (uint)fldcfg.accuracy, NULL, fldcfg.cached_blocks * fldcfg.block_size);
        }

        template <typename LayoutT>
        static inline void layout(type * & store, typename type::pointer & baseptr, LayoutT const & layout) {
          store->resize(layout.sizes[0], layout.sizes[1], layout.sizes[2]);
          baseptr = &((*store)[0]);
        }
      };
      
    }

    namespace utils {
      template<typename T, typename = void>
      struct has_exclude : std::false_type { };

      template<typename T>
      struct has_exclude<T, decltype(std::declval<T>().exclude, void())> : std::true_type { };

      // Provide interface compatible with Kripke's layout (typename LayoutT)
      struct split_layout_t {
        std::vector<size_t> sizes;
      };
    }

    // ElementType:
    //     Used to determine the `underlying_type` when FieldStorage<ConfigT> is defined to inherit from FieldStorageBase<underlying_type>
    //     Needed as Select<ConfigT> is referenced later and introduces additional complexity

    template <
      typename ConfigT,
      bool = std::is_base_of<::Kripke::Config::Field::StaticBase, ConfigT>::value
    >
    struct ElementType;

    template <typename StorageT>
    struct ElementType<StorageT, false> {
      using type = StorageT;
    };

    template <typename ConfigT>
    struct ElementType<ConfigT, true> {
      using type = typename ConfigT::type;
    };

    // RAW:
    //    This corresponds to the traditional contiguous storage found in C/C++
    //    Defined aside because two paths: ConfigT is native type or it uses the `raw` engine

    template <typename CONFIG, typename ELEMENT>
    struct RAW {
      using config_type = CONFIG;
      using element_type = ELEMENT;

      using storage_type = void; // only used with ZFP
      using storage_pointer = element_type *;
      using element_pointer = element_type *;
      using element_reference = element_type &;

      using view_pointer_type = element_type *;

      static inline storage_pointer alloc(size_t size, ::Kripke::Config::Field::Runtime const &) {
        return new element_type[size];
      }

      static inline void free(storage_pointer & store, element_pointer &) {
        delete [] store;
      }

      template <typename LayoutT>
      static inline void layout(storage_pointer & store, element_pointer & baseptr, LayoutT const &, ::Kripke::Config::Field::Runtime const &) {
        baseptr = store;
      }
    };

    // UserDefined:
    //     Generic helper for user-defined storage engine

    template < typename ConfigT, size_t N, ::Kripke::Config::Field::StorageEngine = ConfigT::engine >
    struct UserDefined;

    template < typename ConfigT, size_t N >
    struct UserDefined<ConfigT, N, ::Kripke::Config::Field::StorageEngine::raw> : public RAW<ConfigT, typename ConfigT::type> {};

    // Select:
    //     Entry-point of the StorageHelper. Distinguish between traditional and encapsulated configuration

    template < typename ConfigT, size_t N, bool = std::is_base_of<::Kripke::Config::Field::StaticBase, ConfigT>::value >
    struct Select;

    template <typename StorageT, size_t N>
    struct Select<StorageT, N, false> : public RAW<StorageT, StorageT> {};

    template <typename ConfigT, size_t N>
    struct Select<ConfigT, N, true> : public UserDefined<ConfigT, N> {};

    template <
      typename ConfigT,
      size_t N,
      template <typename, size_t> class ArraySelectorT,
      bool = (ConfigT::exclude > 0)
    >
    struct UserArray;

    template <typename ConfigT, size_t N, template <typename, size_t> class ArraySelectorT>
    struct UserArray<ConfigT, N, ArraySelectorT, false> {
      using config_type = ConfigT;

      using element_type = typename config_type::type;

      using array_selector = ArraySelectorT<element_type, N>;

      using storage_type = typename array_selector::type;
      using storage_pointer = storage_type *;
      using element_pointer = typename storage_type::pointer;
      using element_reference = typename storage_type::reference;

      using view_pointer_type = ViewConfig::array<storage_type>;

      static inline storage_pointer alloc(size_t size, ::Kripke::Config::Field::Runtime const & fldcfg) {
        return array_selector::alloc(size, fldcfg);
      }

      static inline void free(storage_pointer & store, element_pointer &) {
        delete store;
      }

      template <typename LayoutT>
      static inline void layout(storage_pointer & store, element_pointer & baseptr, LayoutT const & layout, ::Kripke::Config::Field::Runtime const &) {
        static_assert(LayoutT::n_dims == N, "Sizes of the layout and storage do not match.");
        array_selector::layout(store, baseptr, layout);
      }
    };

    template <typename ConfigT, size_t N, template <typename, size_t> class ArraySelectorT>
    struct UserArray<ConfigT, N, ArraySelectorT, true> {

      using config_type = ConfigT;

      using element_type = typename config_type::type;

      constexpr static size_t user_dims = N - config_type::exclude;
      static_assert( user_dims > 1 && user_dims < 4 , "Invalid number of ARC dimensions requested for a mixed ARC array.");

      using array_selector = ArraySelectorT<element_type, user_dims>;

      using storage_type = typename array_selector::type;
      using storage_pointer = storage_type **;
      using element_pointer = typename storage_type::pointer *;
      using element_reference = typename storage_type::reference;

      using view_pointer_type = ViewConfig::array_of_array<storage_type, N, config_type::exclude, config_type::fast_dims>;

      static inline storage_pointer alloc(size_t, ::Kripke::Config::Field::Runtime const &) {
        return nullptr;
      }

      static inline void free(storage_pointer & store, element_pointer & baseptr) {
        delete [] store;   // FIXME leaking ARC arrays
        delete [] baseptr; // FIXME leaking ARC pointers
        // TODO need to know the size of the array of array
      }

      static inline element_pointer get_pointer(storage_pointer store) {
        return store;
      }

      template <typename LayoutT>
      static inline void layout(storage_pointer & store, element_pointer & baseptr, LayoutT const & layout, ::Kripke::Config::Field::Runtime const & fldcfg) {
        static_assert(LayoutT::n_dims == N, "Sizes of the layout and storage do not match.");

        utils::split_layout_t slow_layout, fast_layout;

        // Split between slow and fast dimensions and compute the product of each split.
        size_t slow_size = 1;
        size_t fast_size = 1;
        for (size_t i = 0; i < N; i++) {
          if (i < view_pointer_type::num_slow_dims) {
            slow_size *= layout.sizes[i];
            slow_layout.sizes.push_back(layout.sizes[i]);
          } else {
            fast_size *= layout.sizes[i];
            fast_layout.sizes.push_back(layout.sizes[i]);
          }
        }

        if (config_type::fast_dims) {
          store = new storage_type*[slow_size];
          baseptr = new typename storage_type::pointer[slow_size];
          for (size_t i = 0; i < slow_size; i++) {
            store[i] = array_selector::alloc(fast_size, fldcfg);
            array_selector::layout(store[i], baseptr[i], fast_layout);
          }
        } else {
          store = new storage_type*[fast_size];
          baseptr = new typename storage_type::pointer[fast_size];
          for (size_t i = 0; i < fast_size; i++) {
            store[i] = array_selector::alloc(slow_size, fldcfg);
            array_selector::layout(store[i], baseptr[i], slow_layout);
          }
        }
      }
    };

#if defined(KRIPKE_USE_ZFP)
    template <typename ConfigT, size_t N>
    struct UserDefined<ConfigT, N, ::Kripke::Config::Field::StorageEngine::zfp> : public UserArray<ConfigT, N, ArraySelector::ZFP> {};

    template <typename ConfigT, size_t N>
    struct UserDefined<ConfigT, N, ::Kripke::Config::Field::StorageEngine::arc_p> : public UserArray<ConfigT, N, ArraySelector::ARC_P> {};

    template <typename ConfigT, size_t N>
    struct UserDefined<ConfigT, N, ::Kripke::Config::Field::StorageEngine::arc_a> : public UserArray<ConfigT, N, ArraySelector::ARC_A> {};
#endif
  }
#endif

  template<typename ELEMENT>
  class FieldStorageBase : public Kripke::Core::DomainVar {
    public:
      explicit FieldStorageBase(Kripke::Core::Set const &spanned_set) :
        m_set(&spanned_set)
      {}

      virtual ~FieldStorageBase() {}

      RAJA_INLINE
      Kripke::Core::Set const &getSet() const {
        return *m_set;
      }

      /**
       * Returns the number of elements in this subdomain.
       */
      RAJA_INLINE
      size_t size(Kripke::SdomId sdom_id) const {
        size_t chunk_id = Kripke::Core::DomainVar::m_subdomain_to_chunk[*sdom_id];
        return m_chunk_to_size[chunk_id];
      }

      /**
       * Returns the total number of elements in all subdomains.
       */
      RAJA_INLINE
      size_t size() const {
        size_t size = 0;
        for(size_t chunk_id = 0; chunk_id < m_chunk_to_size.size(); ++chunk_id)
          size += m_chunk_to_size[chunk_id];
        return size;
      }

      /**
       * Returns the storage size in this subdomain.
       */
      RAJA_INLINE
      size_t storage_size(Kripke::SdomId sdom_id) const {
//        size_t chunk_id = Kripke::Core::DomainVar::m_subdomain_to_chunk[*sdom_id];
//        return m_chunk_to_storage_size[chunk_id].first;
        return 0;
      }

      /**
       * Returns the total storage size in all subdomains.
       */
      RAJA_INLINE
      size_t storage_size() const {
//        size_t r = 0;
//        for(size_t chunk_id = 0; chunk_id < m_chunk_to_storage_size.size(); ++chunk_id)
//          r += m_chunk_to_storage_size[chunk_id].first;
//        return r;
        return 0;
      }

      /**
       * Returns the storage size in this subdomain.
       */
      RAJA_INLINE
      size_t cache_size(Kripke::SdomId sdom_id) const {
//        size_t chunk_id = Kripke::Core::DomainVar::m_subdomain_to_chunk[*sdom_id];
//        return m_chunk_to_storage_size[chunk_id].second;
        return 0;
      }

      /**
       * Returns the total storage size in all subdomains.
       */
      RAJA_INLINE
      size_t cache_size() const {
//        size_t r = 0;
//        for(size_t chunk_id = 0; chunk_id < m_chunk_to_storage_size.size(); ++chunk_id)
//          r += m_chunk_to_storage_size[chunk_id].second;
//        return r;
        return 0;
      }

    protected:
      Kripke::Core::Set const * m_set;
      std::vector<size_t> m_chunk_to_size;
//    std::vector< std::pair<size_t, size_t> > m_chunk_to_storage_size;
  };

  /**
   * Base class for Field which provides storage allocation
   */
#if defined(KRIPKE_DECOUPLE_STORAGE_FROM_ELEMENTS)
  template<typename ConfigT, unsigned int N>
  class FieldStorage : public FieldStorageBase< typename StorageHelper::ElementType<ConfigT>::type > {
#else
  template<typename ELEMENT>
  class FieldStorage : public FieldStorageBase< ELEMENT > {
#endif
    public:

#if defined(KRIPKE_DECOUPLE_STORAGE_FROM_ELEMENTS)
      using StorageHelperType = StorageHelper::Select<ConfigT, N>;
      using ConfigType = typename StorageHelperType::config_type;
      using StoragePtr = typename StorageHelperType::storage_pointer;
      using ElementType = typename StorageHelperType::element_type;
      using ElementPtr = typename StorageHelperType::element_pointer;
      using ElementRef = typename StorageHelperType::element_reference;
      using ViewPtr = typename StorageHelperType::view_pointer_type;
#else
      using ElementType = ELEMENT;
      using StoragePtr = ElementType *;
      using ElementRef = ElementType;
#if defined(KRIPKE_USE_CHAI)
      using ElementPtr = chai::ManagedArray<ElementType>;
#else
      using ElementPtr = ElementType *;
#endif
      using ViewPtr = ElementPtr;
#endif

      using Layout1dType = RAJA::TypedLayout<RAJA::Index_type, camp::tuple<RAJA::Index_type>>;
      using View1dType = RAJA::View<ElementType, Layout1dType, ViewPtr>;

      using Base = Kripke::Core::FieldStorageBase< ElementType >;

      explicit FieldStorage(Kripke::Core::Set const &spanned_set, ::Kripke::Config::Field::Runtime const & fldcfg) :
        Base(spanned_set)
      {

        // initialize our decomposition to match that of the specified set
        Kripke::Core::DomainVar::setup_initChunks(spanned_set);

        // allocate all of our chunks, and create layouts for each one
        size_t num_chunks = Kripke::Core::DomainVar::m_chunk_to_subdomain.size();
        Base::m_chunk_to_size.resize(num_chunks, 0);
#if defined(KRIPKE_USE_CHAI)
        m_chunk_to_data.resize(num_chunks);
#else
        m_chunk_to_data.resize(num_chunks, nullptr);
#if defined(KRIPKE_DECOUPLE_STORAGE_FROM_ELEMENTS)
        m_chunk_to_base_ptr.resize(num_chunks, nullptr);
#endif
#endif

        for(size_t chunk_id = 0;chunk_id < num_chunks;++ chunk_id){

          // Get the size of the subdomain from the set
          SdomId sdom_id(Kripke::Core::DomainVar::m_chunk_to_subdomain[chunk_id]);
          size_t sdom_size = spanned_set.size(sdom_id);

          Base::m_chunk_to_size[chunk_id] = sdom_size;
#if defined(KRIPKE_USE_CHAI)
          m_chunk_to_data[chunk_id].allocate(sdom_size, chai::CPU,
              [=](chai::Action action, chai::ExecutionSpace space, size_t bytes){
                /*printf("CHAI[%s, %d]: ", BaseVar::getName().c_str(), (int)chunk_id);
                switch(action){
                case chai::ACTION_ALLOC: printf("ALLOC "); break;
                case chai::ACTION_FREE: printf("FREE  "); break;
                case chai::ACTION_MOVE: printf("MOVE  "); break;
                default: printf("UNKNOWN ");
                }

                switch(space){
                case chai::CPU: printf("CPU "); break;
#ifdef KRIPKE_USE_CUDA
                case chai::GPU: printf("GPU  "); break;
#endif
                default: printf("UNK ");
                }

                printf("%lu bytes\n", (unsigned long) bytes);
*/
              }

          );
#elif defined(KRIPKE_DECOUPLE_STORAGE_FROM_ELEMENTS)
          m_chunk_to_data[chunk_id] = StorageHelperType::alloc(sdom_size, fldcfg);
#else
          m_chunk_to_data[chunk_id] = new ElementType[sdom_size];
#endif
        }
      }

      virtual ~FieldStorage(){
#if defined(KRIPKE_USE_CHAI)
#elif defined(KRIPKE_DECOUPLE_STORAGE_FROM_ELEMENTS)
        for(size_t i = 0; i < m_chunk_to_data.size(); i++) {
          StorageHelperType::free(m_chunk_to_data[i], m_chunk_to_base_ptr[i]);
        }
#else
        for(auto i : m_chunk_to_data){
          delete[] i;
        }
#endif
      }

      // Dissallow copy construction
#if defined(KRIPKE_DECOUPLE_STORAGE_FROM_ELEMENTS)
      FieldStorage(FieldStorage<ConfigType, N> const &) = delete;
#else
      FieldStorage(FieldStorage<ElementType> const &) = delete;
#endif


      RAJA_INLINE
      View1dType getView1d(Kripke::SdomId sdom_id) const {
        size_t chunk_id = Kripke::Core::DomainVar::m_subdomain_to_chunk[*sdom_id];

#if defined(KRIPKE_DECOUPLE_STORAGE_FROM_ELEMENTS)
        ElementPtr ptr = m_chunk_to_base_ptr[chunk_id];
#else
        ElementPtr ptr = m_chunk_to_data[chunk_id];
#endif
        size_t sdom_size = Base::m_chunk_to_size[chunk_id];

        return View1dType(ptr, Layout1dType(sdom_size));
      }

      RAJA_INLINE
      StoragePtr getData(Kripke::SdomId sdom_id) const {

        KRIPKE_ASSERT(*sdom_id < (int)Kripke::Core::DomainVar::m_subdomain_to_chunk.size(),
            "sdom_id(%d) >= num_subdomains(%d)",
            (int)*sdom_id,
            (int)(int)Kripke::Core::DomainVar::m_subdomain_to_chunk.size());
        size_t chunk_id = Kripke::Core::DomainVar::m_subdomain_to_chunk[*sdom_id];


#if defined(KRIPKE_USE_CHAI)
        // use pointer conversion to get host pointer
        ElementType *ptr = m_chunk_to_data[chunk_id];
        // return host pointer
        return(ptr);
#else
        return  m_chunk_to_data[chunk_id];
#endif
      }

#if 0
      RAJA_INLINE
      void setRate() const {
        // TODO change current rate/acc/prec
      }
#endif

    protected:
#if defined(KRIPKE_DECOUPLE_STORAGE_FROM_ELEMENTS)
      std::vector<StoragePtr> m_chunk_to_data;
      std::vector<ElementPtr> m_chunk_to_base_ptr;
#else
      std::vector<ElementPtr> m_chunk_to_data;
#endif
  };

  /**
   * Defines a multi-dimensional data field defined over a Set
   */
  template<typename ELEMENT, typename ... IDX_TYPES>
#if defined(KRIPKE_DECOUPLE_STORAGE_FROM_ELEMENTS)
  class Field : public Kripke::Core::FieldStorage<ELEMENT, sizeof...(IDX_TYPES)> {
#else
  class Field : public Kripke::Core::FieldStorage<ELEMENT> {
#endif
    public:

#if defined(KRIPKE_DECOUPLE_STORAGE_FROM_ELEMENTS)
      using Base = Kripke::Core::FieldStorage<ELEMENT, sizeof...(IDX_TYPES)>;
#else
      using Base = Kripke::Core::FieldStorage<ELEMENT>;
#endif

      using ElementType = typename Base::ElementType;
      using ElementPtr = typename Base::ElementPtr;
      using ElementRef = typename Base::ElementRef;
      using ViewPtr = typename Base::ViewPtr;
#if defined(KRIPKE_DECOUPLE_STORAGE_FROM_ELEMENTS)
      using StorageHelperType = typename Base::StorageHelperType;
      using ConfigType = typename Base::ConfigType;
      using StoragePtr = typename Base::StoragePtr;
#endif

      static constexpr size_t NumDims = sizeof...(IDX_TYPES);

      using DefaultLayoutType = RAJA::TypedLayout<RAJA::Index_type, camp::tuple<IDX_TYPES...>>;

      using DefaultViewType = RAJA::View<ElementType, DefaultLayoutType, ViewPtr>;

      template<typename Order>
      Field(Kripke::Core::Set const &spanned_set, Order, ::Kripke::Config::Field::Runtime const & fldcfg = ::Kripke::Config::Field::Runtime()) :
        Base(spanned_set, fldcfg)
      {

        KRIPKE_ASSERT(NumDims == spanned_set.getNumDimensions(),
            "Number of dimensions must match between Field<%d> and Set<%d>\n",
            (int)NumDims, (int)spanned_set.getNumDimensions());

        auto perm = LayoutInfo<Order, IDX_TYPES...>::getPermutation();

        // create layouts for each chunk
        size_t num_chunks = Kripke::Core::DomainVar::m_chunk_to_subdomain.size();
        m_chunk_to_layout.resize(num_chunks);
        for(size_t chunk_id = 0;chunk_id < num_chunks;++ chunk_id){

          // Create a layout using dim sizes from the Set, and permutation
          // defined by the layout function
          SdomId sdom_id(Kripke::Core::DomainVar::m_chunk_to_subdomain[chunk_id]);
          std::array<RAJA::Index_type, NumDims> sizes;
          for(size_t dim = 0;dim < NumDims;++ dim){
            sizes[dim] = spanned_set.dimSize(sdom_id, dim);
          }

          RAJA::Layout<NumDims, RAJA::Index_type> &layout =
              m_chunk_to_layout[chunk_id];
          layout = RAJA::make_permuted_layout<NumDims,RAJA::Index_type>(sizes, perm);

#if defined(KRIPKE_DECOUPLE_STORAGE_FROM_ELEMENTS)
          StorageHelperType::layout(
              Base::m_chunk_to_data[chunk_id], Base::m_chunk_to_base_ptr[chunk_id], layout, fldcfg
          );
#endif
        }
      }

      virtual ~Field(){

      }



      RAJA_INLINE
      DefaultViewType getView(Kripke::SdomId sdom_id) const {

        size_t chunk_id = Kripke::Core::DomainVar::m_subdomain_to_chunk[*sdom_id];

#if defined(KRIPKE_DECOUPLE_STORAGE_FROM_ELEMENTS)
        ElementPtr ptr = Base::m_chunk_to_base_ptr[chunk_id];
#else
        ElementPtr ptr = Base::m_chunk_to_data[chunk_id];
#endif
        auto layout = m_chunk_to_layout[chunk_id];

        return DefaultViewType(ptr, layout);
      }


      template<typename Order>
      RAJA_INLINE
      auto getViewOrder(Kripke::SdomId sdom_id) const ->
        ViewType<Order, ElementType, ViewPtr, IDX_TYPES...>
      {
        size_t chunk_id = Kripke::Core::DomainVar::m_subdomain_to_chunk[*sdom_id];

#if defined(KRIPKE_DECOUPLE_STORAGE_FROM_ELEMENTS)
        ElementPtr ptr = Base::m_chunk_to_base_ptr[chunk_id];
#else
        ElementPtr ptr = Base::m_chunk_to_data[chunk_id];
#endif

        using LInfo = LayoutInfo<Order, IDX_TYPES...>;
        using LType = typename LInfo::Layout;

        LType layout = RAJA::make_stride_one<LInfo::stride_one_dim>(m_chunk_to_layout[chunk_id]);

        return ViewType<Order, ElementType, ViewPtr, IDX_TYPES...>(ptr, layout);
      }



      RAJA_INLINE
      void dump() const {
        printf("Field<>:\n");
        printf("  name:  %s\n", BaseVar::getName().c_str());
        printf("  m_set: %p\n", Base::Base::m_set);

        printf("  m_chunk_to_size: ");
        for(auto x : Base::Base::m_chunk_to_size){printf("%lu ", (unsigned long)x);}
        printf("\n");

#if !defined(KRIPKE_USE_CHAI)
        printf("  m_chunk_to_data: ");
        for(auto x : Base::m_chunk_to_data){printf("%p ", x);}
        printf("\n");
#endif

        for(size_t chunk_id = 0;chunk_id < Base::m_chunk_to_data.size();++ chunk_id){

          SdomId sdom_id(Kripke::Core::DomainVar::m_chunk_to_subdomain[chunk_id]);

          ElementType *ptr = Base::getData(sdom_id);

          printf("Chunk %d Data: ", (int)chunk_id);
          for(size_t i = 0;i < Base::Base::m_chunk_to_size[chunk_id];++ i){
            printf(" %e", ptr[i]);
          }
          printf("\n");
        }

        Kripke::Core::DomainVar::dump();
      }

    protected:
      std::vector<DefaultLayoutType> m_chunk_to_layout;
  };

} } // namespace

#endif

