#pragma once

#include <string>
#include <type_traits>

#include "HyperbolicUpdate_base.h"
#include "mpi/GhostCommunicator_partial_blocks.h"
#ifdef PATCH_LAMBDA
#undef PATCH_LAMBDA
#endif
#ifdef CELL_LAMBDA
#undef CELL_LAMBDA
#endif
#include "foreach_cell/AMRBlockForeachCell_scratch.h"
#include "foreach_cell/ForeachCell_utils.h"

namespace dyablo {
namespace{
using CellIndex     = ForeachCell::CellIndex;
using FieldAccessor = UserData::FieldAccessor;
using offset_t      = typename CellIndex::offset_t;
using PatchForeachCell = AMRBlockForeachCell_scratch;
using PatchArray = PatchForeachCell::CellArray_patch;

inline std::string format_octant_bucket_bound(uint32_t nbOcts)
{
  if (nbOcts % 1024u == 0u)
    return std::to_string(nbOcts / 1024u) + "k";
  return std::to_string(nbOcts);
}

inline std::string make_hydro_patch_bucket_label(uint32_t nbOcts)
{
  if (nbOcts < 1024u)
    return "HydroE::Patch_octs_lt1k";

  uint64_t lower = 1024u;
  uint64_t upper = 2048u;
  while (nbOcts >= upper)
  {
    lower = upper;
    upper *= 2u;
  }

  return "HydroE::Patch_octs_" + format_octant_bucket_bound(static_cast<uint32_t>(lower))
       + "_" + format_octant_bucket_bound(static_cast<uint32_t>(upper));
}

}// namespace
}// namespace dyablo

namespace dyablo {

namespace impl{
namespace {
template <typename T, typename = int, typename=int>
struct HasDensityAndPressure : std::false_type { };

template <typename T>
struct HasDensityAndPressure <T, decltype((void) T::rho, 0), decltype((void) T::p, 0)> : std::true_type { };
}
}

/**
 * @brief Euler update algorithm
 * 
 * @tparam State the type of state to treat
 */
template<typename Policy>
class Hyperbolic_euler : public HyperbolicUpdate {
  static_assert( is_HyperbolicPolicy_v<Policy>,
  "Policy must be wrapped in HyperbolicPolicy_base");

public:
  using PrimState = typename Policy::PrimState;
  using ConsState = typename Policy::ConsState;
  static constexpr int default_hydro_patch_replay_pass_count = 100;

public:
  Hyperbolic_euler(
          ConfigMap& configMap,
          ForeachCell& foreach_cell,
          Timers& timers) 
  : foreach_cell_patch(foreach_cell.get_amr_mesh(), configMap),
    foreach_cell(foreach_cell),
    timers(timers),
    policy_params(Policy::getParams(configMap)),
    ndim(configMap.getValue<int>("mesh", "ndim", 3)),
    smallr( configMap.getValue<real_t>("hydro","smallr", 1e-10) ),
    smallp( configMap.getValue<real_t>("hydro","smallp", 1e-10) ),
    slope_enabled( configMap.getValue<bool>("hydro","slope_enabled", true) ),
    hydro_patch_replay_pass_count( configMap.getValue<int>(
      "hydro",
      "hydro_patch_replay_pass_count",
      default_hydro_patch_replay_pass_count) )
  {
    DYABLO_ASSERT_HOST_RELEASE(
      hydro_patch_replay_pass_count > 0,
      "hydro/hydro_patch_replay_pass_count must be strictly positive, got "
      << hydro_patch_replay_pass_count
    );
  }

  /**
   * @brief Solves hydro for one step using the euler method
   * 
   * @param U the input/output global array
   * @param scalar_data input scalar data
   */
  void update( UserData& U, ScalarSimulationData& scalar_data)
  {
    real_t dt = scalar_data.get<real_t>("dt");
    int ndim = this->ndim;
    bool slope_enabled = this->slope_enabled;

    const Policy policy( this->policy_params, scalar_data ); 
    Timers& timers = this->timers; 
    ForeachCell& foreach_cell = this->foreach_cell;
    PatchForeachCell& foreach_cell_patch = this->foreach_cell_patch;

    FieldAccessor Uin = policy.getUin(U);
    FieldAccessor Uout = policy.getUout(U);
    
    timers.get("HyperbolicUpdate_euler").start();

    ForeachCell::CellMetaData cellmetadata = foreach_cell.getCellMetaData();

    ForeachCell::SearchMode_neighbor search_neighbor( this->foreach_cell.get_amr_mesh().getLightOctree(), ForeachCell::SearchMode_neighbor::ORIGIN );
    ForeachCell::SearchMode_local search_local( ForeachCell::SearchMode_local::ASSERT );

    int nb_ghosts = slope_enabled ? 2 : 1;
    const bool kokkos_tuning_enabled = Kokkos::Tools::Experimental::have_tuning_tool();

    auto initialize_patch_output = [&]()
    {
      foreach_cell.foreach_cell( "HydroE::Init",
        Uout.getShape(),
        KOKKOS_LAMBDA(const CellIndex &iCell)
      {
        ConsState uC = policy.getConsState(Uin, iCell);
        policy.setConsState(Uout, iCell, uC);
      });

      foreach_cell.foreach_ghost_cell( "HydroE::Ghost0",
        Uout.getShape(),
        KOKKOS_LAMBDA(const CellIndex &iCell)
      {
        ConsState empty_state{};
        policy.setConsState(Uout, iCell, empty_state);
      });
    };

    const uint32_t nbOcts = foreach_cell_patch.get_amr_mesh().getNumOctants();
    const std::string hydro_patch_label = make_hydro_patch_bucket_label(nbOcts);

    auto run_hydro_patch = [&]()
    {
      PatchArray::Ref Qpatch_ = foreach_cell_patch.reserve_patch_tmp("Qpatch", nb_ghosts, nb_ghosts, (ndim == 3)?nb_ghosts:0, State_traits<PrimState>::nvars);

      foreach_cell_patch.foreach_patch( hydro_patch_label,
        KOKKOS_LAMBDA( const PatchForeachCell::Patch& patch )
      {
        PatchArray Qpatch = patch.allocate_tmp(Qpatch_);

        patch.foreach_cell( "HydroE::QPatch", Qpatch,
          [&]( const CellIndex& iCell_Qpatch )
        {
          ForeachCell::SearchMode_neighbor search_neighbor_origin( cellmetadata.getLightOctree(), ForeachCell::SearchMode_neighbor::ORIGIN );
          CellIndex iCell_Uin = Uin.getShape().convert_index(iCell_Qpatch, search_neighbor_origin);
          int level_diff = iCell_Uin.level_diff();
          ConsState u = {};
          if (iCell_Uin.is_boundary())
            u = policy.getBoundaryValue(Uin, iCell_Uin, cellmetadata);
          else if (level_diff < 0) {
            int subcell_count =
            foreach_sibling(ndim, iCell_Uin, search_neighbor_origin,
              [&](const CellIndex& iCell_neigh) {
                ConsState uloc = policy.getConsState(Uin, iCell_neigh);
                u += uloc;
              });
            u /= subcell_count;
          }
          else
            u = policy.getConsState(Uin, iCell_Uin);

          const PrimState q = policy.consToPrim( u );
          policy.setPrimState( Qpatch, iCell_Qpatch, q );
        });

        patch.foreach_cell( "HydroE::Flux", Uout.getShape(),
          [&]( const CellIndex& iCell_Uout )
        {
          ForeachCell::SearchMode_neighbor search_neighbor( cellmetadata.getLightOctree(), ForeachCell::SearchMode_neighbor::CLOSEST );

          // Return Slope at position iCell
          auto get_slope = [&](const CellIndex &iCell_Uin, const CellIndex &iCell_Qpatch, ComponentIndex3D dir)
          {
            if(!slope_enabled)
              return PrimState{};

            const PrimState qC = policy.getPrimState(Qpatch, iCell_Qpatch );
            offset_t off_m{}; off_m[dir] = -1;
            const PrimState qL = policy.getPrimState(Qpatch, iCell_Qpatch + off_m);
            offset_t off_p{}; off_p[dir] =  1;
            const PrimState qR = policy.getPrimState(Qpatch, iCell_Qpatch + off_p);

            //!\ Neighbor cells in Qpatch are averaged cells -> size iCell_L != size iCell_Qpatch_L
            CellIndex iCell_L = iCell_Uin.getNeighbor(off_m, search_neighbor);
            CellIndex iCell_R = iCell_Uin.getNeighbor(off_p, search_neighbor);

            // Getting the length right and left
            // Smaller -> use averaged same-size cell -> 1*dx
            // Bigger -> same-size cell in Qpatch has vame value as actual bigger cell -> dx/2 + dx
            constexpr real_t sizes[] = {1.0, 1.0, 1.5};
            const real_t dL = sizes[iCell_L.level_diff()+1];
            const real_t dR = sizes[iCell_R.level_diff()+1];

            // Computing minmod slope for the direction
            PrimState slope = policy.compute_slope( qL, qC, qR, dL, dR);
            return slope;
          }; // get_slope


          auto process_dir = [&](const CellIndex &iCell_Uin, const CellIndex &iCell_Qpatch, ComponentIndex3D dir) {
            // Getting centered value and slope
            PrimState qC0 = policy.getPrimState( Qpatch, iCell_Qpatch );
            PrimState slope_C = get_slope(iCell_Uin, iCell_Qpatch, dir);
            real_t size_C = cellmetadata.getCellSize(iCell_Uin)[dir];

            real_t dim_fac = (ndim == 2 ? 0.5 : 0.25);

            // Compute left side flux
            ConsState fluxL {};
            {
              PrimState qC = qC0 - 0.5 * slope_C;

              offset_t off_m{};
              off_m[dir] = -1;
              const CellIndex iCell_Uin_m = iCell_Uin.getNeighbor(off_m, search_neighbor);
              if( iCell_Uin_m.is_boundary() )
              {
                fluxL = policy.getBoundaryFlux(Uin, iCell_Uin_m, qC, cellmetadata);
              }
              else
              {
                int Ldiff = iCell_Uin_m.level_diff();
                if (Ldiff >= 0)
                {
                  CellIndex iCell_Qpatch_m = iCell_Qpatch + off_m;
                  PrimState qL0 = policy.getPrimState( Qpatch, iCell_Qpatch_m );
                  PrimState slope_L = get_slope(iCell_Uin_m, iCell_Qpatch_m, dir);
                  real_t size_L = cellmetadata.getCellSize(iCell_Uin_m)[dir];

                  // Reconstructing
                  PrimState qL = qL0 + 0.5 * slope_L;

                  // Solving
                  fluxL = policy.riemann_solver(qL, qC, dir);

                  // Adding flux to the neighbor if it is bigger
                  if (Ldiff == 1)
                  {
                    ConsState du_n = fluxL * - dim_fac * dt / size_L;
                    policy.atomic_addConsState(Uout, iCell_Uin_m, du_n);
                  }
                } // If smaller we skip
              }
            }

            // Compute right side flux
            ConsState fluxR {};
            {
              PrimState qC = qC0 + 0.5 * slope_C;

              offset_t off_p{};
              off_p[dir] = 1;
              const CellIndex iCell_Uin_p = iCell_Uin.getNeighbor(off_p, search_neighbor);
              if( iCell_Uin_p.is_boundary() )
              {
                fluxR = policy.getBoundaryFlux(Uin, iCell_Uin_p, qC, cellmetadata);
              }
              else
              {
                int Rdiff = iCell_Uin_p.level_diff();
                if (Rdiff >= 0)
                {
                  CellIndex iCell_Qpatch_p = iCell_Qpatch + off_p;
                  PrimState qR0 = policy.getPrimState( Qpatch, iCell_Qpatch_p );
                  PrimState slope_R = get_slope(iCell_Uin_p, iCell_Qpatch_p, dir);
                  real_t size_R = cellmetadata.getCellSize(iCell_Uin_p)[dir];

                  // Reconstructing
                  PrimState qR = qR0 - 0.5 * slope_R;

                  // Solving
                  fluxR = policy.riemann_solver(qC, qR, dir);

                  // Adding flux to the neighbor if it is bigger
                  if (Rdiff == 1)
                  {
                    ConsState du_n = fluxR * dim_fac * dt / size_R;
                    policy.atomic_addConsState(Uout, iCell_Uin_p, du_n);
                  }
                }
              }
            }

            ConsState du = (fluxL-fluxR) * dt / size_C;
            return du;
          };


          CellIndex iCell_Qpatch = Qpatch.getShape().convert_index(iCell_Uout, search_local);

          ConsState du{};
          du += process_dir(iCell_Uout, iCell_Qpatch, IX);
          du += process_dir(iCell_Uout, iCell_Qpatch, IY);
          if (ndim == 3)
             du += process_dir(iCell_Uout, iCell_Qpatch, IZ);
          policy.atomic_addConsState(Uout, iCell_Uout, du);
        });
      });
    };

    const int pass_count = kokkos_tuning_enabled ? hydro_patch_replay_pass_count : 1;
    for (int iPass = 0; iPass < pass_count; ++iPass)
    {
      initialize_patch_output();
      run_hydro_patch();
    }
      

    // Reducing the ghosts to accumulate the flux in the data arrays 
    int ghost_count = 1;
    GhostCommunicator_partial_blocks ghost_comm ( 
      foreach_cell.get_amr_mesh(),
      Uout.getShape(),
      ghost_count );
    ghost_comm.reduce_ghosts( Uout );
    
    if constexpr ( Policy::has_postProcess() )
    {
      foreach_cell.foreach_cell( "HydroE::Post", Uout.getShape(),
        KOKKOS_LAMBDA(  const ForeachCell::CellIndex& iCell)
      {
        ConsState u = policy.getConsState(Uout, iCell);
        ConsState u_pp = policy.postProcess( u );
        policy.setConsState( Uout, iCell, u_pp );
      });
    }
    
    policy.printWarnings();

    timers.get("HyperbolicUpdate_euler").stop();
  }

private:
  PatchForeachCell foreach_cell_patch;
  ForeachCell& foreach_cell;
  
  Timers& timers;  
  typename Policy::Params policy_params;

  int ndim;
  real_t smallr, smallp;
  bool slope_enabled;
  int hydro_patch_replay_pass_count;
};

} // namespace dyablo
