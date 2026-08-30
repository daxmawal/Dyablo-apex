#include "utils/misc/Dyablo_assert.h"
#include "foreach_cell/ForeachCell.h"
#include "foreach_cell/ForeachCell_utils.h"
#include "hyperbolic/policy/HyperbolicPolicy_Hydro.h"
#include "hyperbolic/scheme/Hyperbolic_euler_nopatch_replay.h"

#include <Kokkos_Core.hpp>
#include <krepe/replayer.hpp>

#include <cstdint>
#include <iostream>

namespace {

using HydroPolicy = dyablo::HyperbolicPolicy_Hydro;
using ForeachCell = dyablo::ForeachCell;
using CellIndex = ForeachCell::CellIndex;
using CellMetaData = ForeachCell::CellMetaData;
using FieldAccessor = dyablo::HyperbolicEulerNoPatchReplay::FieldAccessor;
using ReplayData = dyablo::HyperbolicEulerNoPatchReplay::Data<HydroPolicy>;
using PrimState = HydroPolicy::PrimState;
using ConsState = HydroPolicy::ConsState;
using offset_t = CellIndex::offset_t;
using BlockData = dyablo::AMRBlockForeachCell_CData;
using LightOctree = dyablo::LightOctree;
using dyablo::foreach_smaller_neighbor;

HydroPolicy make_dummy_policy()
{
  using BoundaryConditions =
      dyablo::HyperbolicPolicy_BoundaryConditions_Hydro_dynamic;
  using DefaultBoundaryConditions =
      dyablo::HyperbolicPolicy_BoundaryConditions_Hydro_Default;
  using Slope =
      dyablo::HyperbolicPolicy_Slope_dynamic<dyablo::HyperbolicPolicy_State_Hydro>;
  using MinmodSlope =
      dyablo::HyperbolicPolicy_Slope_minmod<dyablo::HyperbolicPolicy_State_Hydro>;

  const Kokkos::Array<BoundaryConditionType, 3> bc{
      BC_ABSORBING, BC_ABSORBING, BC_ABSORBING};

  HydroPolicy::Params params{
      dyablo::HyperbolicPolicy_Hydro_Params{3, 1.4, 1e-10, 1e-10, 1e-10},
      BoundaryConditions::Params{DefaultBoundaryConditions::Params{bc, bc}},
      Slope::Params{MinmodSlope::Params{}}};

  return HydroPolicy(params, dyablo::ScalarSimulationData{});
}

CellMetaData make_dummy_cell_metadata()
{
  const BlockData block_data{
      3, 1, 1, 1,
      0, 0, 0,
      1, 1, 1,
      1};
  return CellMetaData{block_data, LightOctree{}};
}

} // namespace

int main(int argc, char* argv[])
{
  krepe::ScopeGuard replay_scope(argc, argv);
  Kokkos::ScopeGuard kokkos_scope(argc, argv);

  ReplayData replay_data{
      ForeachCell::IterationSpace_fullArray{ForeachCell::CellArray_shape{}},
      make_dummy_cell_metadata(), true, make_dummy_policy(), FieldAccessor{},
      3, real_t{0}, FieldAccessor{}};

  // Keep the complete replayed computation here for debugging and autotuning.
  auto update_functor = KOKKOS_LAMBDA(uint32_t index)
  {
    const CellIndex iCell = replay_data.get_cell(index);
    const auto& cellmetadata = replay_data.cellmetadata;
    const bool slope_enabled = replay_data.slope_enabled;
    const auto& policy = replay_data.policy;
    const auto& Uin = replay_data.Uin;
    const int ndim = replay_data.ndim;
    const real_t dt = replay_data.dt;
    const auto& Uout = replay_data.Uout;

    ForeachCell::SearchMode_neighbor search_neighbor( cellmetadata.getLightOctree(), ForeachCell::SearchMode_neighbor::CLOSEST );

    // Return Slope at position iCell
    auto get_slope = [&](const CellIndex &iCell, ComponentIndex3D dir) {
      if(!slope_enabled)
          return PrimState{};

      auto get_neighbor_prim_value = [&]( const CellIndex& iCell_n, const CellIndex::offset_t& off )
      {
        ConsState u {};
        // Getting left value
        int level_diff = iCell_n.level_diff();
        if (iCell_n.is_boundary())
          u = policy.getBoundaryValue(Uin, iCell_n, cellmetadata);
        else if (level_diff < 0) {
          int subcell_count =
          foreach_smaller_neighbor(ndim, iCell_n, off, search_neighbor,
            [&](const CellIndex& iCell_neigh) {
              ConsState uloc = policy.getConsState(Uin, iCell_neigh);
              u += uloc;
            });
          u /= subcell_count;
        }
        else
          u = policy.getConsState(Uin, iCell_n);

        return policy.consToPrim(u);
      };

      ConsState uC = policy.getConsState(Uin, iCell);
      const PrimState qC = policy.consToPrim( uC );
      offset_t off_m{}; off_m[dir] = -1;
      CellIndex iCell_L = iCell.getNeighbor(off_m, search_neighbor);
      const PrimState qL = get_neighbor_prim_value(iCell_L, off_m);
      offset_t off_p{}; off_p[dir] =  1;
      CellIndex iCell_R = iCell.getNeighbor(off_p, search_neighbor);
      const PrimState qR = get_neighbor_prim_value(iCell_R, off_p);

      // Getting the length right and left
      constexpr real_t sizes[] = {0.75, 1.0, 1.5};
      const real_t dL = sizes[iCell_L.level_diff()+1];
      const real_t dR = sizes[iCell_R.level_diff()+1];

      // Computing minmod slope for the direction
      PrimState slope = policy.compute_slope( qL, qC, qR, dL, dR);
      return slope;
    }; // get_slope


    auto process_dir = [&](const CellIndex &iCell, ComponentIndex3D dir) {
       // Getting centered value and slope
      ConsState uC = policy.getConsState( Uin, iCell );
      PrimState qC0 = policy.consToPrim(uC);
      PrimState slope_C = get_slope(iCell, dir);
      real_t size_C = cellmetadata.getCellSize(iCell)[dir];

      real_t dim_fac = (ndim == 2 ? 0.5 : 0.25);

      // Compute left side flux
      ConsState fluxL {};
      {
        PrimState qC = qC0 - 0.5 * slope_C;

        offset_t off_m{};
        off_m[dir] = -1;
        const CellIndex iCell_m = iCell.getNeighbor(off_m, search_neighbor);
        if( iCell_m.is_boundary() )
        {
          fluxL = policy.getBoundaryFlux(Uin, iCell_m, qC, cellmetadata);
        }
        else
        {
          int Ldiff = iCell_m.level_diff();
          if (Ldiff >= 0)
          {
            ConsState uL = policy.getConsState( Uin, iCell_m );
            PrimState qL0 = policy.consToPrim(uL);
            PrimState slope_L = get_slope(iCell_m, dir);
            real_t size_L = cellmetadata.getCellSize(iCell_m)[dir];

            // Reconstructing
            PrimState qL = qL0 + 0.5 * slope_L;

            // Solving
            fluxL = policy.riemann_solver(qL, qC, dir);

            // Adding flux to the neighbor if it is bigger
            if (Ldiff == 1)
            {
              ConsState du_n = fluxL * - dim_fac * dt / size_L;
              policy.atomic_addConsState(Uout, iCell_m, du_n);
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
        const CellIndex iCell_p = iCell.getNeighbor(off_p, search_neighbor);
        if( iCell_p.is_boundary() )
        {
          fluxR = policy.getBoundaryFlux(Uin, iCell_p, qC, cellmetadata);
        }
        else
        {
          int Rdiff = iCell_p.level_diff();
          if (Rdiff >= 0)
          {
            ConsState uR = policy.getConsState( Uin, iCell_p );
            PrimState qR0 = policy.consToPrim(uR);
            PrimState slope_R = get_slope(iCell_p, dir);
            real_t size_R = cellmetadata.getCellSize(iCell_p)[dir];

            // Reconstructing
            PrimState qR = qR0 - 0.5 * slope_R;

            // Solving
            fluxR = policy.riemann_solver(qC, qR, dir);

            // Adding flux to the neighbor if it is bigger
            if (Rdiff == 1)
            {
              ConsState du_n = fluxR * dim_fac * dt / size_R;
              policy.atomic_addConsState(Uout, iCell_p, du_n);
            }
          }
        }
      }

      ConsState du = (fluxL-fluxR) * dt / size_C;
      return du;
    };

    ConsState du{};
    du += process_dir(iCell, IX);
    du += process_dir(iCell, IY);
    if (ndim == 3)
      du += process_dir(iCell, IZ);
    policy.atomic_addConsState(Uout, iCell, du);
  };

  krepe::parallel_for(
      "Hyperbolic_euler::update",
      Kokkos::RangePolicy<>(0, replay_data.cell_count()),
      update_functor);
  Kokkos::fence();

  std::cout << "replay kernel completed\n";
  return 0;
}
