#include "DyabloSession.hpp"
#include "amr/AMRmesh.h"
#include "utils/misc/Dyablo_assert.h"
#include "foreach_cell/ForeachCell.h"
#include "foreach_cell/ForeachCell_utils.h"
#include "hyperbolic/policy/HyperbolicPolicy_Hydro.h"
#include "user_data/FieldAccessor.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_StdAlgorithms.hpp>
#include <krepe/replayer.hpp>

#include <cstdint>
#include <iostream>
#include <stdexcept>

namespace {

using HydroPolicy = dyablo::HyperbolicPolicy_Hydro;
using ForeachCell = dyablo::ForeachCell;
using CellIndex = ForeachCell::CellIndex;
using CellMetaData = ForeachCell::CellMetaData;
using FieldAccessor = dyablo::UserData::FieldAccessor;
using PrimState = HydroPolicy::PrimState;
using ConsState = HydroPolicy::ConsState;
using offset_t = CellIndex::offset_t;
using dyablo::foreach_smaller_neighbor;

} // namespace

int main(int argc, char* argv[])
try
{
  krepe::ScopeGuard replay_scope(argc, argv);
  dyablo::DyabloSession dyablo_session(argc, argv);

  // Construct valid placeholders, KREPE restores their captured values.
  dyablo::AMRmesh mesh{3, {false, false, false}, 0, 0};
  CellMetaData cellmetadata({}, mesh);
  bool slope_enabled = true;
  const HydroPolicy policy({}, {});
  FieldAccessor Uin{};
  int ndim = 0;
  real_t dt = 0;
  // core/src/foreach_cell/AMRBlockForeachCell_CellArray.h #641
  FieldAccessor Uout{};

  const auto f = KOKKOS_LAMBDA(const CellIndex& iCell)
  {
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

  // Placeholders restored by KREPE, together with the execution policy.
  const ForeachCell::IterationSpace_fullArray iter_space(ForeachCell::CellArray_shape{});
  uint32_t bx = 0;
  uint32_t by = 0;
  uint32_t nbCellsPerBlock = 0;

  // This is the outer lambda captured in AMRBlockForeachCell::foreach_cell.
  krepe::parallel_for("Hyperbolic_euler::update", 0, KOKKOS_LAMBDA(uint32_t index)
  {
    uint32_t iOct = index/nbCellsPerBlock;
    index = index%nbCellsPerBlock;

    uint32_t k = index/(bx*by);
    uint32_t j = (index - k*bx*by)/bx;
    uint32_t i = index - j*bx - k*bx*by;

    CellIndex iCell = iter_space.getCellIndex(iOct, i, j, k);
    f(iCell);
  });

  Kokkos::fence();
  using MemorySpace = Kokkos::DefaultExecutionSpace::memory_space;
  const auto same_value = KOKKOS_LAMBDA(real_t expected, real_t actual) {
    return actual == expected;
  };
  for (const auto& allocation : krepe::get_allocations<MemorySpace>("UserData_fields"))
  {
    krepe::compare_views<real_t*, MemorySpace>(allocation,
      [&allocation, &same_value](const auto& expected, const auto& actual) {
        if (!Kokkos::Experimental::equal(Kokkos::DefaultExecutionSpace{}, expected, actual, same_value))
          throw std::runtime_error("Output mismatch in " + allocation.label);
      });
  }
  std::cout << "KREPE output verification PASSED\n";
  return 0;
}
catch (const std::exception& error)
{
  std::cerr << "Replay verification failed: " << error.what() << '\n';
  return 1;
}
