#pragma once

#include "HyperbolicUpdate_base.h"
#include "mpi/GhostCommunicator_partial_blocks.h"
#include "foreach_cell/ForeachCell_utils.h"
#if defined(DYABLO_USE_KKF_REPLAY_FUNCTOR)
#include "Hyperbolic_euler_nopatch_replay.h"
#include <krepe/extractor.hpp>
#endif

namespace dyablo {
namespace{
using CellIndex     = ForeachCell::CellIndex;
using FieldAccessor = UserData::FieldAccessor;
using offset_t      = typename CellIndex::offset_t;

enum VarIndex_gravity {IGX, IGY, IGZ};

}// namespace

/**
 * @brief Euler update algorithm
 * 
 * @tparam State the type of state to treat
 */
template<typename Policy>
class HyperbolicUpdate_euler_nopatch : public HyperbolicUpdate {
  static_assert( is_HyperbolicPolicy_v<Policy>,
  "Policy must be wrapped in HyperbolicPolicy_base");

public:
  using PrimState = typename Policy::PrimState;
  using ConsState = typename Policy::ConsState;

public:
  HyperbolicUpdate_euler_nopatch(
          ConfigMap& configMap,
          ForeachCell& foreach_cell,
          Timers& timers) 
  : foreach_cell(foreach_cell),
    timers(timers),
    policy_params(Policy::getParams(configMap)),
    ndim(configMap.getValue<int>("mesh", "ndim", 3)),
    smallr( configMap.getValue<real_t>("hydro","smallr", 1e-10) ),
    smallp( configMap.getValue<real_t>("hydro","smallp", 1e-10) ),
    slope_enabled( configMap.getValue<bool>("hydro","slope_enabled", true) )
  { }

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

    FieldAccessor Uin = policy.getUin(U);
    FieldAccessor Uout = policy.getUout(U);
    
    timers.get("HyperbolicUpdate_euler").start();

    ForeachCell::CellMetaData cellmetadata = foreach_cell.getCellMetaData();

    // Initializing output array 
    // TODO : remove this and copy Uin->Uout in timeloop or field creation logic
    foreach_cell.foreach_cell( "Hyperbolic_euler::init",
      Uout.getShape(),
      KOKKOS_LAMBDA(const CellIndex &iCell) 
    {
      ConsState uC = policy.getConsState(Uin, iCell);
      policy.setConsState(Uout, iCell, uC);
    });

    // Setting the ghosts to 0 to accumulate fluxes
    foreach_cell.foreach_ghost_cell( "Hyperbolic_euler::resetting_ghosts",
      Uout.getShape(),
      KOKKOS_LAMBDA(const CellIndex &iCell) 
    {
      ConsState empty_state{};
      policy.setConsState(Uout, iCell, empty_state);
    });

    // Iterate over cells
#if defined(DYABLO_USE_KKF_REPLAY_FUNCTOR)
    using ReplayData = HyperbolicEulerNoPatchReplay::Data<Policy>;
    ReplayData replay_data{
      ForeachCell::IterationSpace_fullArray{Uout.getShape()},
      cellmetadata, slope_enabled, policy, Uin, ndim, dt, Uout};

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
#else
    foreach_cell.foreach_cell( "Hyperbolic_euler::update",
      Uout.getShape(),
      KOKKOS_LAMBDA(const CellIndex &iCell)
    {
#endif
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
      
#if defined(DYABLO_USE_KKF_REPLAY_FUNCTOR)
    };
    krepe::parallel_for(
      "Hyperbolic_euler::update",
      Kokkos::RangePolicy<>(0, replay_data.cell_count()),
      update_functor);
#else
    });
#endif

    // Reducing the ghosts to accumulate the flux in the data arrays 
    int ghost_count = 1;
    GhostCommunicator_partial_blocks ghost_comm ( 
      foreach_cell.get_amr_mesh(),
      Uout.getShape(),
      ghost_count );
    ghost_comm.reduce_ghosts( Uout );

    if constexpr ( Policy::has_postProcess() )
    {
      foreach_cell.foreach_cell( "HyperbolicUpdate::post-process", Uout.getShape(),
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
  ForeachCell& foreach_cell;
  
  Timers& timers;  
  typename Policy::Params policy_params;

  int ndim;
  real_t smallr, smallp;
  bool slope_enabled;
};

} // namespace dyablo

