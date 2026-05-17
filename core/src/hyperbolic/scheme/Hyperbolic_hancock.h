#pragma once

#include "HyperbolicUpdate_base.h"
#include "mpi/GhostCommunicator_partial_blocks.h"
#include "foreach_cell/ForeachCell_utils.h"
#include "utils/misc/Dyablo_assert.h"

namespace dyablo {
namespace{
using CellIndex     = ForeachCell::CellIndex;
using FieldAccessor = UserData::FieldAccessor;
using offset_t      = typename CellIndex::offset_t;

enum VarIndex_gravity {IGX, IGY, IGZ};

}// namespace

template<typename Policy>
class Hyperbolic_hancock : public HyperbolicUpdate {
  static_assert( is_HyperbolicPolicy_v<Policy>,
  "Policy must be wrapped in HyperbolicPolicy_base");

public:
  using PrimState = typename Policy::PrimState;
  using ConsState = typename Policy::ConsState;
  using ComputeFn = void (Hyperbolic_hancock::*)(UserData&, ScalarSimulationData&);

public:
  Hyperbolic_hancock(
          ConfigMap& configMap,
          ForeachCell& foreach_cell,
          Timers& timers) 
  : foreach_cell(foreach_cell),
    timers(timers),
    policy_params(Policy::getParams(configMap)),
    ndim(configMap.getValue<int>("mesh", "ndim", 3)),
    compute_fn((ndim == 3) ? &Hyperbolic_hancock::compute<3> : &Hyperbolic_hancock::compute<2>),
    gamma0( configMap.getValue<real_t>("hydro","gamma0", 1.4) ),
    smallr( configMap.getValue<real_t>("hydro","smallr", 1e-10) ),
    smallp( configMap.getValue<real_t>("hydro","smallp", 1e-10) )
  { }

  /**
   * @brief Solves hydro for one step using the euler method
   * 
   * @param U the input/output global array
   * @param scalar_data input scalar data
   */

  void update( UserData& U, ScalarSimulationData& scalar_data) override 
  {
    (this->*compute_fn)(U, scalar_data);
  }

  template<int NDim>
  void compute( UserData& U, ScalarSimulationData& scalar_data)
  {
    const real_t dt = scalar_data.get<real_t>("dt");
    const real_t gamma0 = this->gamma0;

    const Policy policy( this->policy_params, scalar_data ); 
    Timers& timers = this->timers; 
    ForeachCell& foreach_cell = this->foreach_cell;

    const FieldAccessor Uin = policy.getUin(U);
    FieldAccessor Uout = policy.getUout(U);
    
    timers.get("Hyperbolic_hancock").start();

    const ForeachCell::CellMetaData cellmetadata = foreach_cell.getCellMetaData();

    // Initializing output array 
    // TODO : remove this and copy Uin->Uout in timeloop or field creation logic
    foreach_cell.foreach_cell( "Hyperbolic_euler::init",
      Uout.getShape(),
      KOKKOS_LAMBDA(const CellIndex &iCell) 
    {
      const ConsState uC = policy.getConsState(Uin, iCell);
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

    // Create abstract temporary ghosted arrays for patches 
    using PatchArray = ForeachCell::CellArray_patch;
    constexpr uint32_t nbFields = State_traits<PrimState>::nvars;
    const PatchArray::Ref Qpatch_ = foreach_cell.reserve_patch_tmp("Qpatch", 2, 2, (NDim == 3)?2:0, nbFields);
    const PatchArray::Ref HalfStep_ = foreach_cell.reserve_patch_tmp("HalfStep", 1, 1, (NDim==3)?1:0, nbFields);
    const PatchArray::Ref SlopesX_ = foreach_cell.reserve_patch_tmp("SlopesX", 1, 1, (NDim==3)?1:0, nbFields);
    const PatchArray::Ref SlopesY_ = foreach_cell.reserve_patch_tmp("SlopesY", 1, 1, (NDim==3)?1:0, nbFields);
    const PatchArray::Ref SlopesZ_ = (NDim == 3) ? foreach_cell.reserve_patch_tmp("SlopesZ", 1, 1, 1, nbFields) : PatchArray::Ref{};

    const ForeachCell::SearchMode_local search_local( ForeachCell::SearchMode_local::ASSERT );

    // Iterate over cells
    foreach_cell.foreach_patch( "Hyperbolic_euler::update",
      PATCH_LAMBDA(const ForeachCell::Patch& patch)
    {
      const PatchArray Qpatch = patch.allocate_tmp(Qpatch_);

      patch.foreach_cell( Qpatch, 
        CELL_LAMBDA( const CellIndex& iCell_Qpatch )
      {
        const ForeachCell::SearchMode_neighbor search_neighbor_origin( cellmetadata.getLightOctree(), ForeachCell::SearchMode_neighbor::ORIGIN );
        const CellIndex iCell_Uin = Uin.getShape().convert_index(iCell_Qpatch, search_neighbor_origin);
        const int level_diff = iCell_Uin.level_diff();
        ConsState u = {};
        if (iCell_Uin.is_boundary())
          u = policy.getBoundaryValue(Uin, iCell_Uin, cellmetadata);
        else if (level_diff < 0) {
          const int subcell_count = 
          foreach_sibling(NDim, iCell_Uin, search_neighbor_origin,
            [&](const CellIndex& iCell_neigh) {
              const ConsState uloc = policy.getConsState(Uin, iCell_neigh);
              u += uloc;
            });
          u /= subcell_count;
        }
        else
          u = policy.getConsState(Uin, iCell_Uin);
        
        const PrimState q = policy.consToPrim( u );
        policy.setPrimState( Qpatch, iCell_Qpatch, q );
      });

      const PatchArray SlopesX = patch.allocate_tmp(SlopesX_);
      const PatchArray SlopesY = patch.allocate_tmp(SlopesY_);
      const PatchArray SlopesZ = (NDim == 3) ? patch.allocate_tmp(SlopesZ_) : PatchArray{};
      const PatchArray HalfStep = patch.allocate_tmp(HalfStep_);

      patch.foreach_cell( HalfStep.getShape(),
        CELL_LAMBDA(const CellIndex& iCell_tmp)
      {
        // Return Slope at position iCell
        const auto compute_slope = [&](const CellIndex &iCell_Uin, const CellIndex &iCell_Qpatch, ComponentIndex3D dir) 
        {        
          const PrimState qC = policy.getPrimState(Qpatch, iCell_Qpatch );
          offset_t off_m{}; off_m[dir] = -1;
          const PrimState qL = policy.getPrimState(Qpatch, iCell_Qpatch + off_m); 
          offset_t off_p{}; off_p[dir] =  1;
          const PrimState qR = policy.getPrimState(Qpatch, iCell_Qpatch + off_p); 
        
          // //!\ Neighbor cells in Qpatch are averaged cells -> size iCell_L != size iCell_Qpatch_L
          // CellIndex iCell_L = iCell_Uin.getNeighbor_ghost(off_m, Uout.getShape());
          // CellIndex iCell_R = iCell_Uin.getNeighbor_ghost(off_p, Uout.getShape());   

          // // Getting the length right and left
          // // Smaller -> use averaged same-size cell -> 1*dx
          // // Bigger -> same-size cell in Qpatch has vame value as actual bigger cell -> dx/2 + dx             
          // constexpr real_t sizes[] = {1.0, 1.0, 1.5}; 
          // const real_t dL = sizes[iCell_L.level_diff()+1];
          // const real_t dR = sizes[iCell_R.level_diff()+1];  
          const real_t dL = 1;
          const real_t dR = 1;

          // Computing minmod slope for the direction
          const PrimState slope = policy.compute_slope( qL, qC, qR, dL, dR);
          return slope;
        }; // compute_slope

        const auto compute_half_step = [&]  ( const PrimState q,
                              PrimState sx,
                              PrimState sy,
                              PrimState sz,
                              const real_t dtdx, const real_t dtdy, const real_t dtdz )
        {
          // retrieve variations = dx * slopes
          sx*=0.5;
          sy*=0.5;
          sz*=0.5;

          PrimState half_step{};
          if( NDim == 3 )
          {
            half_step.rho = q.rho + (-q.u * sx.rho - sx.u * q.rho) * dtdx 
                                  + (-q.v * sy.rho - sy.v * q.rho) * dtdy 
                                  + (-q.w * sz.rho - sz.w * q.rho) * dtdz;
            half_step.u   = q.u + (-q.u * sx.u - sx.p / q.rho) * dtdx 
                                + (-q.v * sy.u) * dtdy 
                                + (-q.w * sz.u) * dtdz;
            half_step.v   = q.v + (-q.u * sx.v) * dtdx 
                                + (-q.v * sy.v - sy.p / q.rho) * dtdy 
                                + (-q.w * sz.v) * dtdz;
            half_step.w   = q.w + (-q.u * sx.w) * dtdx  
                                + (-q.v * sy.w) * dtdy 
                                + (-q.w * sz.w - sz.p / q.rho) * dtdz;
            half_step.p   = q.p + (-q.u * sx.p - sx.u * gamma0 * q.p) * dtdx 
                                + (-q.v * sy.p - sy.v * gamma0 * q.p) * dtdy 
                                + (-q.w * sz.p - sz.w * gamma0 * q.p) * dtdz;
          }
          else
          {
            half_step.rho = q.rho + (-q.u * sx.rho - sx.u * q.rho) * dtdx 
                                  + (-q.v * sy.rho - sy.v * q.rho) * dtdy;
            half_step.u   = q.u + (-q.u * sx.u - sx.p / q.rho) * dtdx 
                                + (-q.v * sy.u) * dtdy;
            half_step.v   = q.v + (-q.u * sx.v) * dtdx 
                                + (-q.v * sy.v - sy.p / q.rho) * dtdy;
            half_step.p   = q.p + (-q.u * sx.p - sx.u * gamma0 * q.p) * dtdx 
                                + (-q.v * sy.p - sy.v * gamma0 * q.p) * dtdy;
          }
          return half_step;
        };
        const ForeachCell::SearchMode_neighbor search_neighbor_origin( cellmetadata.getLightOctree(), ForeachCell::SearchMode_neighbor::ORIGIN );
        const CellIndex iCell_Uin = Uin.getShape().convert_index(iCell_tmp, search_neighbor_origin);

        if( iCell_Uin.is_valid() && iCell_Uin.level_diff() >= 0 )
        { //Compute slopes only inside of domain and skip smaller neighbors
          const CellIndex iCell_Qpatch = Qpatch.getShape().convert_index(iCell_tmp, search_neighbor_origin);

          const PrimState q = policy.getPrimState(Qpatch, iCell_Qpatch );

          const auto size = cellmetadata.getCellSize(iCell_Uin);

          const PrimState sx = compute_slope(iCell_Uin, iCell_Qpatch, IX);
          const PrimState sy = compute_slope(iCell_Uin, iCell_Qpatch, IY);
          const PrimState sz = (NDim == 3) ? compute_slope(iCell_Uin, iCell_Qpatch, IZ) : PrimState{};

          const PrimState q_half = compute_half_step( q, 
                                        sx, sy, sz, 
                                        dt/size[IX], dt/size[IY], dt/size[IZ]);

          policy.setPrimState( SlopesX, iCell_tmp, sx );
          policy.setPrimState( SlopesY, iCell_tmp, sy );
          if(NDim == 3)
            policy.setPrimState( SlopesZ, iCell_tmp, sz );
          policy.setPrimState( HalfStep, iCell_tmp, q_half );
        }
      });

      patch.foreach_cell( Uout.getShape(),
        CELL_LAMBDA(const CellIndex& iCell)
      {
        const ForeachCell::SearchMode_neighbor search_neighbor( cellmetadata.getLightOctree(), ForeachCell::SearchMode_neighbor::CLOSEST );

        const auto process_dir = [&](const CellIndex &iCell_U, ComponentIndex3D dir) {
          const auto get_slope = [&](const CellIndex &iCell_tmp, ComponentIndex3D dir)
          {
            if( dir==IX )
              return policy.getPrimState( SlopesX, iCell_tmp );
            else if( dir==IY )
              return policy.getPrimState( SlopesY, iCell_tmp );
            else
              return policy.getPrimState( SlopesZ, iCell_tmp );
          };       
          
          // Getting centered value and slope
          const CellIndex iCell_tmp = HalfStep.getShape().convert_index( iCell_U, search_local );
          const PrimState slope_C = get_slope(iCell_tmp, dir);       
          const PrimState qC_half = policy.getPrimState( HalfStep, iCell_tmp );
          const auto size_C = cellmetadata.getCellSize(iCell_U);

          constexpr real_t dim_fac = (NDim == 2 ? 0.5 : 0.25);
          ConsState du_dir {};
          const real_t fac_C = dt / size_C[dir];

          // Compute left side flux
          {
            ConsState flux {};
            const PrimState qC = qC_half - 0.5 * slope_C;

            offset_t off_m{}; 
            off_m[dir] = -1;
            const CellIndex iCell_m_U = iCell_U.getNeighbor(off_m, search_neighbor);
            if( iCell_m_U.is_boundary() )
            {
              flux = policy.getBoundaryFlux(Uin, iCell_m_U, qC, cellmetadata);
            }
            else
            {  
              const int Ldiff = iCell_m_U.level_diff();
              if (Ldiff >= 0) 
              {       
                
                const CellIndex iCell_m_tmp = iCell_tmp + off_m;
                const PrimState slope_L = get_slope(iCell_m_tmp, dir);

                const PrimState qL_half = policy.getPrimState( HalfStep, iCell_m_tmp );

                // Reconstructing
                const PrimState qL = qL_half + 0.5 * slope_L;

                // Solving
                flux = policy.riemann_solver(qL, qC, dir);
                
                // Adding flux to the neighbor if it is bigger
                if (Ldiff == 1) 
                {
                  const auto size_L = cellmetadata.getCellSize(iCell_m_U);
                  const ConsState du_n = flux * - dim_fac * dt / size_L[dir];
                  policy.atomic_addConsState(Uout, iCell_m_U, du_n);
                }
              } // If smaller we skip
            }
            du_dir += flux * fac_C;
          }

          // Compute right side flux
          {     
            ConsState flux {};
            const PrimState qC = qC_half + 0.5 * slope_C;

            offset_t off_p{}; 
            off_p[dir] = 1;
            const CellIndex iCell_p_U = iCell_U.getNeighbor(off_p, search_neighbor);
            if( iCell_p_U.is_boundary() )
            {
              flux = policy.getBoundaryFlux(Uin, iCell_p_U, qC, cellmetadata);
            }
            else
            {
              const int Rdiff = iCell_p_U.level_diff();
              if (Rdiff >= 0) 
              {
                const CellIndex iCell_p_tmp = iCell_tmp + off_p;
                const PrimState slope_R = get_slope(iCell_p_tmp, dir);

                const PrimState qR_half = policy.getPrimState( HalfStep, iCell_p_tmp );

                // Reconstructing
                const PrimState qR = qR_half - 0.5 * slope_R;

                // Solving
                flux = policy.riemann_solver(qC, qR, dir);

                // Adding flux to the neighbor if it is bigger
                if (Rdiff == 1)
                {
                  const auto size_R = cellmetadata.getCellSize(iCell_p_U);
                  const ConsState du_n = flux * dim_fac * dt / size_R[dir];
                  policy.atomic_addConsState(Uout, iCell_p_U, du_n);
                }          
              }
            }
            du_dir -= flux * fac_C;
          } 

          return du_dir;
        };

        ConsState du{};
        du += process_dir(iCell, IX);
        du += process_dir(iCell, IY);
        if(NDim == 3)
          du += process_dir(iCell, IZ);
        policy.atomic_addConsState(Uout, iCell, du);
      });     
    });

    // Reducing the ghosts to accumulate the flux in the data arrays 
    const int ghost_count = 1;
    const GhostCommunicator_partial_blocks ghost_comm ( 
      foreach_cell.get_amr_mesh(),
      Uout.getShape(),
      ghost_count );
    ghost_comm.reduce_ghosts( Uout );

    if constexpr ( Policy::has_postProcess() )
    {
      foreach_cell.foreach_cell( "HyperbolicUpdate::post-process", Uout.getShape(),
        KOKKOS_LAMBDA(  const ForeachCell::CellIndex& iCell)
      {
        const ConsState u = policy.getConsState(Uout, iCell);
        const ConsState u_pp = policy.postProcess( u );
        policy.setConsState( Uout, iCell, u_pp );
      });
    }

    policy.printWarnings();

    timers.get("Hyperbolic_hancock").stop();
  }

private:
  ForeachCell& foreach_cell;
  
  Timers& timers;  
  const typename Policy::Params policy_params;

  const int ndim;
  const ComputeFn compute_fn;
  const real_t gamma0, smallr, smallp;
};

} // namespace dyablo
