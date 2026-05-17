#pragma once

#include "HyperbolicUpdate_base.h"
#include "mpi/GhostCommunicator_partial_blocks.h"
#include "foreach_cell/ForeachCell_utils.h"

namespace dyablo {
namespace{
using CellIndex     = ForeachCell::CellIndex;
using FieldAccessor = UserData::FieldAccessor;
using offset_t      = typename CellIndex::offset_t;
using PatchArray    = ForeachCell::CellArray_patch;

enum VarIndex_gravity {IGX, IGY, IGZ};

template<ComponentIndex3D DIR, int NDIM, typename Policy>
KOKKOS_INLINE_FUNCTION
typename Policy::PrimState get_slope_dir(
  const Policy& policy,
  const PatchArray& SlopesX,
  const PatchArray& SlopesY,
  const PatchArray& SlopesZ,
  const CellIndex& iCell_tmp)
{
  static_assert(DIR == IX || DIR == IY || DIR == IZ);
  static_assert(NDIM == 2 || NDIM == 3);

  if constexpr (DIR == IX)
    return policy.getPrimState( SlopesX, iCell_tmp );
  else if constexpr (DIR == IY)
    return policy.getPrimState( SlopesY, iCell_tmp );
  else
  {
    static_assert(NDIM == 3);
    return policy.getPrimState( SlopesZ, iCell_tmp );
  }
}

template<ComponentIndex3D DIR, int NDIM, typename Policy>
KOKKOS_INLINE_FUNCTION
typename Policy::ConsState process_dir(
  const CellIndex& iCell_U,
  const Policy& policy,
  const FieldAccessor& Uin,
  const FieldAccessor& Uout,
  const PatchArray& SlopesX,
  const PatchArray& SlopesY,
  const PatchArray& SlopesZ,
  const PatchArray& HalfStep,
  const ForeachCell::CellMetaData& cellmetadata,
  const ForeachCell::SearchMode_local& search_local,
  const ForeachCell::SearchMode_neighbor& search_neighbor,
  real_t dt)
{
  using PrimState = typename Policy::PrimState;
  using ConsState = typename Policy::ConsState;

  static_assert(DIR == IX || DIR == IY || DIR == IZ);
  static_assert(NDIM == 2 || NDIM == 3);
  static_assert(DIR != IZ || NDIM == 3);

  // Getting centered value and slope
  CellIndex iCell_tmp = HalfStep.getShape().convert_index( iCell_U, search_local );
  PrimState slope_C = get_slope_dir<DIR, NDIM>(policy, SlopesX, SlopesY, SlopesZ, iCell_tmp);
  PrimState qC_half = policy.getPrimState( HalfStep, iCell_tmp );
  auto size_C = cellmetadata.getCellSize(iCell_U);

  constexpr real_t dim_fac = (NDIM == 2 ? 0.5 : 0.25);
  ConsState du_dir {};
  const real_t fac_C = dt / size_C[DIR];

  // Compute left side flux
  {
    ConsState flux {};
    PrimState qC = qC_half - 0.5 * slope_C;

    offset_t off_m{};
    off_m[DIR] = -1;
    const CellIndex iCell_m_U = iCell_U.getNeighbor(off_m, search_neighbor);
    if( iCell_m_U.is_boundary() )
    {
      flux = policy.getBoundaryFlux(Uin, iCell_m_U, qC, cellmetadata);
    }
    else
    {
      int Ldiff = iCell_m_U.level_diff();
      if (Ldiff >= 0)
      {
        const CellIndex iCell_m_tmp = iCell_tmp + off_m;
        PrimState slope_L = get_slope_dir<DIR, NDIM>(policy, SlopesX, SlopesY, SlopesZ, iCell_m_tmp);

        PrimState qL_half = policy.getPrimState( HalfStep, iCell_m_tmp );

        // Reconstructing
        PrimState qL = qL_half + 0.5 * slope_L;

        // Solving
        flux = policy.riemann_solver(qL, qC, DIR);

        // Adding flux to the neighbor if it is bigger
        if (Ldiff == 1)
        {
          auto size_L = cellmetadata.getCellSize(iCell_m_U);
          ConsState du_n = flux * - dim_fac * dt / size_L[DIR];
          policy.atomic_addConsState(Uout, iCell_m_U, du_n);
        }
      } // If smaller we skip
    }
    du_dir += flux * fac_C;
  }

  // Compute right side flux
  {
    ConsState flux {};
    PrimState qC = qC_half + 0.5 * slope_C;

    offset_t off_p{};
    off_p[DIR] = 1;
    const CellIndex iCell_p_U = iCell_U.getNeighbor(off_p, search_neighbor);
    if( iCell_p_U.is_boundary() )
    {
      flux = policy.getBoundaryFlux(Uin, iCell_p_U, qC, cellmetadata);
    }
    else
    {
      int Rdiff = iCell_p_U.level_diff();
      if (Rdiff >= 0)
      {
        const CellIndex iCell_p_tmp = iCell_tmp + off_p;
        PrimState slope_R = get_slope_dir<DIR, NDIM>(policy, SlopesX, SlopesY, SlopesZ, iCell_p_tmp);

        PrimState qR_half = policy.getPrimState( HalfStep, iCell_p_tmp );

        // Reconstructing
        PrimState qR = qR_half - 0.5 * slope_R;

        // Solving
        flux = policy.riemann_solver(qC, qR, DIR);

        // Adding flux to the neighbor if it is bigger
        if (Rdiff == 1)
        {
          auto size_R = cellmetadata.getCellSize(iCell_p_U);
          ConsState du_n = flux * dim_fac * dt / size_R[DIR];
          policy.atomic_addConsState(Uout, iCell_p_U, du_n);
        }
      }
    }
    du_dir -= flux * fac_C;
  }

  return du_dir;
}

}// namespace

template<typename Policy>
class Hyperbolic_hancock : public HyperbolicUpdate {
  static_assert( is_HyperbolicPolicy_v<Policy>,
  "Policy must be wrapped in HyperbolicPolicy_base");

public:
  using PrimState = typename Policy::PrimState;
  using ConsState = typename Policy::ConsState;

public:
  Hyperbolic_hancock(
          ConfigMap& configMap,
          ForeachCell& foreach_cell,
          Timers& timers) 
  : foreach_cell(foreach_cell),
    timers(timers),
    policy_params(Policy::getParams(configMap)),
    ndim(configMap.getValue<int>("mesh", "ndim", 3)),
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
  void update( UserData& U, ScalarSimulationData& scalar_data)
  {
    real_t dt = scalar_data.get<real_t>("dt");
    real_t gamma0 = this->gamma0;
    int ndim = this->ndim;

    const Policy policy( this->policy_params, scalar_data ); 
    Timers& timers = this->timers; 
    ForeachCell& foreach_cell = this->foreach_cell;

    FieldAccessor Uin = policy.getUin(U);
    FieldAccessor Uout = policy.getUout(U);
    
    timers.get("Hyperbolic_hancock").start();

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

    // Create abstract temporary ghosted arrays for patches 
    uint32_t nbFields = State_traits<PrimState>::nvars;
    PatchArray::Ref Qpatch_ = foreach_cell.reserve_patch_tmp("Qpatch", 2, 2, (ndim == 3)?2:0, nbFields);
    PatchArray::Ref HalfStep_ = foreach_cell.reserve_patch_tmp("HalfStep", 1, 1, (ndim==3)?1:0, nbFields);
    PatchArray::Ref SlopesX_ = foreach_cell.reserve_patch_tmp("SlopesX", 1, 1, (ndim==3)?1:0, nbFields);
    PatchArray::Ref SlopesY_ = foreach_cell.reserve_patch_tmp("SlopesY", 1, 1, (ndim==3)?1:0, nbFields);
    PatchArray::Ref SlopesZ_;
    if( ndim == 3 )
      SlopesZ_ = foreach_cell.reserve_patch_tmp("SlopesZ", 1, 1, 1, nbFields);

    ForeachCell::SearchMode_neighbor search_neighbor( this->foreach_cell.get_amr_mesh().getLightOctree(), ForeachCell::SearchMode_neighbor::ORIGIN );
    ForeachCell::SearchMode_local search_local( ForeachCell::SearchMode_local::ASSERT );

    // Iterate over cells
    foreach_cell.foreach_patch( "Hyperbolic_euler::update",
      PATCH_LAMBDA(const ForeachCell::Patch& patch)
    {
      PatchArray Qpatch = patch.allocate_tmp(Qpatch_);

      patch.foreach_cell( Qpatch, 
        CELL_LAMBDA( const CellIndex& iCell_Qpatch )
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

      PatchArray SlopesX = patch.allocate_tmp(SlopesX_);
      PatchArray SlopesY = patch.allocate_tmp(SlopesY_);
      PatchArray SlopesZ;
      if( ndim == 3 )
        SlopesZ = patch.allocate_tmp(SlopesZ_);
      PatchArray HalfStep = patch.allocate_tmp(HalfStep_);

      patch.foreach_cell( HalfStep.getShape(),
        CELL_LAMBDA(const CellIndex& iCell_tmp)
      {
        // Return Slope at position iCell
        auto compute_slope = [&](const CellIndex &iCell_Uin, const CellIndex &iCell_Qpatch, ComponentIndex3D dir) 
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
          real_t dL = 1, dR = 1;

          // Computing minmod slope for the direction
          PrimState slope = policy.compute_slope( qL, qC, qR, dL, dR);
          return slope;
        }; // compute_slope

        auto compute_half_step = [&]  ( PrimState q,
                              PrimState sx,
                              PrimState sy,
                              PrimState sz,
                              real_t dtdx, real_t dtdy, real_t dtdz )
        {
          // retrieve variations = dx * slopes
          sx*=0.5;
          sy*=0.5;
          sz*=0.5;

          PrimState half_step{};
          if( ndim == 3 )
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
        ForeachCell::SearchMode_neighbor search_neighbor_origin( cellmetadata.getLightOctree(), ForeachCell::SearchMode_neighbor::ORIGIN );
        CellIndex iCell_Uin = Uin.getShape().convert_index(iCell_tmp, search_neighbor_origin);

        if( iCell_Uin.is_valid() && iCell_Uin.level_diff() >= 0 )
        { //Compute slopes only inside of domain and skip smaller neighbors
          CellIndex iCell_Qpatch = Qpatch.getShape().convert_index(iCell_tmp, search_neighbor_origin);

          const PrimState q = policy.getPrimState(Qpatch, iCell_Qpatch );

          auto size = cellmetadata.getCellSize(iCell_Uin);

          PrimState sx = compute_slope(iCell_Uin, iCell_Qpatch, IX);
          PrimState sy = compute_slope(iCell_Uin, iCell_Qpatch, IY);
          PrimState sz {};
          if(ndim == 3)
            sz = compute_slope(iCell_Uin, iCell_Qpatch, IZ);

          PrimState q_half = compute_half_step( q, 
                                        sx, sy, sz, 
                                        dt/size[IX], dt/size[IY], dt/size[IZ]);

          policy.setPrimState( SlopesX, iCell_tmp, sx );
          policy.setPrimState( SlopesY, iCell_tmp, sy );
          if(ndim == 3)
            policy.setPrimState( SlopesZ, iCell_tmp, sz );
          policy.setPrimState( HalfStep, iCell_tmp, q_half );
        }
      });

      patch.foreach_cell( Uout.getShape(),
        CELL_LAMBDA(const CellIndex& iCell)
      {
        ForeachCell::SearchMode_neighbor search_neighbor( cellmetadata.getLightOctree(), ForeachCell::SearchMode_neighbor::CLOSEST );

        ConsState du{};
        if (ndim == 3)
        {
          du += process_dir<IX, 3>(iCell, policy, Uin, Uout, SlopesX, SlopesY, SlopesZ, HalfStep, cellmetadata, search_local, search_neighbor, dt);
          du += process_dir<IY, 3>(iCell, policy, Uin, Uout, SlopesX, SlopesY, SlopesZ, HalfStep, cellmetadata, search_local, search_neighbor, dt);
          du += process_dir<IZ, 3>(iCell, policy, Uin, Uout, SlopesX, SlopesY, SlopesZ, HalfStep, cellmetadata, search_local, search_neighbor, dt);
        }
        else
        {
          du += process_dir<IX, 2>(iCell, policy, Uin, Uout, SlopesX, SlopesY, SlopesZ, HalfStep, cellmetadata, search_local, search_neighbor, dt);
          du += process_dir<IY, 2>(iCell, policy, Uin, Uout, SlopesX, SlopesY, SlopesZ, HalfStep, cellmetadata, search_local, search_neighbor, dt);
        }
        policy.atomic_addConsState(Uout, iCell, du);
      });     
    });

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

    timers.get("Hyperbolic_hancock").stop();
  }

private:
  ForeachCell& foreach_cell;
  
  Timers& timers;  
  typename Policy::Params policy_params;

  int ndim;
  real_t gamma0, smallr, smallp;
};

} // namespace dyablo
