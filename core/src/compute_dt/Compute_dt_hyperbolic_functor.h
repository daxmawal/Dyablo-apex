#pragma once

#include "compute_dt/Compute_dt_base.h"

#include "hyperbolic/policy/HyperbolicPolicy_GLMMHD.h"
#include "hyperbolic/policy/HyperbolicPolicy_Hydro.h"

#include <type_traits>

namespace dyablo {

template <typename Policy>
class ComputeDtHyperbolicFunctor
{
public:
  using CellIndex = ForeachCell::CellIndex;
  using PrimState = typename Policy::PrimState;
  using ConsState = typename Policy::ConsState;

  ComputeDtHyperbolicFunctor( int ndim,
                              real_t gamma0,
                              const ForeachCell::CellMetaData& cells,
                              const UserData::FieldAccessor& Uin,
                              const Policy& policy )
    : ndim(ndim),
      gamma0(gamma0),
      cells(cells),
      Uin(Uin),
      policy(policy)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const CellIndex& iCell, real_t& inv_dt_update ) const
  {
    auto cell_size = cells.getCellSize(iCell);
    real_t dx = cell_size[IX];
    real_t dy = cell_size[IY];
    real_t dz = cell_size[IZ];

    ConsState uLoc = policy.getConsState(Uin, iCell);
    PrimState qLoc = policy.consToPrim(uLoc);

    const real_t cs = sqrt(qLoc.p * gamma0 / qLoc.rho);

    real_t vx = cs + FABS(qLoc.u);
    real_t vy = cs + FABS(qLoc.v);
    real_t vz = (ndim==2)? 0 : cs + FABS(qLoc.w);

    inv_dt_update = FMAX( inv_dt_update, vx/dx + vy/dy + vz/dz );

    // TODO : Find a BETTER way to do this !
    if constexpr (std::is_same_v<PrimState, HyperbolicPolicy_PrimGLMMHDState>) {
      const real_t Bx = qLoc.Bx;
      const real_t By = qLoc.By;
      const real_t Bz = qLoc.Bz;
      const real_t gr = cs*cs*qLoc.rho;
      const real_t Bt2 [] = {By*By+Bz*Bz,
                             Bx*Bx+Bz*Bz,
                             Bx*Bx+By*By};
      const real_t B2 = Bx*Bx + By*By + Bz*Bz;
      const real_t cf1 = gr-B2;
      const real_t V [] = {qLoc.u, qLoc.v, qLoc.w};
      const real_t D [] = {dx, dy, dz};

      real_t cmax = 0.0;
      for (int i=0; i < ndim; ++i) {
        const real_t cf2 = gr + B2 + sqrt(cf1*cf1 + 4.0*gr*Bt2[i]);
        const real_t cf = sqrt(0.5 * cf2 / qLoc.rho);

        cmax += (cf + Kokkos::abs(V[i])) / D[i];
      }
      inv_dt_update = FMAX(inv_dt_update, cmax);
    }
  }

private:
  int ndim;
  real_t gamma0;
  ForeachCell::CellMetaData cells;
  UserData::FieldAccessor Uin;
  Policy policy;
};

} // namespace dyablo
