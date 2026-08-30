#pragma once

#include "user_data/FieldAccessor.h"

#include <cstdint>

namespace dyablo::HyperbolicEulerNoPatchReplay {

using CellIndex = ForeachCell::CellIndex;
using CellMetaData = ForeachCell::CellMetaData;
using FieldAccessor = UserData_Impl::UserData_FieldAccessor_impl<false>;

template<typename Policy>
struct Data
{
  ForeachCell::IterationSpace_fullArray iter_space;
  CellMetaData cellmetadata;
  bool slope_enabled;
  Policy policy;
  FieldAccessor Uin;
  int ndim;
  real_t dt;
  FieldAccessor Uout;

  uint32_t cell_count() const
  {
    return iter_space.bx() * iter_space.by() * iter_space.bz()
         * iter_space.iOct_count();
  }

  KOKKOS_INLINE_FUNCTION
  CellIndex get_cell(uint32_t index) const
  {
    const uint32_t bx = iter_space.bx();
    const uint32_t by = iter_space.by();
    const uint32_t nbCellsPerBlock = bx * by * iter_space.bz();
    const uint32_t iOct = index / nbCellsPerBlock;
    index %= nbCellsPerBlock;

    const uint32_t k = index / (bx * by);
    const uint32_t j = (index - k * bx * by) / bx;
    const uint32_t i = index - j * bx - k * bx * by;
    return iter_space.getCellIndex(iOct, i, j, k);
  }
};

} // namespace dyablo::HyperbolicEulerNoPatchReplay
