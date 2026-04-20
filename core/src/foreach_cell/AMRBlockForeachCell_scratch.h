#pragma once

#include "foreach_cell/AMRBlockForeachCell.h"

#define PATCH_LAMBDA KOKKOS_LAMBDA
#define CELL_LAMBDA [&]
#define SCRATCH_LEVEL 1

namespace dyablo {


namespace AMRBlockForeachCell_scratch_impl{

using namespace AMRBlockForeachCell_CellArray_impl;
using CData = AMRBlockForeachCell_CData;
using policy_t = Kokkos::TeamPolicy<>;

class PatchManager;

class CellArray_patch : public CellArray_base<false,false>
{
public:
  using Ref = CellArray_patch;

  KOKKOS_INLINE_FUNCTION
  CellArray_patch(){};
  
  KOKKOS_INLINE_FUNCTION
  CellArray_patch( const Shape_t& s )
  : CellArray_base( s )
  {}
};

/**
 * Represents a group of cells to be iterated upon
 * This is an abstract interface for hyerarchical parallelism on cell arrays 
 * It enables the allocation of temporary arrays to store intermediate results 
 * for the current patch
 **/
class Patch{
  friend PatchManager;
public:
  struct PData{
    const CData cdata;
    uint32_t iOct;
    policy_t::member_type team_member;
  };

private: 
  PData pdata;

public:
  KOKKOS_INLINE_FUNCTION
  Patch( const PData& pdata )
    : pdata(pdata)
  {}

  /**
   * Apply the user-defined function f to every cell of the patch
   * @param iter_space : the iCell parameter in f will take every valid position inside iter_space
   * @param f : a const CellIndex& iCell -> void function 
   *            This is usually a lambda that class CellArray_patch;reads and modify CellArrays at position iCell
   **/
  template <typename Function>
  KOKKOS_INLINE_FUNCTION
  void foreach_cell(const CellArray_shape& iter_space, const Function& f) const
  {
    uint32_t bx = iter_space.bx;
    uint32_t by = iter_space.by;
    uint32_t bz = iter_space.bz;

    uint32_t iOct = pdata.iOct;
    Kokkos::parallel_for( Kokkos::TeamThreadRange(pdata.team_member, bz*by),
      [&]( uint32_t kj )
    {
      uint32_t k = kj/by;
      uint32_t j = kj%by;
      Kokkos::parallel_for( Kokkos::ThreadVectorRange(pdata.team_member, bx),
      [&]( uint32_t i )
      {
        CellIndex iCell = {{iOct,false}, i, j, k, bx, by, bz};
        f( iCell );
      });
    });

    pdata.team_member.team_barrier();
  }

  template <typename Function>
  KOKKOS_INLINE_FUNCTION
  void foreach_cell(const char* /*kernel_name*/, const CellArray_shape& iter_space, const Function& f) const
  {
    foreach_cell(iter_space, f);
  }

  KOKKOS_INLINE_FUNCTION
  CellArray_patch allocate_tmp( const CellArray_patch::Ref& array_ref ) const
  {
    CellArray_patch res = array_ref;
    const auto& s = array_ref.shape;
    res.U = CellArray_patch::View_t(this->pdata.team_member.team_scratch(SCRATCH_LEVEL), s.bx*s.by*s.bz, s.nbFields, 1);
    return res;
  }
};

class PatchManager{
public:
  using Patch = AMRBlockForeachCell_scratch_impl::Patch;
  using CellArray_patch = AMRBlockForeachCell_scratch_impl::CellArray_patch;

private:
  const CData cdata;
  const AMRmesh& pmesh;
  uint32_t scratch_size = 0;

public:
  PatchManager(const CData& cdata, const AMRmesh& pmesh)
  : cdata(cdata), pmesh(pmesh)
  {}

  CellArray_patch::Ref reserve_patch_tmp(std::string name, int gx, int gy, int gz, int nvars)
  {
    const CData& cdata = this->cdata;
    uint32_t bx = cdata.bx+2*gx;
    uint32_t by = cdata.by+2*gy;
    uint32_t bz = cdata.bz+2*gz;
    DYABLO_ASSERT_HOST_RELEASE( cdata.ndim != 2 || bz==1, "bz should be 1 in 2D" );

    scratch_size += CellArray_patch::View_t::shmem_size(bx*by*bz, nvars, 1);
    CellArray_patch::Shape_t shape{bx, by, bz, (uint32_t)nvars, 1};
    return CellArray_patch::Ref(shape);
  }  

  auto make_patch_policy(uint32_t nbOcts, const std::string& kernel_name) const
  {
    auto make_fixed_policy = [&](int team_size, int vector_length)
    {
      return policy_t(nbOcts, team_size, vector_length)
        .set_scratch_size(SCRATCH_LEVEL, Kokkos::PerTeam(this->scratch_size));
    };

    // Fixed CUDA launch settings replayed from apex_converged_tuning.0.yaml.
    if( kernel_name == "HydroE::Patch_octs_4k_8k" )
      return make_fixed_policy(16, 2);

    if( kernel_name == "HydroE::Patch_octs_8k_16k" )
      return make_fixed_policy(16, 2);

    if( kernel_name == "HydroE::Patch_octs_16k_32k" )
      return make_fixed_policy(16, 4);

    if( kernel_name == "HydroE::Patch_octs_32k_64k" )
      return make_fixed_policy(16, 4);

    return policy_t(nbOcts, Kokkos::AUTO(), Kokkos::AUTO())
      .set_scratch_size(SCRATCH_LEVEL, Kokkos::PerTeam(this->scratch_size));
  }
  
  template <typename Function>
  void foreach_patch(const std::string& kernel_name, const Function& f)
  {
    const CData& cdata = this->cdata;
    uint32_t nbOcts = pmesh.getNumOctants();
    auto policy = make_patch_policy(nbOcts, kernel_name);

    Kokkos::parallel_for( kernel_name,
      policy,
      KOKKOS_LAMBDA( policy_t::member_type team_member )
    {
      uint32_t iOct = team_member.league_rank();
      f( Patch({cdata, iOct, team_member}) ); 
    });

    scratch_size = 0;
  }
};

} // namespace AMRBlockForeachCell_scratch_impl

using AMRBlockForeachCell_scratch = AMRBlockForeachCell_impl<AMRBlockForeachCell_scratch_impl::PatchManager>;


} // namespace dyablo
