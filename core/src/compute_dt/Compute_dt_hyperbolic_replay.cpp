#include "compute_dt/Compute_dt_hyperbolic_functor.h"

#include <kernel_replayer.hpp>

#include "amr/AMRmesh.h"
#include "utils/config/ConfigMap.h"
#include "utils/mpi/GlobalMpiSession.h"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

constexpr const char* ini_filename = "test_sod_3D.ini";
constexpr std::uint32_t nb_octs = 1796;
constexpr real_t expected_dt = 0.000593442;
constexpr real_t dt_tolerance = 1e-9;

ConfigMap load_config(const std::string& filename) {
  std::ifstream input(filename);
  if (!input) {
    throw std::runtime_error("Could not open .ini file: " + filename);
  }

  std::ostringstream buffer;
  buffer << input.rdbuf();
  return ConfigMap(buffer.str());
}

} // namespace

int main(int argc, char* argv[]) {
  cexa::kernel_replayer::ScopeGuard replay_scope(argc, argv);
  dyablo::GlobalMpiSession mpi_session(&argc, &argv);
  Kokkos::ScopeGuard kokkos_scope(argc, argv);

  ConfigMap configMap = load_config(ini_filename);

  using Policy = dyablo::HyperbolicPolicy_Hydro;
  dyablo::ScalarSimulationData scalar_data;
  typename Policy::Params params = Policy::getParams(configMap);
  Policy policy{params, scalar_data};

  const auto mesh_params = dyablo::AMRmesh::parse_parameters(configMap);
  dyablo::AMRmesh mesh(mesh_params.dim, mesh_params.periodic,
                        mesh_params.level_min, mesh_params.level_max,
                        mesh_params.coarse_grid_size);
  dyablo::ForeachCell foreach_cell(mesh, configMap);

  dyablo::ForeachCell::CellArray_shape shape{
      configMap.getValue<std::uint32_t>("amr", "bx", 0),
      configMap.getValue<std::uint32_t>("amr", "by", 0),
      configMap.getValue<std::uint32_t>("amr", "bz", 1),
      static_cast<std::uint32_t>(Policy::ConsState::getFieldsInfo().size()),
      nb_octs,
      0,
      0,
      0};

  dyablo::UserData::FieldAccessor Uin;
  dyablo::ComputeDtHyperbolicFunctor<Policy> functor(
      foreach_cell.getDim(), params.policy_params.gamma0,
      foreach_cell.getCellMetaData(), Uin, policy);

  real_t inv_dt = 0;
  foreach_cell.reduce_cell("compute_dt", shape, functor,
                            Kokkos::Max<real_t>(inv_dt));
  Kokkos::fence();

  const real_t cfl = configMap.getValue<real_t>("dt", "hydro_cfl", 0.5);
  const real_t dt = cfl / inv_dt;
  if (std::abs(dt - expected_dt) > dt_tolerance) {
    std::cerr << "compute_dt replay mismatch: got dt=" << dt
              << ", expected " << expected_dt
              << " +/- " << dt_tolerance << '\n';
    return 1;
  }

  std::cout << "compute_dt replay complete: cells="
            << static_cast<std::uint64_t>(shape.bx) * shape.by * shape.bz *
                  nb_octs
            << " inv_dt=" << inv_dt << " dt=" << dt << '\n';
  return 0;
}
