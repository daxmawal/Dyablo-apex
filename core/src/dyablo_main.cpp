#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <string>
#include <vector>

#include "DyabloSession.hpp"
#include "utils/config/ConfigMap.h"
#include "DyabloTimeLoop.h"

namespace {

struct DyabloLaunchOptions
{
  std::vector<char*> forwarded_argv;
  std::string input_file;
  bool has_hydro_patch_replay_pass_count = false;
  int hydro_patch_replay_pass_count = 0;
  bool has_hydro_patch_launch_policy = false;
  std::string hydro_patch_launch_policy;
};

bool starts_with(const std::string& value, const std::string& prefix)
{
  return value.rfind(prefix, 0) == 0;
}

std::string to_lower_copy(std::string value)
{
  std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c)
  {
    return static_cast<char>(std::tolower(c));
  });
  return value;
}

void print_usage(const char* program_name)
{
  const char* exe_name = program_name != nullptr ? program_name : "./dyablo";

  std::cout << "Error : no input file" << std::endl;
  std::cout << "Usage:" << std::endl;
  std::cout << "  " << exe_name
            << " [--kokkos-***=*]"
            << " [--dyablo-hydro-patch-replay-pass-count=N]"
            << " [--dyablo-hydro-patch-policy=optimized|auto]"
            << " input_file.ini" << std::endl;
}

int parse_positive_int_option(const std::string& value, const std::string& option_name)
{
  std::size_t parsed_length = 0;
  int parsed_value = 0;
  try
  {
    parsed_value = std::stoi(value, &parsed_length);
  }
  catch(...)
  {
    DYABLO_ASSERT_HOST_RELEASE(false, "Invalid value for " << option_name << " : `" << value << "`");
  }

  DYABLO_ASSERT_HOST_RELEASE(
    parsed_length == value.size(),
    "Invalid value for " << option_name << " : `" << value << "`"
  );
  DYABLO_ASSERT_HOST_RELEASE(
    parsed_value > 0,
    option_name << " must be strictly positive, got `" << value << "`"
  );
  return parsed_value;
}

std::string normalize_hydro_patch_launch_policy(const std::string& value)
{
  const std::string normalized_value = to_lower_copy(value);
  DYABLO_ASSERT_HOST_RELEASE(
    normalized_value == "optimized" || normalized_value == "auto",
    "Invalid value for --dyablo-hydro-patch-policy : `" << value
    << "`. Expected `optimized` or `auto`."
  );
  return normalized_value;
}

DyabloLaunchOptions parse_launch_options(int argc, char** argv)
{
  DyabloLaunchOptions options;
  options.forwarded_argv.reserve(argc + 1);

  if( argc > 0 )
    options.forwarded_argv.push_back(argv[0]);

  for( int i = 1; i < argc; ++i )
  {
    const std::string arg(argv[i]);

    auto read_option_value = [&](const std::string& option_name) -> std::string
    {
      DYABLO_ASSERT_HOST_RELEASE(i + 1 < argc, "Missing value after " << option_name);
      ++i;
      return argv[i];
    };

    if( starts_with(arg, "--dyablo-hydro-patch-replay-pass-count=") )
    {
      const std::string value = arg.substr(std::string("--dyablo-hydro-patch-replay-pass-count=").size());
      options.hydro_patch_replay_pass_count = parse_positive_int_option(
        value,
        "--dyablo-hydro-patch-replay-pass-count"
      );
      options.has_hydro_patch_replay_pass_count = true;
    }
    else if( arg == "--dyablo-hydro-patch-replay-pass-count" )
    {
      const std::string value = read_option_value("--dyablo-hydro-patch-replay-pass-count");
      options.hydro_patch_replay_pass_count = parse_positive_int_option(
        value,
        "--dyablo-hydro-patch-replay-pass-count"
      );
      options.has_hydro_patch_replay_pass_count = true;
    }
    else if( starts_with(arg, "--dyablo-hydro-patch-policy=") )
    {
      const std::string value = arg.substr(std::string("--dyablo-hydro-patch-policy=").size());
      options.hydro_patch_launch_policy = normalize_hydro_patch_launch_policy(value);
      options.has_hydro_patch_launch_policy = true;
    }
    else if( arg == "--dyablo-hydro-patch-policy" )
    {
      const std::string value = read_option_value("--dyablo-hydro-patch-policy");
      options.hydro_patch_launch_policy = normalize_hydro_patch_launch_policy(value);
      options.has_hydro_patch_launch_policy = true;
    }
    else
    {
      options.forwarded_argv.push_back(argv[i]);
      if( options.input_file.empty() && !arg.empty() && arg[0] != '-' )
        options.input_file = arg;
    }
  }

  options.forwarded_argv.push_back(nullptr);
  return options;
}

} // namespace

int main(int argc, char *argv[])
{
  using namespace dyablo;
  auto launch_options = parse_launch_options(argc, argv);
  argc = static_cast<int>(launch_options.forwarded_argv.size()) - 1;
  argv = launch_options.forwarded_argv.data();

  if( launch_options.input_file.empty() )
  {
    print_usage(argc > 0 ? argv[0] : nullptr);
    return EXIT_FAILURE;
  }

  DyabloSession mpi_session(argc, argv);

  /*
   * read parameter file and initialize a ConfigMap object
   */
  ConfigMap configMap = ConfigMap::broadcast_parameters(launch_options.input_file);
  if( launch_options.has_hydro_patch_replay_pass_count )
    configMap.setValue("hydro", "hydro_patch_replay_pass_count", launch_options.hydro_patch_replay_pass_count);
  if( launch_options.has_hydro_patch_launch_policy )
    configMap.setValue("hydro", "hydro_patch_launch_policy", launch_options.hydro_patch_launch_policy);
  if( configMap.hasValue("units","time") )
  { // Set code units
    auto unit_time = configMap.getValue<Units::Time>("units", "time");
    auto unit_length = configMap.getValue<Units::Length>("units", "length");
    Units::Mass unit_mass = Units::kg();
    DYABLO_ASSERT_HOST_RELEASE( !(configMap.hasValue("units","mass") && configMap.hasValue("units","density")), "Parsing units in .ini : cannot set code density and mass et the same time" );
    if( configMap.hasValue("units","density") )
    {
      auto unit_density = configMap.getValue<Units::Density>("units", "density");
      unit_mass = unit_density * unit_length.pow<3>();
    }
    else
    {
      unit_mass = configMap.getValue<Units::Mass>("units", "mass");
    }
    
    Units::code_units_init( Units::UnitSystem(
        unit_time,
        unit_length,
        unit_mass
      ));
  }
  DyabloTimeLoop simulation( configMap );

  simulation.run();

  return EXIT_SUCCESS;

} // end main
