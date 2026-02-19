#include <cstdlib>
#include <fstream>
#include <iostream>
#include <simgrid/s4u.hpp>
#include <string>
#include <vector>

struct PlatformConfig {
  int num_hosts;
  std::string net_bw;
  std::string net_lat;
  std::string gpu_bw;
  std::string gpu_lat;
  std::string hostfile_path;

  static std::string require_env_var(const char *name) {
    const char *val = std::getenv(name);
    if (!val) {
      std::cerr << "Error: Environment variable '" << name
                << "' is required but not set." << std::endl;
      std::exit(1);
    }
    return std::string(val);
  }

  static PlatformConfig load() {
    PlatformConfig cfg;
    cfg.num_hosts = std::stoi(require_env_var("PLATFORM_NUM_HOSTS"));
    cfg.net_bw = require_env_var("PLATFORM_NET_BW");
    cfg.net_lat = require_env_var("PLATFORM_NET_LAT");
    cfg.gpu_bw = require_env_var("PLATFORM_GPU_BW");
    cfg.gpu_lat = require_env_var("PLATFORM_GPU_LAT");
    cfg.hostfile_path = require_env_var("PLATFORM_HOSTFILE");
    return cfg;
  }
};

extern "C" void load_platform(simgrid::s4u::Engine &e) {
  PlatformConfig cfg = PlatformConfig::load();
  std::cout << "[S4U Plugin] Loading Platform: " << cfg.num_hosts << " hosts."
            << std::endl;

  auto *root = e.get_netzone_root();
  simgrid::s4u::NetZone *zone = nullptr;

  if (root) {
    zone = root->add_netzone_floyd("world");
  } else {
    zone = e.get_netzone_root()->add_netzone_floyd("world");
  }

  auto *switch_router = zone->add_router("switch");

  for (int i = 0; i < cfg.num_hosts; ++i) {
    std::string host_name = "host-" + std::to_string(i);
    std::string gpu_name = "gpu-" + std::to_string(i);
    std::string net_link_name = "link-net-" + std::to_string(i);
    std::string gpu_link_name = "link-gpu-" + std::to_string(i);

    auto *host = zone->add_host(host_name, "1f");
    auto *gpu = zone->add_host(gpu_name, "1f");

    // Network (Host <-> Switch)
    auto *net_link =
        zone->add_link(net_link_name, cfg.net_bw)->set_latency(cfg.net_lat);

    std::vector<const simgrid::s4u::Link *> net_route = {net_link};
    zone->add_route(host->get_netpoint(), switch_router, net_route);

    // GPU Connection (PCIe)
    auto *gpu_link =
        zone->add_link(gpu_link_name, cfg.gpu_bw)
            ->set_latency(cfg.gpu_lat)
            ->set_sharing_policy(simgrid::s4u::Link::SharingPolicy::SHARED);

    std::vector<simgrid::s4u::LinkInRoute> gpu_route;
    gpu_route.push_back(simgrid::s4u::LinkInRoute(gpu_link));

    zone->add_route(host->get_netpoint(), gpu->get_netpoint(), gpu_route, true);
  }

  zone->seal();
}

int main(int argc, char *argv[]) {
  PlatformConfig cfg = PlatformConfig::load();

  std::cout << "[Generator] Writing hostfile to: " << cfg.hostfile_path
            << std::endl;
  std::ofstream hf(cfg.hostfile_path);
  if (!hf.is_open()) {
    std::cerr << "Error: Could not open " << cfg.hostfile_path
              << " for writing." << std::endl;
    return 1;
  }

  for (int i = 0; i < cfg.num_hosts; ++i) {
    hf << "host-" << i << std::endl;
  }
  hf.close();

  std::cout << "[Generator] Generated " << cfg.num_hosts << " hosts."
            << std::endl;
  return 0;
}
