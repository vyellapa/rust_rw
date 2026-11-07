#include <cassert>
#include <cstdint>
#include <random>
#include <tuple>
#include <vector>
#include "../../delly/src/probeorder.h"

struct SV { int32_t chr; };
struct Region {
  int32_t bppos;
  uint32_t id;
  uint8_t bpPoint;
  bool operator<(Region const& other) const { return bppos < other.bppos; }
  bool operator==(Region const& other) const {
    return std::tie(bppos, id, bpPoint) == std::tie(other.bppos, other.id, other.bpPoint);
  }
};

int main() {
  std::mt19937 random(42);
  for (uint32_t n : {1u, 2u, 12u, 13u, 24u, 128u, 256u}) {
    std::vector<SV> svs;
    for (uint32_t id = 0; id < n; ++id) svs.push_back({int32_t(id % 3)});
    std::vector<Region> serial;
    for (int32_t chr = 0; chr < 3; ++chr)
      for (uint32_t id = 0; id < n; ++id)
        if (svs[id].chr == chr)
          for (uint8_t side = 0; side < 2; ++side)
            serial.push_back({int32_t(100 + id % 2), id, side});
    auto expected = serial;
    std::sort(expected.begin(), expected.end());
    for (int iteration = 0; iteration < 100; ++iteration) {
      auto parallel = serial;
      std::shuffle(parallel.begin(), parallel.end(), random);
      torali::sortProbeRegions(parallel, svs);
      assert(parallel == expected);
    }
  }
}
