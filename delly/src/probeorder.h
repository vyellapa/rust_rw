#ifndef PROBEORDER_H
#define PROBEORDER_H

#include <algorithm>
#include <tuple>

namespace torali {

  // Probe generation visits chromosomes, SV IDs, then breakpoint sides. Restore
  // that serial insertion order after parallel collection, then use the original
  // position-only std::sort. This preserves the single-threaded tie permutation.
  template<typename TRegions, typename TSVs>
  inline void sortProbeRegions(TRegions& regions, TSVs const& svs) {
    std::sort(regions.begin(), regions.end(), [&](auto const& a, auto const& b) {
      return std::make_tuple(svs[a.id].chr, a.id, a.bpPoint) <
             std::make_tuple(svs[b.id].chr, b.id, b.bpPoint);
    });
    std::sort(regions.begin(), regions.end());
  }

}

#endif
