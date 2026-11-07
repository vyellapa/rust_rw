// Flat C interface over edlib for Rust FFI.
#include <cstdlib>
#include <cstring>
#include "edlib.h"

extern "C" {

  struct DlEdlibResult {
    int status;
    int editDistance;
    int endLocation0;     // endLocations[0] or -1
    int alignmentLength;
    unsigned char* alignment;  // owned, release with dl_edlib_free
  };

  void dl_edlib_align(const char* q, int ql, const char* t, int tl, int k, int mode, int task, DlEdlibResult* out) {
    EdlibAlignResult r = edlibAlign(q, ql, t, tl, edlibNewAlignConfig(k, (EdlibAlignMode) mode, (EdlibAlignTask) task, NULL, 0));
    out->status = r.status;
    out->editDistance = r.editDistance;
    out->endLocation0 = (r.numLocations > 0 && r.endLocations != NULL) ? r.endLocations[0] : -1;
    out->alignmentLength = r.alignmentLength;
    out->alignment = NULL;
    if ((r.alignment != NULL) && (r.alignmentLength > 0)) {
      out->alignment = (unsigned char*) std::malloc(r.alignmentLength);
      std::memcpy(out->alignment, r.alignment, r.alignmentLength);
    }
    edlibFreeAlignResult(r);
  }

  void dl_edlib_free(unsigned char* p) { std::free(p); }
}
