#include <algorithm>
#include <cassert>

#include "mp/flat/redef/encodings.h"

namespace mp {

/// ZZI encoding, Huchette & Vielma '22
class ZZI_Encoding {
public:
  /// Return column k
  std::vector<double>
  GetExtendedColumn(int r, int k, int v0, int v1) {
    if (v1 > size_d()+1) {
      ExtendTo(v1);
      assert(r == size_r());
    }
    std::vector<double> result;
    assert((int)C_.at(k).size() > v1);
    result.assign(C_.at(k).begin()+v0, C_.at(k).begin()+v1+1);
    auto d_this = 1<<r;
    if (d_this+1 == v1)      // v1 is the copied last line
      result.back() = result[result.size()-2];
    return result;
  }


protected:
  /// Return r
  int size_r() const {
    auto r = int(C_.size()-1);
    assert(size_d()==1<<r);
    return r;
  }
  /// Return d==2^r.
  /// We can access up to line d+1
  int size_d() const { return int(C_.at(1).size()-2); }
  /// Extend matrix up to given v (so that d<v)
  void ExtendTo(int v) {
    while (size_d()<v-1)
      Extend_1Duplication();
  }
  /// One binary extension
  void Extend_1Duplication() {
    auto r_old = size_r();
    auto r = r_old+1;
    auto d_old = size_d();
    auto d = d_old*2;
    C_.resize(r+1);
    for (int k=1; k<=r; ++k)
      C_.at(k).resize(d+2);
    for (int k=1; k<r; ++k) {
      auto& C_k = C_.at(k);
      for (int v=d_old+1; v<=d; ++v) {
        C_k[v] = C_k[v-d_old] + C_k[d_old];
      }
    }
    auto& C_r = C_.at(r);
    std::fill(C_r.begin()+d_old+1, C_r.end(), 1.0);
  }


private:
  /// Index 0 not used for first index (k).
  /// Has space for line d+1 for convenience.
  /// Initialize with C^1
  std::vector< std::vector<double> > C_
  { { {}, {0.0, 0.0, 1.0, 0.0} } };
};

P_ZZI_Encoding MakeZZIEncoding() {
  return std::make_shared<ZZI_Encoding>();
}

std::vector<double>
GetExtendedColumn(ZZI_Encoding& zzi,
                  int r, int k, int v0, int v1) {
  return zzi.GetExtendedColumn(r, k, v0, v1);
}

}  // namespace mp
