#ifndef MP_FLAT_CVT_SDP_H
#define MP_FLAT_CVT_SDP_H

#include <cmath>
#include <vector>
#include <string_view>
#include "mp/format.h"
#include <cassert>

#include "mp/flat/constr_keeper.h"
#include "mp/valcvt-link.h"
#include "mp/error.h"

namespace mp {

/**
 * Functor to filter SDP constructs.
 *
 * We support the syntax with function kth_diag() from
 * Benson, H., Vanderbei, R.
 * Solving Problems with Semidefinite and Related Constraints
 * Using Interior-Point Methods for Nonlinear Programming.
 * Math. Program., Ser. B 95, 279–302 (2003).
 *
 * @param MCType the main converter.
 */
template <class MCType>
class SDPConverter : public MCKeeper<MCType> {
public:
  virtual ~SDPConverter() { }

  /// Constructor
  SDPConverter(MCType& mc) :
      MCKeeper<MCType>(mc) { }

  /// Run conversions
  void Run() {
    if (this->MC().GetNumberOfAddable(
            (const CallConstraint*)nullptr)
        && this->MC().template
           UserAcceptsExprForCon<SDPDotProdConstraint>()) {
      ScanCalls();
      ScanCallInvocations();
      ScanDotProducts();
    }
  }

protected:
  /// Filter actuall kth_diag() calls
  void ScanCalls() {
    auto& ck = this->MC().GetConstraintKeeper(
        (CallConstraint*)nullptr );
    ck.ForEachActive(
        [this](const CallConstraint& cc, int ) {
          const auto& func =
              this->MC().GetFunction(cc.GetParameters()[0]);
          if (std::string_view("kth_diag") == func.Name())
            RegisterCall(cc);
          return false;
        });
  }

  /// Check where kth_diag()'s are called from
  /// and, if correct,
  /// disable these containing constraints
  void ScanCallInvocations() {

  }

  /// Filter dot products
  void ScanDotProducts() {

  }

protected:
  /// Register individual call
  void RegisterCall(const CallConstraint& cc) {
    CheckCallFormat(cc);
    RegisterArgs(cc);
  }

  void CheckCallFormat(const CallConstraint& cc) {
    const auto& args = cc.GetArguments();
    nn_ = (int)args.size();
    n_ = (int)std::round(std::sqrt(1+8*nn_)-1)/2;
    MP_ASSERT_ALWAYS(nn_ == n_*(n_+1)/2,
                     "An SDP variable should be declared with\n"
                     "SDP_x {1..n}:\n"
                     "   kth_diag({i in 1..n, j in i..n} x[i,j]) >= eps;");
    MP_ASSERT_ALWAYS(n_ == this->MC().VarUsage(cc.GetResultVar()),
                     "Constraint\n"
                     "   kth_diag({i in 1..n, j in i..n} x[i,j]) >= eps;\n"
                     "should appear n times\n"
                     "for any n*n SDP variable");
  }

  void RegisterArgs(const CallConstraint& cc) {
    const auto& args = cc.GetArguments();
    auto sdp_index = this->MC().AddSDPVar(n_);
    int row=0, col=0;
    for (auto v: args) {  // unpack vector into L
      RegisterScalarElement(sdp_index, v, row, col);
      if (n_ == ++row)
      { row = ++col; }
    }
  }

  void RegisterScalarElement(int sdp_i, int v, int row, int col) {
    MP_ASSERT_ALWAYS(this->MC().is_free(v),
                     fmt::format(
                         "Variable _svar[{}] is bounded\nand cannot be in an SDP matrix.\n"
                         "To bound an element of an SDP matrix,\n"
                         "equate it to a (bounded) scalar variable\nor any expression",
                         v+1));
    MP_ASSERT_ALWAYS(!this->MC().HasInitExpression(v),
                     fmt::format(
                         "Element [{}, {}] of SDP matrix {}\nis an expression.\n"
                         "Only proper variables are allowed",
                         row+1, col+1, sdp_i+1));
    auto sdpvi = this->MC().GetSDPVarInfo(v);
    MP_ASSERT_ALWAYS(sdpvi.sdp_var_index_ != sdp_i,
                     fmt::format(
                         "Variable _svar[{}] is used repeatedly\nin the same SDP matrix",
                         v+1));
    MP_ASSERT_ALWAYS(sdpvi.sdp_var_index_ < 0,
                     fmt::format(
                         "Variable _svar[{}] is used in\nmultiple SDP matrices",
                         v+1));
    this->MC().SetSDPVarInfo(v, {sdp_i, row, col});
    // @todo add links for scalar<->sdp variable values
    // just before scalar->scalar
  }

private:
  /// current SDP var sizes
  int n_{}, nn_{};
};

}  // namespace mp

#endif // MP_FLAT_CVT_SDP_H
