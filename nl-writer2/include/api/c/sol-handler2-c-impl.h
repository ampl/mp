/**
 * "Implementation" of SOLHandler2_C:
 * A C++ class "wrapping" it
 * in order to interface it for NLReader2
 */
#ifndef SOLHANDLER2CIMPL_H
#define SOLHANDLER2CIMPL_H

#include "api/c/sol-handler2-c.h"

#include "mp/sol-handler2.h"

namespace mp {

/// Declare VecReader<>
template <class Value>
class VecReader;

/// Wrap SOLHandler2_C into a C++ class,
/// in order to interface it for NLReader2
class SOLHandler2_C_Impl
    : public SOLHandler2 {
public:
  /// Construct
  SOLHandler2_C_Impl(NLW2_SOLHandler2_C* psh2)
    : solh2_c_(*psh2) { }

  /** The NLHeader used to write the NL file. */
  NLHeader Header() const {
    assert(SH().Header);
    auto h_c = SH().Header(SH().p_user_data_);
    NLHeader hdr;
    *(NLProblemInfo_C*)(&hdr) = h_c.pi;
    *(NLInfo_C*)(&hdr) = h_c.nli;

    return hdr;
  }

  /** Receive solve message.
   *  The message always ends with '\n'.
   *
   *  @param nbs: number of backspaces
   *  in the original solve message.
   *  So many characters should be skipped
   *  from the message if printed straightaway.
   *  AMPL solver drivers can supply the message
   *  with initial backspaces to indicate
   *  that so many characters should be skipped
   *  when printing. For example, if the driver prints
   *  MINOS 5.51:
   *  and exits, and the message starts with that again,
   *  this part should be skipped.
   */
  void OnSolveMessage(const char* s, int nbs) {
    assert(SH().OnSolveMessage);
    SH().OnSolveMessage(SH().p_user_data_, s, nbs);
  }

  /**
   * Can be ignored by external systems.
   * @return non-zero to stop solution input.
   */
  int OnAMPLOptions(const AMPLOptions& ao) {
    assert(SH().OnAMPLOptions);

    AMPLOptions_C ao_c;
    ao_c.n_options_ = (int)ao.options_.size();
    std::copy(ao.options_.begin(), ao.options_.end(),
              ao_c.options_);
    ao_c.has_vbtol_ = ao.has_vbtol_;
    ao_c.vbtol_ = ao.vbtol_;

    return SH().OnAMPLOptions(SH().p_user_data_, ao_c);
  }

  /**
   * Dual values for algebraic constraints,
   * if provided in the solution.
   * Number of values <= NumAlgCons().
   * Implementation:
   *
   *   duals.reserve(rd.Size());
   *   while (rd.Size())
   *     duals.push_back(rd.ReadNext());
   */
  template <class VecReader>
  void OnDualSolution(VecReader& vr) {
    assert(SH().OnDualSolution);

    static_assert(
    std::is_same<VecReader, mp::VecReader<double> >::value,
        "internal error: "
        "SOLHandler2_C_Impl::OnDualSolution() requires "
        "a VecReader<double>");

    if (auto nvals = vr.Size()) {
      SH().OnDualSolution(SH().p_user_data_, nvals, &vr);
    }
  }

  /**
   * Variable values, if provided.
   * Number of values <= NumVars().
   */
  template <class VecReader>
  void OnPrimalSolution(VecReader& vr) {
    assert(SH().OnPrimalSolution);

    static_assert(
    std::is_same<VecReader, mp::VecReader<double> >::value,
        "internal error: "
        "SOLHandler2_C_Impl::OnPrimalSolution() requires "
        "a VecReader<double>");

    if (auto nvals = vr.Size()) {
      SH().OnPrimalSolution(SH().p_user_data_, nvals, &vr);
    }
  }

  /**
   * Receive notification of the objective index
   * used by the driver (solver option 'objno').
   */
  void OnObjno(int on) {
    assert(SH().OnObjno);
    SH().OnObjno(SH().p_user_data_, on);
  }

  /**
   * Receive notification of the solve code.
   */
  void OnSolveCode(int sc) {
    assert(SH().OnSolveCode);
    SH().OnSolveCode(SH().p_user_data_, sc);
  }

  /**
   * OnIntSuffix().
   *
   * For constraints, can include values for
   * logical constraints (after algebraic.)
   * Sparse representation - can be empty
   * (i.e., all values zero.)
   *
   * const auto& si = sr.SufInfo();
   * int kind = si.Kind();
   * int nmax = nitems_max[kind & 3];
   * const std::string& name = si.Name();
   * const std::string& table = si.Table();
   * while (sr.Size()) {
   *   std::pair<int, int> val = sr.ReadNext();
   *   if (val.first<0 || val.first>=nmax) {
   *     sr.SetError(mp::SOL_Read_Bad_Suffix,
   *       "bad suffix element index");
   *     return;
   *   }
   *   suf[val.first] = val.second;
   * }
   * if (mp::SOL_Read_OK == sr.ReadResult())    // Can check
   *   RegisterSuffix(kind, name, table, suf);
   */
//  template <class SuffixReader>
//  void OnIntSuffix(SuffixReader& ) { }

  /**
   * Same as OnIntSuffix(), but
   * sr.ReadNext() returns pair<int, double>
   */
//  template <class SuffixReader>
//  void OnDblSuffix(SuffixReader& ) { }

protected:
  const NLW2_SOLHandler2_C& SH() const { return solh2_c_; }

private:
  /// Just store copy
  const NLW2_SOLHandler2_C solh2_c_;
};

}  // namespace mp

#endif // SOLHANDLER2CIMPL_H
