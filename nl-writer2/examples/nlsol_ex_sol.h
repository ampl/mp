/**
 * SOLHanlder2 implementation for the nlsol_ex example.
 *
 */

#ifndef NLSOL_EX_SOL_H
#define NLSOL_EX_SOL_H

#include <string>
#include <cstring>

#include "mp/sol-handler2.h"
#include "nlsol_ex_mdl.h"

class ExampleSOLHandler2
    : public mp::SOLHandler2 {
public:
  /// Construct
  ExampleSOLHandler2(ExampleModel& mdl)
    : mdl_(mdl) { }

	/** The NLHeader used to write the NL file. */
	mp::NLHeader Header() const { return Model().Header(); }

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
   *  this part should be skipped just now.
   */
  void OnSolveMessage(const char* s, int nbs) {
    auto n = std::strlen(s);
    if (nbs) {
      if (nbs > n)
        nbs = n;
      s += nbs;
      n -= nbs;
    }
    if (n > 0) {
      std::printf("%s", s);
      std::fflush(stdout);
    }
  }

  /**
   * Can be ignored by external systems.
   * @return non-zero to stop solution input.
   */
  int OnAMPLOptions(const AMPLOptions& ) { return 0; }

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
  void OnDualSolution(VecReader& rd) {
    Model().sol_y_.reserve(rd.Size());
    for (int i=0; rd.Size(); ++i)
      Model().sol_y_.push_back( rd.ReadNext() );
  }

  /**
   * Variable values, if provided.
   * Number of values <= NumVars().
   */
  template <class VecReader>
  void OnPrimalSolution(VecReader& rd) {
    Model().sol_x_.reserve(rd.Size());
    for (int i=0; rd.Size(); ++i)
      Model().sol_x_.push_back( rd.ReadNext() );
  }

  /**
   * Receive notification of the objective index
   * used by the driver (solver option 'objno').
   */
  void OnObjno(int on) { Model().objno_=on; }

  /**
   * Receive notification of the solve code.
   */
  void OnSolveCode(int sr) { Model().solve_result_ = sr; }

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
  template <class SuffixReader>
  void OnIntSuffix(SuffixReader& sr) {
    const auto& si = sr.SufInfo();
    int kind = si.Kind();
    const std::string& name = si.Name();
    assert(0 == (kind & 4));
    auto& suf_val = Model().suf_out_[std::make_pair(name, kind)];
    while (sr.Size()) {
      std::pair<int, int> val = sr.ReadNext();
      if (suf_val.size() <= val.first)
        suf_val.resize(val.first+1);
      suf_val[val.first] = val.second;
    }
  }

  /**
   * Same as OnIntSuffix(), but
   * sr.ReadNext() returns pair<int, double>
   */
  template <class SuffixReader>
  void OnDblSuffix(SuffixReader& sr) {
    const auto& si = sr.SufInfo();
    int kind = si.Kind();
    const std::string& name = si.Name();
    assert(0 != (kind & 4));
    auto& suf_val = Model().suf_out_[std::make_pair(name, kind)];
    while (sr.Size()) {
      std::pair<int, double> val = sr.ReadNext();
      if (suf_val.size() <= val.first)
        suf_val.resize(val.first+1);
      suf_val[val.first] = val.second;
    }
  }


  //////////////////////////////////////////////
  /// Print solution.
  /// @param stub: stub filename used.
  void PrintSolution(const std::string& stub);


protected:
  const ExampleModel& Model() const { return mdl_; }
  ExampleModel& Model() { return mdl_; }

private:
  ExampleModel& mdl_;
};


void ExampleSOLHandler2::PrintSolution(
    const std::string& stub) {
  printf(
        "\n     ********** SOLUTION (%s.sol) ***********\n",
        stub.c_str());
  printf("%s\n", "Duals:");
  for (auto y: Model().sol_y_)
    printf("%.17g\t", y);
  printf("\n%s\n", "Primals:");
  for (auto x: Model().sol_x_)
    printf("%.17g\t", x);

  printf("\nObjno used: %d, solve_result_num: %d\n",
         Model().objno_, Model().solve_result_);

  printf("\n%s\n", "Suffixes:");
  for (auto suf: Model().suf_out_) {
    printf("  %s  [kind=%d]: ",
           suf.first.first.c_str(), suf.first.second);
    for (auto v: suf.second)
      printf("\t%.17g", v);
    printf("\n");
  }
  printf("%s\n", "");
}

#endif // NLSOL_EX_SOL_H
