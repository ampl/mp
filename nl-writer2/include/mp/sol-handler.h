/*
 mp::SOLHandler.
 Interface for a solution handler to be used with mp::SOLReader2.
 mp::SOLReader2 is a zero-overhead reader
 implemented with inline template code.

 AMPL SOL is a format for representing solutions of optimization problems
 such as linear, quadratic, nonlinear, complementarity
 and constraint programming problems in discrete or continuous variables,
 used by AMPL. It is not documented yet.

 Usage: recommended via mp::NLSOL class.

 If you need to read solutions separately, proceed as follows:

   MySOLHandler handler;
   mp::NLUtils utils;
   auto status = mp::ReadSOLFile(
          "filestub", handler, utils);
   if (mp::SOL_Read_OK != status.first) {
     printf(status.second.c_str());
   }

 where handler is an object that receives information on solution
 components. Below is an interface for such a handler.

 Copyright (C) 2023 AMPL Optimization, Inc.

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc. disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Gleb Belov
 */

#ifndef SOLHandler_H
#define SOLHandler_H

#include <vector>
#include <cassert>

#include "mp/nl-header.h"

namespace mp {

/**
  \rst
    SOLHandler: reads solution details on request
    via provided callback objects.
    See the examples folder.

    `~mp::SOLHandler` can be used as a base class for other handlers,
    or just be an interface prototype.
  \endrst
 */
class SOLHandler {
public:
	/** The NLHeader used to write the NL file. */
	NLHeader Header() const { return {}; }

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
  void OnSolveMessage(const char* s, int nbs) { }

  /**
   * AMPL internal options.
   */
  struct AMPLOptions {
    /// Option values
    std::vector<long> options_;
    /// Has vbtol?
    bool has_vbtol_;
    /// Vbtol
    double vbtol_;
  };

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
    while (rd.Size())
      rd.ReadNext();       // read & forget by default
  }

  /**
   * Variable values, if provided.
   * Number of values <= NumVars().
   */
  template <class VecReader>
  void OnPrimalSolution(VecReader& rd) {
    while (rd.Size())
      rd.ReadNext();       // read & forget by default
  }

  /**
   * Receive notification of the objective index
   * used by the driver (solver option 'objno'-1).
   */
  void OnObjno(int ) { }

  /**
   * Receive notification of the solve code.
   * Solve result codes docu:
   * https://mp.ampl.com/features-guide.html#solve-result-codes
   */
  void OnSolveCode(int ) { }

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
   *     sr.SetError(NLW2_SOLRead_Bad_Suffix,
   *       "bad suffix element index");
   *     return;
   *   }
   *   suf[val.first] = val.second;
   * }
   * if (NLW2_SOLRead_OK == sr.ReadResult())    // Can check
   *   RegisterSuffix(kind, name, table, suf);
   */
  template <class SuffixReader>
  void OnIntSuffix(SuffixReader& sr) {
    while (sr.Size())
      sr.ReadNext();       // read & forget by default
  }

  /**
   * Same as OnIntSuffix(), but
   * sr.ReadNext() returns pair<int, double>
   */
  template <class SuffixReader>
  void OnDblSuffix(SuffixReader& sr) {
    while (sr.Size())
      sr.ReadNext();       // read & forget by default
  }

};

}  // namespace mp

#endif // SOLHandler_H
