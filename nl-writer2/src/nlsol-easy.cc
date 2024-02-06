/**
 NL Solver "Easy", for special model classes

 Copyright (C) 2023 AMPL Optimization Inc.

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Gleb Belov
 */

#include <cstring>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>

#include "mp/nlsol-easy.h"
#include "mp/nlsol.h"
#include "mp/nl-opcodes.h"

extern "C"
NLW2_NLOptionsBasic NLW2_MakeNLOptionsBasic_Default() {
  NLW2_NLOptionsBasic result;
  result.n_text_mode_ = 0;
  result.want_nl_comments_ = 0;
  result.flags_ = 1;      // want out suffixes

  return result;
}

namespace mp {

/// Specialize NLFeeder2 for NLSOL_Easy
class NLFeeder2_Easy
    : public NLFeeder2<NLFeeder2_Easy, void*> {
public:
  /// Construct
  NLFeeder2_Easy(NLModel_Easy& nls, NLW2_NLOptionsBasic opts)
    : nlme_(nls), nlopt_(opts) { Init(); }

  /// NL header
  NLHeader Header() const { return header_; }

  /// NL comments?
  bool WantNLComments() const { return nlopt_.want_nl_comments_; }

  ///////////////////// 2. OBJECTIVES /////////////////////
  const char* ObjDescription(int i)
  { assert(0==i); return NLME().ObjName(); }

  /** Provide type of objective \a i.
     *  0 - minimization;
     *  1 - maximization. */
  int ObjType(int i)
  { assert(0==i); return NLME().ObjSense(); }

  /** Feed gradient for objective \a i.
   *  Should include entries for all potentially
   *  nonzero elements (sparsity pattern).
   *
   *  Implementation skeleton:
   *      if (obj_grad[i].size()) {
   *        auto svw = svwf.MakeVectorWriter(obj_grad[i].size());
   *        for (size_t j=0; j<obj_grad.size(); ++j)
   *          svw.Write(obj_grad[j].var_index, obj_grad[j].coef);
   *      }
   */
  template <class ObjGradWriterFactory>
  void FeedObjGradient(int , ObjGradWriterFactory& svwf) {
    if (header_.num_obj_nonzeros) {
      auto svw = svwf.MakeVectorWriter(header_.num_obj_nonzeros);
      auto c = NLME().ObjCoefficients();
      for (int j=0; j<header_.num_vars; ++j)
        if (obj_grad_supp_[j])
          svw.Write(VPerm(j), c[j]);
    }
  }

  /** Feed nonlinear expression of objective \a i.
   *
   *  The default implementation below feeds constant 0
   *  (linear models.)
   *
   *  Implementation example:
   *      ew.EPut(obj_root_expr[i]);
   *
   *  Details of ObjExprWriter: see NLWriter2. */
  template <class ObjExprWriter>
  void FeedObjExpression(int , ObjExprWriter& ew) {
    auto Q = NLME().Hessian();
    auto c0 = NLME().ObjOffset();
    if (!Q.num_nz_) {
      ew.NPut(c0);             // linear, write the offset
    } else {
      bool if_offset = (c0);
      int num_el = Q.num_nz_ + if_offset;
      auto sumw = ew.OPutN(nl::SUM, num_el);
      if (if_offset)
        sumw.NPut(c0);
      int pos_end = Q.num_nz_;
      for (auto i=NLME().NumCols(); i--; ) {
        for (auto pos=Q.start_[i]; pos!=pos_end; ++pos) {
          auto coef = 0.5 * Q.value_[pos];
          auto prod1 = sumw.OPut2(nl::MUL);
          prod1.NPut(coef);
          auto prod2 = prod1.OPut2(nl::MUL);
          prod2.VPut(VPerm(i));             // x
          prod2.VPut(VPerm(Q.index_[pos])); // y
        }
        pos_end = Q.start_[i];
      }
    }
  }


  ///////////////////// 4. VARIABLE BOUNDS /////////////////////
  /** Bounds for variables (except defined variables).
   *  Use +-inf for missing lower and/or upper bounds.
     *  Note that variable type is given by variable ordering,
   *  see NLHeader.
   *
   *  Implementation skeleton:
   *      for (int i = 0; i < hdr.num_vars; i++)
   *        vbw.WriteLbUb(lb[i], ub[i]);
   */
  template <class VarBoundsWriter>
  void FeedVarBounds(VarBoundsWriter& vbw) {
    auto vars = NLME().ColData();
    for (int i = 0; i < header_.num_vars; i++)
      vbw.WriteLbUb(VPermRev(vars.lower_[i]),
                    VPermRev(vars.upper_[i]));
  }


  ///////////////// 5. CONSTRAINT BOUNDS & COMPLEMENTARITY ///////

  /** Bounds/complementarity for all algebraic constraints
   *  (\a num_algebraic_cons).
   *
   *  Implementation skeleton:
   *      for (int j=0; j<hdr.num_algebraic_cons; j++) {
   *        AlgConRange bnd;
   *        if (compl_var && compl_var[j]) {
   *          j = compl_var[j]-1;
   *          bnd.k = 0;
   *          if (vlb[j] > negInfinity)
   *            bnd.k = 1;
   *          if (vub[j] < Infinity)
   *            bnd.k |= 2;
   *          assert(bnd.k);
   *          bnd.cvar = j;
   *        } else {
   *          bnd.L = clb[j];
   *          bnd.U = cub[j];
   *        }
   *        cbw.WriteAlgConRange(bnd);
   *      }
   */
  template <class ConBoundsWriter>
  void FeedConBounds(ConBoundsWriter& cbw) {
    for (int j=0; j<header_.num_algebraic_cons; j++) {
      AlgConRange bnd;
      bnd.L = NLME().RowLowerBounds()[j];
      bnd.U = NLME().RowUpperBounds()[j];
      cbw.WriteAlgConRange(bnd);
    }
  }


  ///////////////////// 6. CONSTRAINTS /////////////////////
  /** Description of constraint \a i
   *    (\a i in 0..num_algebraic_cons+num_logical_cons-1).
   *  With WantNLComments()==true, this is
     *  written to text-format NL as a comment. */
  const char* ConDescription(int i) {
    return (NLME().RowNames() && NLME().RowNames()[i])
        ? NLME().RowNames()[i] : "";
  }

  /** Feed the linear part of algebraic constraint \a i.
    * For smooth solvers, should contain entries for all
    * potential nonzeros (Jacobian sparsity pattern).
    *
    * For QP, only linear constraints.
    *
    *  Implementation skeleton:
    *      if (con_grad[i].size()) {
    *        auto sv = svw.MakeVectorWriter(con_grad[i].size());
    *        for (size_t j=0; j<con_grad.size(); ++j)
    *          sv.Write(con_grad[j].var_index, con_grad[j].coef);
    *      }
    */
  template <class ConLinearExprWriterFactory>
  void FeedLinearConExpr(int i, ConLinearExprWriterFactory& svw) {
    auto A = NLME().GetA();
    assert(!A.num_nz_ || NLW2_MatrixFormatColwise == A.format_);
    assert(NLME().NumRows() == A.num_row_);
    auto start = A.start_[i];
    int end = (i < A.num_row_-1) ? A.start_[i+1] : A.num_nz_;
    if (start!=end) {
      auto sv = svw.MakeVectorWriter(end-start);
      for (auto pos=start; pos!=end; ++pos)
        sv.Write(VPerm(A.index_[pos]), A.value_[pos]);
    }
  }


  ///////////////////// 11. COLUMN SIZES /////////////////////

  /** Jacobian column sizes (including potential nonzeros).
     *  Should feed LP column sizes
     *  for all but the last variable.
   *
   *  Implementation skeleton:
   *      if (WantColumnSizes())
   *        for (int i=0; i < num_vars+num_rand_vars-1; ++i)
   *          csw.Write(col_size[i]);
   */
  template <class ColSizeWriter>
  void FeedColumnSizes(ColSizeWriter& csw) {
    if (WantColumnSizes())
      for (int i=0; i < header_.num_vars-1; ++i)
        csw.Write(col_sizes_[VPermRev(i)]);
  }


  ///////////////////// 12. INITIAL GUESSES /////////////////////
  /** Initial primal guesses.
   *
   *  Implementation:
   *      if (ini_guess.size()) {
   *        auto ig = igw.MakeVectorWriter(ini_guess.size());
   *        for (size_t i=0; i<ini_guess.size(); ++i)
   *          ig.Write(i, ini_guess[i]);
   *      }
   */
  template <class IGWriter>
  void FeedInitialGuesses(IGWriter& ) { }

  /** Initial dual guesses. */
  template <class IDGWriter>
  void FeedInitialDualGuesses(IDGWriter& ) { }


  ///////////////////// 13. SUFFIXES /////////////////////
  /** Feed suffixes.
     *
     *  For constraints, assume ordering:
     *  first algebraic, then logical.
   *
   *  Implementation:
   *      while (....) {
   *        auto sw = swf.StartIntSuffix(  // or ...DblSuffix
   *          suf_name, kind, n_nonzeros);
   *        for (int i=0; i<n_nonzeros; ++i)
   *          sw.Write(index[i], value[i]);
   *      }
     */
  template <class SuffixWriterFactory>
  void FeedSuffixes(SuffixWriterFactory& ) { }


  //////////////////// 14. ROW/COLUMN NAMES ETC /////////////////////
  /** FeedRowAndObjNames:
   *  Provide constraint, then objective names.
   *  Name information is optional.
   *
   *  Implementation:
   *      if ((output_desired) && wrt)
   *        for (i: ....)
   *          wrt << name[i].c_str();
     */
  template <class RowObjNameWriter>
  void FeedRowAndObjNames(RowObjNameWriter& wrt) {
    if (NLME().RowNames() && wrt) {
      for (int i=0; i<NLME().NumRows(); ++i)
        wrt << NLME().RowNames()[i];
      wrt << NLME().ObjName();
    }
  }

  /** Provide variable names. */
  template <class ColNameWriter>
  void FeedColNames(ColNameWriter& wrt) {
    if (NLME().ColNames() && wrt) {
      for (int i=0; i<NLME().NumCols(); ++i)
        wrt << NLME().ColNames()[i];
    }
  }


protected:
  const NLModel_Easy& NLME() const { return nlme_; }
  NLModel_Easy& NLME() { return nlme_; }

  /// Reorder variables, compute statistics
  void Init() {
    FillNonlinearVars();
    PermuteVars();           // use permutation from here
    FillObjNonzeros();
    FillColSizes();

    FillHeader();            // Some info already filled
  }

  /// Fill which vars are nonlinear in obj
  void FillNonlinearVars() {
    nlv_obj_.resize(NLME().NumCols());
    auto Q=NLME().Hessian();
    if (Q.num_nz_) {
      ++header_.num_nl_objs;                // STATS
      for (auto i=Q.num_nz_; i--; ) {
        assert(i<nlv_obj_.size());
        nlv_obj_[i] = true;
        ++header_.num_nl_vars_in_objs;      // STATS
      }
    }
  }

  /// Permute variables.
  /// For QP it should be:
  /// - nonlinear continuous, int;
  /// - linear cont, bin, int.
  void PermuteVars() {
    auto vars = NLME().ColData();
    var_perm_.resize(NLME().NumCols());
    assert(nlv_obj_.size());               // was filled
    for (auto i=var_perm_.size(); i--; ) {
      var_perm_[i] = {-2*(int)nlv_obj_[i], i};  // Prioritize nonlinear
      if (vars.type_ && vars.type_[i]) {        // int var
        ++var_perm_[i].first;
        if (nlv_obj_[i])
          ++header_.num_nl_integer_vars_in_objs; // STATS
        else {
          if (0.0!=std::fabs(vars.lower_[i]) || 1.0!=vars.upper_[i]) {
            ++var_perm_[i].first;           // For linear, first binary
            ++header_.num_linear_integer_vars;   // STATS
          } else
            ++header_.num_linear_binary_vars;    // STATS
        }
      }
    }
    std::stable_sort(var_perm_.begin(), var_perm_.end());
    // Create reverse mapping
    for (auto i=var_perm_.size(); i--; )
      var_perm_[var_perm_[i].second].first = i;
  }
  using VarInfo = std::pair<int, int>;

  /// Direct permutation
  int VPerm(int i) const
  { assert(i<(int)var_perm_.size()); return var_perm_[i].second; }
  /// Reverse permutation
  int VPermRev(int i) const
  { assert(i<(int)var_perm_.size()); return var_perm_[i].first; }

  void FillObjNonzeros() {
    obj_grad_supp_.resize(NLME().NumCols());
    // Linear part
    for (auto i=NLME().NumCols(); i--; )
      obj_grad_supp_[i] = (NLME().ObjCoefficients()[i]);
    // QP part
    auto Q = NLME().Hessian();
    assert(Q.num_row_ == NLME().NumCols());
    int pos_end = Q.num_nz_;
    for (auto i=NLME().NumCols(); i--; ) {
      for (auto pos=Q.start_[i]; pos!=pos_end; ++pos) {
        obj_grad_supp_[i] = true;             // x
        obj_grad_supp_[Q.index_[pos]] = true; // y
      }
      pos_end = Q.start_[i];
    }
    header_.num_obj_nonzeros
        = std::accumulate(obj_grad_supp_.begin(),
                          obj_grad_supp_.end(), 0.0);
  }

  void FillColSizes() {
    auto A = NLME().GetA();
    col_sizes_.resize(NLME().NumCols());
    for (auto pos=A.num_nz_; pos--; )
      ++col_sizes_[A.index_[pos]];
  }

  void FillHeader() {
    header_.format = nlopt_.n_text_mode_
        ? NL_FORMAT_TEXT : NL_FORMAT_BINARY;
    header_.flags = nlopt_.flags_;

    header_.max_var_name_len = 0;
    if (auto pn = NLME().ColNames())
      for (auto i=NLME().NumCols(); i--; )
        if (pn[i]
            && header_.max_var_name_len<(int)std::strlen(pn[i]))
          header_.max_var_name_len=(int)std::strlen(pn[i]);
    header_.max_con_name_len = 0;
    if (auto pn = NLME().RowNames())
      for (auto i=NLME().NumRows(); i--; )
        if (pn[i]
            && header_.max_con_name_len<(int)std::strlen(pn[i]))
          header_.max_con_name_len=(int)std::strlen(pn[i]);
    if (header_.max_con_name_len<(int)std::strlen(NLME().ObjName()))
      header_.max_con_name_len=(int)std::strlen(NLME().ObjName());

    header_.num_algebraic_cons = NLME().NumRows();
    // Need number of Jacobian nonzeros.
    // For linear constraints, it's just those in A.
    header_.num_con_nonzeros = NLME().GetA().num_nz_;
    // (Non)linear binary/int vars,
    // num_obj_nonzeros filled before
    // ...

    header_.num_objs = 1;
    header_.num_ranges = NLME().NumRows();
    header_.num_vars = NLME().NumCols();
    header_.prob_name = NLME().ProbName();
  }

private:
  NLModel_Easy& nlme_;
  NLW2_NLOptionsBasic nlopt_;

  std::vector<bool> nlv_obj_;      // if var nonlinear in obj
  std::vector<VarInfo> var_perm_;

  std::vector<bool> obj_grad_supp_;  // obj sparsity pattern
  std::vector<int> col_sizes_;

  NLHeader header_;
};

std::string NLModel_Easy::WriteNL(
    const std::string &fln, NLW2_NLOptionsBasic opts, NLUtils &ut) {
  NLFeeder2_Easy nlf(*this, opts);
  return WriteNLFile(fln, nlf, ut).second;
}

}  // namespace mp
