/*
 NL Solver, part of implementation.

 Copyright (C) 2024 AMPL Optimization Inc.

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
#include <random>
#include <filesystem>

#if defined(_WIN32) || defined(_WIN64)
  #include <io.h>     // _mktemp[_s]
#elif __APPLE__
  #include <unistd.h> // mkdtemp
#else
  #include <cstdlib>  // mkdtemp
#endif

#include "mp/nl-solver.hpp"
#include "mp/nl-opcodes.h"

extern "C"
NLW2_NLOptionsBasic_C NLW2_MakeNLOptionsBasic_C_Default() {
  NLW2_NLOptionsBasic_C result;
  result.n_text_mode_ = 0;
  result.want_nl_comments_ = 0;
  result.flags_ = 1;      // want out suffixes

  return result;
}

namespace mp {

/// Specialize NLFeeder for NLModel
class NLFeeder_Easy
    : public NLFeeder<NLFeeder_Easy, void*> {
public:
  /// Construct
  NLFeeder_Easy(const NLModel& nls, NLW2_NLOptionsBasic_C opts)
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
          svw.Write(VPerm(j), c ? c[j] : 0.0);
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
      auto pos_end = Q.num_nz_;
      for (auto i=NLME().NumCols(); i--; ) {
        for (auto pos=Q.start_[i]; pos!=pos_end; ++pos) {
          auto coef = 0.5 * Q.value_[pos];
          auto prod1 = sumw.OPut2(nl::MUL);
          prod1.NPut(coef);
          auto prod2 = prod1.OPut2(nl::MUL);
          prod2.VPut(VPerm(i), NLME().ColName(i));   // x
          prod2.VPut(VPerm(Q.index_[pos]),
                     NLME().ColName(Q.index_[pos])); // y
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
      vbw.WriteLbUb(vars.lower_[VPermInv(i)],
                    vars.upper_[VPermInv(i)]);
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
    assert(!A.num_nz_ || NLW2_MatrixFormatRowwise == A.format_);
    assert(NLME().NumRows() == A.num_colrow_);
    auto start = A.start_[i];
    auto end = (i < A.num_colrow_-1) ? A.start_[i+1] : A.num_nz_;
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
        csw.Write(col_sizes_[VPermInv(i)]);
  }


  ///////////////////// 12. INITIAL GUESSES /////////////////////
  /** Initial primal guesses.
   *
   *  Implementation:
   *      if (ini_guess.size()) {
   *        auto ig = igw.MakeVectorWriter(ini_guess.size());
   *        for (size_t i=0; i<ini_guess.size(); ++i)
   *          ig.Write(ini_guess[i].index_, ini_guess[i].value_);
   *      }
   */
  template <class IGWriter>
  void FeedInitialGuesses(IGWriter& igw) {
    auto ini = NLME().Warmstart();
    if (ini.num_) {
      auto ig = igw.MakeVectorWriter(ini.num_);
      for (int i=0; i<ini.num_; ++i)
        ig.Write(VPerm( ini.index_[i] ),
                 ini.value_[i]);
    }
  }

  /** Initial dual guesses. */
  template <class IDGWriter>
  void FeedInitialDualGuesses(IDGWriter& igw) {
    auto ini = NLME().DualWarmstart();
    if (ini.num_) {
      auto ig = igw.MakeVectorWriter(ini.num_);
      for (int i=0; i<ini.num_; ++i)
        ig.Write(ini.index_[i], ini.value_[i]);
    }
  }


  ///////////////////// 13. SUFFIXES /////////////////////
  /** Feed suffixes.
     *
     *  For constraints, assume ordering:
     *  first algebraic, then logical.
   *
   *  Implementation: write all non-0 entries
   *      while (....) {
   *        auto sw = swf.StartIntSuffix(  // or ...DblSuffix
   *          suf_name, kind, n_nonzeros);
   *        for (int i=0; i<n_nonzeros; ++i)
   *          sw.Write(index[i], value[i]);
   *      }
     */
  template <class SuffixWriterFactory>
  void FeedSuffixes(SuffixWriterFactory& swf) {
    for (const auto& suf: NLME().Suffixes()) {
      int nnz=0;
      for (auto v: suf.values_)
        if (v)
          ++nnz;
      bool ifVars = (0==(suf.kind_&3));
      if (suf.kind_ & 4) {
        auto sw = swf.StartDblSuffix(
              suf.name_.c_str(), suf.kind_, nnz);
        for (size_t i=0; i<suf.values_.size(); ++i)
          if (suf.values_[i])
            sw.Write(ifVars ? VPerm(i) : i,
                     suf.values_[i]);
      } else {
        auto sw = swf.StartIntSuffix(
              suf.name_.c_str(), suf.kind_, nnz);
        for (size_t i=0; i<suf.values_.size(); ++i)
          if (suf.values_[i])
            sw.Write(ifVars ? VPerm(i) : i,
                     std::round(suf.values_[i]));
      }
    }
  }


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
        wrt << NLME().ColNames()[VPermInv(i)];
    }
  }

  void ExportPreproData(NLModel::PreprocessData &pd) {
    pd.vperm_.resize(NLME().NumCols());
    pd.vperm_inv_.resize(NLME().NumCols());
    for (auto i=pd.vperm_.size(); i--; ) {
      pd.vperm_[i] = VPerm(i);
      pd.vperm_inv_[i] = VPermInv(i);
    }
  }


protected:
  const NLModel& NLME() const { return nlme_; }
  NLModel& NLME() { return nlme_; }

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
        nlv_obj_[Q.index_[i]] = true;
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
  /// Variables info element
  using VarInfo = std::pair<int, int>;

  /// Var direct permutation
  int VPerm(int i) const
  { assert(i<(int)var_perm_.size()); return var_perm_[i].first; }
  /// Var inverse permutation
  int VPermInv(int i) const
  { assert(i<(int)var_perm_.size()); return var_perm_[i].second; }

  void FillObjNonzeros() {
    obj_grad_supp_.resize(NLME().NumCols());
    // Linear part -- if provided
    if (NLME().ObjCoefficients()) {
      for (auto i=NLME().NumCols(); i--; )
        obj_grad_supp_[i] = (NLME().ObjCoefficients()[i]);
    }
    // QP part
    auto Q = NLME().Hessian();
    if (Q.num_nz_) {
      assert(Q.num_colrow_ == NLME().NumCols());
      auto pos_end = Q.num_nz_;
      for (auto i=NLME().NumCols(); i--; ) {
        for (auto pos=Q.start_[i]; pos!=pos_end; ++pos) {
          obj_grad_supp_[i] = true;             // x
          obj_grad_supp_[Q.index_[pos]] = true; // y
        }
        pos_end = Q.start_[i];
      }
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
  NLModel nlme_;
  NLW2_NLOptionsBasic_C nlopt_;

  std::vector<bool> nlv_obj_;      // if var nonlinear in obj
  std::vector<VarInfo> var_perm_;

  std::vector<bool> obj_grad_supp_;  // obj sparsity pattern
  std::vector<int> col_sizes_;

  NLHeader header_;
};

std::string NLModel::WriteNL(
    const std::string &fln, NLW2_NLOptionsBasic_C opts,
    NLUtils &ut, PreprocessData &pd) {
  NLFeeder_Easy nlf(*this, opts);
  nlf.ExportPreproData(pd);
  return WriteNLFile(fln, nlf, ut).second;
}

double NLModel::ComputeObjValue(const double *x) const {
  double result {obj_c0_};
  for (auto i=NumCols(); i--; )
    result += obj_c_[i] * x[i];
  if (Q_.num_nz_) {
    auto pos_end = Q_.num_nz_;
    for (auto i=NumCols(); i--; ) {
      for (auto pos=Q_.start_[i]; pos!=pos_end; ++pos) {
        result
            += 0.5 * Q_.value_[pos] * x[i] * x[Q_.index_[pos]];
      }
      pos_end = Q_.start_[i];
    }
  }
  return result;
}

NLSolver::NLSolver()
  : p_ut_(&utils_) { }

NLSolver::NLSolver(mp::NLUtils* put)
  : p_ut_(put ? put : &utils_) { }

NLSolver::~NLSolver() { DestroyAutoStub(); }

void NLSolver::InitAutoStub() {
  // init file stub
  std::random_device dev;
  std::mt19937 prng(dev());
  std::uniform_int_distribution<unsigned long> rand(0);
  auto path = std::filesystem::temp_directory_path();
  path /= "nlw2_"; // via '/'
  char rnds[64] = "rndhex";
  std::snprintf(rnds, sizeof(rnds)-1, "%lX", rand(prng));
  path += rnds;    // no '/'

  path += "_XXXXXX";
  pathstr_ = path.string();

#if defined(_WIN32) || defined(_WIN64)
  auto p1 = _mktemp((char*)pathstr_.c_str());
  assert(p1);
  if (!std::filesystem::create_directory(pathstr_))
    Utils().myexit("Could not create temp dir '"
                   + pathstr_ + "'");
#else
  if (!mkdtemp((char*)pathstr_.c_str()))
    Utils().myexit("Could not create a temp dir\n"
                   "from pattern '" + pathstr_ + "'");
#endif
  path = pathstr_;
  // Plus filename
  std::snprintf(rnds, sizeof(rnds)-1, "%lX", rand(prng));
  path /= rnds;
  filestub_ = path.string();
}

void NLSolver::DestroyAutoStub() {
  // delete temp folder if created
  if (pathstr_.size()) {
    std::error_code ec;
    std::filesystem::remove_all(pathstr_, ec);
    if (ec)
      Utils().log_warning("Failed to remove temp dir '%s': %s",
                          pathstr_.c_str(), ec.message().c_str());
  }
}

void NLSolver::SetFileStub(std::string stub) {
  if (stub.size()) {
    filestub_ = stub;
    filestubCustom_ = true;
  }
}


bool NLSolver::LoadModel(const NLModel& mdl) {
  NLFeeder_Easy nlf(mdl, nl_opts_);
  nlf.ExportPreproData(pd_);
  p_nlheader_.reset(new NLHeader(nlf.Header()));
  return LoadModel(nlf);
}

bool NLSolver::Solve(const std::string& solver,
                  const std::string& solver_opts) {
  if (GetFileStub().empty())
    return (err_msg_="NLSolver: provide filestub.", false);
  if (solver.empty())
    return (err_msg_="NLSolver: provide solver.", false);
  auto call = solver
      + ' ' + GetFileStub()
      + " -AMPL "
      + solver_opts;
  if (auto status = std::system(call.c_str()))
    return (err_msg_="NLSolver: call \""
        + call + "\" failed (code "
        + std::to_string(status) + ").", false);
  return true;
}


/// Specialize SOLHandler for NLModel
class SOLHandler_Easy
    : public SOLHandler {
public:
  /// Construct
  SOLHandler_Easy(const NLHeader& h,
                   const NLModel::PreprocessData& pd,
                   NLSolution& sol)
    : header_(h), pd_(pd), sol_(sol) { }

  /** The NLHeader used to write the NL file. */
  NLHeader Header() const { return header_; }

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
    sol_.nbs_ = nbs;
    sol_.solve_message_ = s;
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
  void OnDualSolution(VecReader& rd) {
    sol_.y_.clear();
    sol_.y_.reserve(header_.num_algebraic_cons
                    + header_.num_logical_cons);
    while (rd.Size())         // No permutation
      sol_.y_.push_back(rd.ReadNext());
  }

  /**
   * Variable values, if provided.
   * Number of values <= NumVars().
   */
  template <class VecReader>
  void OnPrimalSolution(VecReader& rd) {
    sol_.x_.clear();
    sol_.x_.resize(header_.num_vars);
    for (int i=0; rd.Size(); ++i)    // Permute
      sol_.x_[pd_.vperm_inv_[i]] = rd.ReadNext();
  }

  /**
   * Receive notification of the solve code.
   * Solve result codes docu:
   * https://mp.ampl.com/features-guide.html#solve-result-codes
   */
  void OnSolveCode(int c) { sol_.solve_result_=c; }

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
  void OnIntSuffix(SuffixReader& sr)
  { OnSuffix(sr); }

  /**
   * Same as OnIntSuffix(), but
   * sr.ReadNext() returns pair<int, double>
   */
  template <class SuffixReader>
  void OnDblSuffix(SuffixReader& sr)
  { OnSuffix(sr); }

protected:
  int NItemsMax(int kind) const {
    switch (kind & 3) {    // Without the Float bit
    case 0: return header_.num_vars;
    case 1: return
          header_.num_algebraic_cons
          + header_.num_logical_cons;
    case 2: return header_.num_objs;
    default: return 1;
    }
  }
  template <class SuffixReader>
  void OnSuffix(SuffixReader& sr) {
    const auto& si = sr.SufInfo();
    int kind = si.Kind();
    int nmax = NItemsMax(kind);
    std::vector<double> values(nmax);
    const std::string& name = si.Name();
    const std::string& table = si.Table();
    bool ifVars = (0==(kind & 3));
    while (sr.Size()) {
      auto val = sr.ReadNext();
      if (val.first<0 || val.first>=nmax) {
        sr.SetError(NLW2_SOLRead_Bad_Suffix,
                    "bad suffix element index");
        return;
      }
      values[
          ifVars ? pd_.vperm_inv_[ val.first ] : val.first]
                                                 = val.second;
    }
    if (NLW2_SOLRead_OK == sr.ReadResult())    // Can check
      sol_.suffixes_.Add({name, table, kind,
                          std::move(values)});
  }


private:
  NLHeader header_;
  const NLModel::PreprocessData& pd_;
  NLSolution& sol_;
};

NLSolution NLSolver::ReadSolution() {
  NLSolution result;
  if (!p_nlheader_)
    return (err_msg_
            ="NLSolver: "
             "can only ReadSolution(void) after loading NLModel",
            result);
  SOLHandler_Easy solh(*p_nlheader_, pd_, result);
  ReadSolution(solh);
  return result;
}

}  // namespace mp
