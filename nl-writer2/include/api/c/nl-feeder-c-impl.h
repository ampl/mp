/*
 C API: "Implementation" of NLFeeder_C:
 A C++ class "wrapping" it
 in order to interface it for NLWriter2.

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
#ifndef NLFeederCIMPL_H
#define NLFeederCIMPL_H

#include <functional>

#include "api/c/nl-feeder-c.h"   // C wrapper

#include "mp/nl-feeder.h"        // C++ base


namespace mp {

/// Implementation:
/// Wrap NLW2_NLFeeder_C into a C++ class,
/// in order to interface it for NLWriter2
class NLW2_NLFeeder_C_Impl
    : public NLFeeder<NLW2_NLFeeder_C_Impl, void*> {
public:
  /// typedef base class
  using Base = NLFeeder<NLW2_NLFeeder_C_Impl, void*>;

  /// Construct
  NLW2_NLFeeder_C_Impl(NLW2_NLFeeder_C* pnlf2)
    : nlf2_c_(*pnlf2) { }

  ///////////////////// 1. NL HEADER AND OPTIONS /////////////////
  /** Provide NLHeader.
   *
   *	This method is called first.
   *
   *  NLHeader summarizes the model and provides some
   *  technical parameters,
   *  such as text/binary NL format. */
  NLHeader Header() {
    assert(NLF().Header);
    auto h_c = NLF().Header(NLF().p_user_data_);
    NLHeader hdr;
    *(NLProblemInfo_C*)(&hdr) = h_c.pi;
    *(NLInfo_C*)(&hdr) = h_c.nli;

    return hdr;
  }

  /// NL comments?
  bool WantNLComments() const { return false; }

  /// The maximum number of significant digits written.
  /// The default value requests full precision, which
  /// might be the shortest representation that, when
  /// converted to binary and properly rounded, will
  /// give exactly the binary value stored in the computer.
  int OutputPrecision() const { return 0; }

  /// Write bounds first?
  /// The default is yes in AMPL, controlled by
  /// (the value of option nl_permute) & 32
  /// (the bit is 0 for yes).
  /// Changing this option is deprecated, see
  /// https://netlib.org/ampl/changes.
  bool WantBoundsFirst() const { return true; }

  /// Want Jacobian column sizes?
  /// Required by some nonlinear solvers.
  /// Options: 0 - none, 1 - cumulative,
  /// 2 - non-cumulative.
  /// This option controls how ColSizeWriter
  /// writes the provided sizes (which should be
  /// non-cumulative).
  int WantColumnSizes() const { return 1; }


  ///////////////////// 2. OBJECTIVES /////////////////////
  /** Description for objective function \a i
   *    (\a i in 0..num_objs-1).
   *  With WantNLComments()==true, this is
   *  written to text-format NL as a comment. */
  const char* ObjDescription(int i) {
    assert(NLF().ObjDescription);
    return NLF().ObjDescription(NLF().p_user_data_, i);
  }

  /** Provide type of objective \a i.
   *  0 - minimization;
   *  1 - maximization. */
  int ObjType(int i) {
    assert(NLF().ObjType);
    return NLF().ObjType(NLF().p_user_data_, i);
  }

  /** Feed gradient for objective \a i.
   *  Should include entries for all potentially
   *  nonzero elements (sparsity pattern).
   *
   *  Implementation skeleton:
   *      if (obj_grad[i].size()) {
   *        auto sv = svwf.MakeVectorWriter(obj_grad[i].size());
   *        for (size_t j=0; j<obj_grad.size(); ++j)
   *          sv.Write(obj_grad[j].var_index, obj_grad[j].coef);
   *      }
   */
  template <class ObjGradWriterFactory>
  void FeedObjGradient(int i, ObjGradWriterFactory& svwf) {
    assert(NLF().ObjGradientNNZ);
    int nnz = NLF().ObjGradientNNZ(NLF().p_user_data_, i);
    if (nnz) {
      auto sv = svwf.MakeVectorWriter(nnz);
      assert(NLF().FeedObjGradient);
      // We need the callback to be callable from C.
      // Due to the NLWriter2 having different types
      // for text vs binary, static functions won't work.
      std::function<void(int, double)> svw
          = [&sv](int i, double v){
        sv.Write(i, v);
      };
      NLF().FeedObjGradient(NLF().p_user_data_, i, &svw);
    }
  }

  /** Feed nonlinear expression of objective \a i.
   *
   *  Implementation example:
   *      ew.EPut(obj_root_expr[i]);
   *
   *  Details of ObjExprWriter: see NLWriter2. */
  template <class ObjExprWriter>
  void FeedObjExpression(int i, ObjExprWriter& ew) {
    Base::FeedObjExpression(i, ew);      // reuse default
  }


  ///////////////////// 3. DEFINED VARIABLES /////////////////////
  /** Defined variables.
   *
   *  Classical NL writes first the defined variables
   *  which are used in several places (constraints and/or
   *  objectives). Defined variables used in a single place
   *  (1 constraint, or 1 objective), are written
   *  just before the expression tree of their usage.
   *
   *  For most solvers, this requirement can be ignored
   *  and this method can return all defined variables
     *  in the first group (for \a i=0).
   *
   *	The method is guaranteed to be called in the following order:
   *		1. For \a i=0;
   *		2. For \a i>0, increasing, before constraint \a (i-1)'s expression;
   *		3. For \a i<0, decreasing, before objective \a (-i-1)'s expression.
   *
   *  @param i:
   *		- For \a i=0, feed a sequence of defined variables
   *			used in several constraints and/or objectives.
   *		- For \a i>0, feed the defined variables used solely
   *			in constraint \a i-1.
   *		- For \a i<0, feed the defined variables used solely
   *			in objective \a -i-1.
   *
   *  Implementation skeleton:
   *      // dvar_index in num_vars..num_vars+num_defvars-1.
   *      for (int dvar_index: dvar_indexes[i]) {
   *        auto dv = dvw.StartDefVar(dvar_index, lin_nnz, name_or_comment);
   *        /////////// Write the linear part:
   *        auto linw = dv.GetLinExprWriter();
   *        for (int i=0; i<lin_nnz; ++i)
   *          linw.Write(linexp_var[i], linexp_coef[i]);
   *        /////////// Write the expression tree:
   *        auto ew = dv.GetExprWriter();
   *        ew.EPut(root_expr);
   *      }
   */
  template <class DefVarWriterFactory>
  void FeedDefinedVariables(int i, DefVarWriterFactory& ) {
    // None yet
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
    assert(NLF().FeedVarBounds);
    std::function<void(double, double)> vbwc
        = [&vbw](double i, double v){
      vbw.WriteLbUb(i, v);
    };
    NLF().FeedVarBounds(NLF().p_user_data_, &vbwc);
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
  void FeedConBounds(ConBoundsWriter& crw) {
    assert(NLF().FeedConBounds);
    std::function<void(NLW2_AlgConRange_C* )> crw_c
        = [&crw](NLW2_AlgConRange_C* bnd_c){
      AlgConRange bnd;
      bnd.k = bnd_c->k;
      bnd.cvar = bnd_c->cvar;
      bnd.L = bnd_c->L;
      bnd.U = bnd_c->U;
      crw.WriteAlgConRange(bnd);
    };
    NLF().FeedConBounds(NLF().p_user_data_, &crw_c);
  }


  ///////////////////// 6. CONSTRAINTS /////////////////////
  /** Description of constraint \a i
   *    (\a i in 0..num_algebraic_cons+num_logical_cons-1).
   *  With WantNLComments()==true, this is
   *  written to text-format NL as a comment. */
  const char* ConDescription(int i) {
    assert(NLF().ConDescription);
    return NLF().ConDescription(NLF().p_user_data_, i);
  }

  /** Feed the linear part of algebraic constraint \a i.
    * For smooth solvers, should contain entries for all
    * potential nonzeros (Jacobian sparsity pattern).
    *
    *  Implementation skeleton:
    *      if (con_grad[i].size()) {
    *        auto sv = svw.MakeVectorWriter(con_grad[i].size());
    *        for (size_t j=0; j<con_grad.size(); ++j)
    *          sv.Write(con_grad[j].var_index, con_grad[j].coef);
    *      }
    */
  template <class ConLinearExprWriterFactory>
  void FeedLinearConExpr(int i, ConLinearExprWriterFactory& clewf) {
    assert(NLF().LinearConExprNNZ);
    int nnz = NLF().LinearConExprNNZ(NLF().p_user_data_, i);
    if (nnz) {
      auto sv = clewf.MakeVectorWriter(nnz);
      assert(NLF().FeedLinearConExpr);
      std::function<void(int, double)> svw
          = [&sv](int i, double v){
        sv.Write(i, v);
      };
      NLF().FeedLinearConExpr(NLF().p_user_data_, i, &svw);
    }
  }

  /** Feed nonlinear expression of constraint \a i.
   *  Algebraic constraints (num_algebraic_cons)
   *  come before logical (num_logical_cons).
   *  For linear constraints, the expression should be
   *  constant 0.
   */
  template <class ConExprWriter>
  void FeedConExpression(int i, ConExprWriter& ew) {
    Base::FeedConExpression(i, ew);
  }


  ///////////////////// 7. EXPRESSIONS /////////////////////
  template <class ExprWriter>
  void FeedExpr(Expr , ExprWriter& ) { }


  ///////////////////// 8. PL-SOS CONSTRAINTS ////////////
  template <class PLSOSWriter>
  void FeedPLSOS(PLSOSWriter& ) { }


  ///////////////////// 9. FUNCTIONS /////////////////////
  FuncDef Function(int ) { return {}; }


  ///////////////////// 10. RANDOM VARIABLES /////////////////////
  template <class RandVarWriterFactory>
  void FeedRandomVariables(RandVarWriterFactory& ) { }


  ///////////////////// 11. COLUMN SIZES /////////////////////

  /** Jacobian column sizes.
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
    assert(NLF().FeedColumnSizes);
    std::function<void(int )> f = [&csw](int cs) {
      csw.Write(cs);
    };
    NLF().FeedColumnSizes(NLF().p_user_data_, &f);
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
  void FeedInitialGuesses(IGWriter& );

  /** Initial dual guesses. */
  template <class IDGWriter>
  void FeedInitialDualGuesses(IDGWriter& );


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
  void FeedSuffixes(SuffixWriterFactory& );


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
    assert(NLF().FeedRowAndObjNames);
    if (NLF().want_row_and_obj_names_ && wrt) {
      std::function<void(const char*)> wrt_c
          = [&wrt](const char* name) {
        wrt << name;
      };
      NLF().FeedRowAndObjNames(NLF().p_user_data_, &wrt_c);
    }
  }

  /** Provide deleted row names.*/
  template <class DelRowNameWriter>
  void FeedDelRowNames(DelRowNameWriter& wrt) {
    assert(NLF().FeedDelRowNames);
    if (NLF().want_del_row_names_ && wrt) {
      std::function<void(const char*)> wrt_c
          = [&wrt](const char* name) {
        wrt << name;
      };
      NLF().FeedDelRowNames(NLF().p_user_data_, &wrt_c);
    }
  }

  /** Provide variable names. */
  template <class ColNameWriter>
  void FeedColNames(ColNameWriter& wrt) {
    assert(NLF().FeedColNames);
    if (NLF().want_col_names_ && wrt) {
      std::function<void(const char*)> wrt_c
          = [&wrt](const char* name) {
        wrt << name;
      };
      NLF().FeedColNames(NLF().p_user_data_, &wrt_c);
    }
  }

  /** Provide unused variable names. */
  template <class UnusedVarNameWriter>
  void FeedUnusedVarNames(UnusedVarNameWriter& wrt) {
    assert(NLF().FeedUnusedVarNames);
    if (NLF().want_unused_var_names_ && wrt) {
      std::function<void(const char*)> wrt_c
          = [&wrt](const char* name) {
        wrt << name;
      };
      NLF().FeedUnusedVarNames(NLF().p_user_data_, &wrt_c);
    }
  }

  /** Provide {fixed variable, extra info} pairs.
   *  This includes defined eliminated variables.
   *
   *  Implementation:
   *      if ((output_desired) && wrt)
   *        for (....)
   *          wrt << typename Writer::StrStrValue
   *          { name[i].c_str(), comment[i].c_str() };
   */
  template <class FixedVarNameWriter>
  void FeedFixedVarNames(FixedVarNameWriter& wrt) {
    assert(NLF().FeedFixedVarNames);
    if (NLF().want_fixed_var_names_ && wrt) {
      std::function<void(const char*, const char*)> wrt_c
          = [&wrt](const char* name, const char* comment) {
        wrt << typename FixedVarNameWriter::StrStrValue
        { name, comment };
      };
      NLF().FeedFixedVarNames(NLF().p_user_data_, &wrt_c);
    }
  }

  /** Provide {obj name, constant term} pairs.
   *
   *  Implementation:
   *      if (wrt)
   *        for (....)
   *          wrt << typename Writer::StrDblValue
   *          { name[i].c_str(), (double)obj_offset[i] };
   */
  template <class ObjOffsetWriter>
  void FeedObjAdj(ObjOffsetWriter& wrt) {
    assert(NLF().FeedObjAdj);
    if (NLF().want_obj_adj_ && wrt) {
      std::function<void(const char*, double)> wrt_c
          = [&wrt](const char* name, double num) {
        wrt << typename ObjOffsetWriter::StrDblValue
        { name, num };
      };
      NLF().FeedObjAdj(NLF().p_user_data_, &wrt_c);
    }
  }


protected:
  const NLW2_NLFeeder_C& NLF() const { return nlf2_c_; }


private:
  /// Just store copy
  const NLW2_NLFeeder_C nlf2_c_;
};

}  // namespace mp

#endif // NLFeederCIMPL_H
