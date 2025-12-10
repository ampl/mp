/*
 mp::NLFeederAbstract.
 Abstract-type interface for mp::NLFeeder.
 This interface is more flexible to use in a dynamically dispatched
 implementation.

 NL is a format for representing optimization problems such as linear,
 quadratic, nonlinear, complementarity and constraint programming problems
 in discrete or continuous variables. It is described in the technical report
 "Writing .nl Files" (https://ampl.github.io/nlwrite.pdf).

 Usage: recommended via mp::NLSOL class.
 See the MP2NL driver for example implementation.


 Copyright (C) 2024 AMPL Optimization, Inc.

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

#ifndef NLFEEDERABSTRACT_H
#define NLFEEDERABSTRACT_H

#include "mp/nl-feeder.h"
#include "mp/nl-feeder-abstract-defs.h"


namespace mp {

/**
    \rst
    NLFeederAbstract: writes model details on request
    via provided callback objects.
    Abstract-type interface for mp::NLFeeder.
    This interface is more flexible to use in a dynamically dispatched
    implementation.

    Implement methods listed in the 1st part, if needed.

    Usage: recommended via mp::NLSOL class.
    See the MP2NL driver for example implementation.


    For the NL format, variables and constraints must have certain order.

    **Variable ordering:**
      first continuous, then integer.
      Some solvers might require more elaborate ordering, see NLHeader.

    **Constraint ordering:**
      first algebraic (including complementarity), then logical.
    Some solvers might require nonlinear constraints first.

    @param: *Impl* is the final CRTP type
    derived from `~mp::NLFeederAbstract`.

    @param: *ExprType* is a type storing expressions from
    methods such as `~FeedExpr`. If not used,
    it can be any default-constructible type.
    \endrst
 */
template <typename Impl, typename ExprType>
class NLFeederAbstract
    : public NLFeeder<Impl, ExprType> {
public:
  /** Base class. */
  using BaseClass = NLFeeder<Impl, ExprType>;

  /** The expression type. */
  using typename BaseClass::Expr;

  ///////////////////////////////////////////////////////////////////////
  /// Methods which can/should be implemented in the derived class.
  ///////////////////////////////////////////////////////////////////////

  ///////////////////// 1. NL HEADER AND OPTIONS /////////////////
  /** Provide NLHeader.
     *
     *	This method is called first.
     *
     *  NLHeader summarizes the model and provides some
     *  technical parameters,
     *  such as text/binary NL format. */
  NLHeader Header() { assert(0); return {}; }

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
  const char* ObjDescription(int i) { return ""; }

  /** Provide type of objective \a i.
     *  - 0 - minimization;
     *  - 1 - maximization. */
  int ObjType(int i) { return {}; }

  /** Feed gradient for objective \a i.
   *  Should include entries for all potentially
   *  nonzero elements (sparsity pattern).
   *
   *  Implementation skeleton:
   *      if (obj_grad[i].size()) {
   *        auto& svw = svwf.MakeVectorWriter(obj_grad[i].size());
   *        for (size_t j=0; j<obj_grad.size(); ++j)
   *          svw.Write(obj_grad[j].var_index, obj_grad[j].coef);
   *      }
   */
  void FeedObjGradient(int i, BasicSparseDblVecWrtFactory& ) { }

  /** Feed nonlinear expression of objective \a i.
   *
   *  The default implementation below feeds constant 0
   *  (linear models.)
   *
   *  Implementation example:
   *      ew.EPut(obj_root_expr[i]);
   */
  void FeedObjExpression(int , BasicNLExprWriter<Expr>& ew)
  { ew.NPut(0.0); }


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
   *        auto& linw = dv.GetLinExprWriter();
   *        for (int i=0; i<lin_nnz; ++i)
   *          linw.Write(linexp_var[i], linexp_coef[i]);
   *        /////////// Write the expression tree:
   *        auto& ew = dv.GetExprWriter();
   *        ew.EPut(root_expr);
   *      }
     */
  void FeedDefinedVariables(int i, BasicNLDefVarWriterFactory<Expr>& ) { }


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
  void FeedVarBounds(BasicNLVarBndWriter& ) { }


  ///////////////// 5. CONSTRAINT BOUNDS & COMPLEMENTARITY ///////
  /// \rst
  /// Algebraic constraint bounds (for a single constraint):
  /// either range (lb, ub),
  /// or complementarity info (k, cvar), when k>0.
  ///
  /// For a complementarity constraint to hold, if cvar is at
  ///	its lower bound, then body >= 0; if cvar is at its upper
  /// bound, then body <= 0;
  ///	and if cvar is strictly between its bounds, then body = 0.
  /// The integer k in a complementarity constraint line indicates
  /// which bounds on cvar are finite: 1 and 3 imply a finite
  /// lower bound; 2 and 3 imply a finite upper bound; 0 (which
  ///	should not occur) would imply no finite bounds, i.e.,
  /// body = 0 must always hold.
  ///
  /// Example:
  ///
  /// .. code-block:: ampl
  ///
  ///    ampl: var x; var y; var z;
  ///	   ampl: s.t. Compl1: x+y >= 3 complements x-z <= 15;
  ///	   ampl: s.t. Compl2: -2 <= 2*y+3*z <= 13 complements 6*z-2*x;
  ///	   ampl: expand;
  ///	   subject to Compl1:
  ///					3 <= x + y
  ///			 complements
  ///					x - z <= 15;
  ///
  ///	   subject to Compl2:
  ///					-2 <= 2*y + 3*z <= 13
  ///			 complements
  ///					-2*x + 6*z;
  ///
  ///	   ampl: solexpand;
  ///	   Nonsquare complementarity system:
  ///					4 complementarities including 2 equations
  ///					5 variables
  ///	   subject to Compl1.L:
  ///					x + y + Compl1$cvar = 0;
  ///
  ///	   subject to Compl1.R:
  ///					-15 + x - z <= 0
  ///			 complements
  ///					Compl1$cvar <= -3;
  ///
  ///	   subject to Compl2.L:
  ///					2*y + 3*z - Compl2$cvar = 0;
  ///
  ///	   subject to Compl2.R:
  ///					-2*x + 6*z
  ///			 complements
  ///					-2 <= Compl2$cvar <= 13;
  ///
  /// \endrst

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
  void FeedConBounds(BasicNLConBndWriter& ) { }


  ///////////////////// 6. CONSTRAINTS /////////////////////
  /** Description of constraint \a i
   *    (\a i in 0..num_algebraic_cons+num_logical_cons-1).
   *  With WantNLComments()==true, this is
     *  written to text-format NL as a comment. */
  const char* ConDescription(int ) { return ""; }

  /** Feed the linear part of algebraic constraint \a i.
    * For smooth solvers, should contain entries for all
    * potential nonzeros (Jacobian sparsity pattern).
    *
    *  Implementation skeleton:
    *      if (con_grad[i].size()) {
    *        auto& sv = svw.MakeVectorWriter(con_grad[i].size());
    *        for (size_t j=0; j<con_grad.size(); ++j)
    *          sv.Write(con_grad[j].var_index, con_grad[j].coef);
    *      }
    */
  void FeedLinearConExpr(int i, BasicSparseDblVecWrtFactory& ) { }

  /** Feed nonlinear expression of constraint \a i.
     *  Algebraic constraints (num_algebraic_cons)
     *  come before logical (num_logical_cons).
     *  For linear constraints, the expression should be
   *  constant 0, which is implemented as default.
     */
  void FeedConExpression(int , BasicNLExprWriter<Expr>& ew)
  { ew.NPut(0.0); }


  ///////////////////// 7. EXPRESSIONS /////////////////////
  /** Feed native expression.
     *  This method is recursively called from NLWriter,
     *  when Feeder uses ExprWriter::EPut().
     *  Feeder should not call this method itself.
   */
  void FeedExpr(Expr e, BasicNLExprWriter<Expr>& ) { }


  ///////////////////// 8. PL-SOS CONSTRAINTS ////////////
  /**
     *  The below feature is for AMPL's internal
     *  linearization of piecewise-linear functions.
     *  For user-definable SOS constraints, use suffixes
     *  .sosno/.ref.
     *
     *  The below is a feeder interface
     *  for .sos/.sosref suffixes.
     *  The feeder can provide 3 sparse vectors:
     *  - .sos for variables:
     *    Each nonzero value defines SOS group number.
     *    Negative means SOS Type 2, positive - SOS Type 1.
     *  - .sos for constraints:
     *    Each nonzero value denotes a constraint used in a
     *    linearization of an SOS. The constraint can be deleted
   *    by the solver driver if using solver's SOS.
     *  - .sosref for variables:
     *    SOS weights. Variables participating in an SOS having
     *    zero weights are involved in linearization and can be
     *    deleted if the solver accepts SOS natively.
   *
   *  Implementation:
   *      auto& sosv = plsos.StartSOSVars(nvsos);
   *      for (int i=0; i<nvsos; ++i)
   *        sosv.Write(i, vsos[i]);
   *      if (ncsos) {
   *        auto& sosc = plsos.StartSOSCons(ncsos);
   *        for ....
   *      }
   *      auto& sosrefv = plsos.StartSOSREFVars(ac->nsosref);
   *      ....
    */
  void FeedPLSOS(BasicPLSOSWriter& ) { }


  /** Function definition. */
  using typename BaseClass::FuncDef;

  /** Provide definition
   *  of function \a i, i=0..num_funcs-1. */
  FuncDef Function(int i) { return {}; }


  ///////////////////// 10. RANDOM VARIABLES /////////////////////
  /// Random variables.
  /// Undocumented feature. SNL2006.
  /// @note Not adapting this one.
  using BaseClass::FeedRandomVariables;


  ///////////////////// 11. COLUMN SIZES /////////////////////

  /** Jacobian column sizes (including potential nonzeros).
     *  Should feed LP/Jacobian column sizes
     *  for all but the last variable.
     *
     *  This is called before writing Jacobian rows.
   *
   *  Implementation skeleton:
   *      if (WantColumnSizes())
   *        for (int i=0; i < num_vars+num_rand_vars-1; ++i)
   *          csw.Write(col_size[i]);
   */
  void FeedColumnSizes(BasicNLColSizeWriter& ) { }


  ///////////////////// 12. INITIAL GUESSES /////////////////////
  /** Initial primal guesses.
   *
   *  Implementation: write all meaningfuls entries (incl. zeros.)
   *      if (ini_guess.size()) {
   *        auto ig = igw.MakeVectorWriter(ini_guess.size());
   *        for (size_t i=0; i<ini_guess.size(); ++i)
   *          ig.Write(ini_guess[i].index_, ini_guess[i].value_);
   *      }
   */
  void FeedInitialGuesses(BasicSparseDblVecWrtFactory& ) { }

  /** Initial dual guesses. */
  void FeedInitialDualGuesses(BasicSparseDblVecWrtFactory& ) { }


  ///////////////////// 13. SUFFIXES /////////////////////
  /** Feed suffixes.
     *
     *  For constraints, assume ordering:
     *  first algebraic, then logical.
   *
   *  Implementation: write all non-0 entries (0 is the default.)
   *      while (....) {
   *        auto sw = swf.StartIntSuffix(  // or ...DblSuffix
   *          suf_name, kind, n_nonzeros);
   *        for (int i=0; i<n_nonzeros; ++i)
   *          sw.Write(index[i], value[i]);
   *      }
     */
  void FeedSuffixes(BasicNLSuffixWriterFactory& ) { }


  //////////////////// 14. ROW/COLUMN NAMES ETC /////////////////////
  /** FeedRowAndObjNames:
   *  Provide constraint, then objective names.
   *  Name information is optional.
   *
   *  Implementation:
   *      if ((output_desired) && wrt)
   *        for (i: ....)
   *          wrt.Write( name[i].c_str() );
     */
  void FeedRowAndObjNames(BasicNLNameFileWriter& wrt) { }

  /** Provide deleted row names.*/
  void FeedDelRowNames(BasicNLNameFileWriter& ) { }

  /** Provide variable names. */
  void FeedColNames(BasicNLNameFileWriter& ) { }

  /** Provide unused variable names. */
  void FeedUnusedVarNames(BasicNLNameFileWriter& ) { }

  /** Provide {fixed variable, extra info} pairs.
     *  This includes defined eliminated variables.
   *
   *  Implementation:
   *      if ((output_desired) && wrt)
   *        for (....)
   *          wrt.Write( name[i].c_str(), comment[i].c_str() );
     */
  void FeedFixedVarNames(BasicNLNameFileWriter& ) { }

  /** Provide {obj name, constant term} pairs.
   *
   *  Implementation:
   *      if (wrt)
   *        for (....)
   *          wrt.Write( name[i].c_str(), (double)obj_offset[i] );
     */
  void FeedObjAdj(BasicNLNameFileWriter& ) { }


  ///////////////////////////////////////////////////////////////////////
  /// Methods which normally should not be reimplemented.
  ///////////////////////////////////////////////////////////////////////

  /// Dispatcher method
  template <class ObjGradWriterFactory>
  void FeedObjGradient(int i, ObjGradWriterFactory& );

  /// Dispatcher method
  template <class ObjExprWriter>
  void FeedObjExpression(int , ObjExprWriter& ew);


  /// Dispatcher method
  template <class DefVarWriterFactory>
  void FeedDefinedVariables(int i, DefVarWriterFactory& );


  /// Dispatcher method
  template <class VarBoundsWriter>
  void FeedVarBounds(VarBoundsWriter& );


  /// Dispatcher method
  template <class ConBoundsWriter>
  void FeedConBounds(ConBoundsWriter& );


  /// Dispatcher method
  template <class ConLinearExprWriterFactory>
  void FeedLinearConExpr(int i, ConLinearExprWriterFactory& );

  /// Dispatcher method
  template <class ConExprWriter>
  void FeedConExpression(int , ConExprWriter& ew);


  /// Dispatcher method
  template <class ExprWriter>
  void FeedExpr(Expr e, ExprWriter& );


  /// Dispatcher method
  template <class PLSOSWriter>
  void FeedPLSOS(PLSOSWriter& );


  /// Dispatcher method
  template <class ColSizeWriter>
  void FeedColumnSizes(ColSizeWriter& );


  /// Dispatcher method
  template <class IGWriter>
  void FeedInitialGuesses(IGWriter& );

  /// Dispatcher method
  template <class IDGWriter>
  void FeedInitialDualGuesses(IDGWriter& );


  /// Dispatcher method
  template <class SuffixWriterFactory>
  void FeedSuffixes(SuffixWriterFactory& );


  /// Dispatcher method
  template <class RowObjNameWriter>
  void FeedRowAndObjNames(RowObjNameWriter& wrt);

  /// Dispatcher method
  template <class DelRowNameWriter>
  void FeedDelRowNames(DelRowNameWriter& );

  /// Dispatcher method
  template <class ColNameWriter>
  void FeedColNames(ColNameWriter& );

  /// Dispatcher method
  template <class UnusedVarNameWriter>
  void FeedUnusedVarNames(UnusedVarNameWriter& );

  /// Dispatcher method
  template <class FixedVarNameWriter>
  void FeedFixedVarNames(FixedVarNameWriter& );

  /// Dispatcher method
  template <class ObjOffsetWriter>
  void FeedObjAdj(ObjOffsetWriter& );

};

}  // namespace mp

#endif // NLFEEDERABSTRACT_H
