/**
 * NLFeeder implementation for the nlsol_ex example.
 *
 */

#ifndef NLSOL_EX_NL_H
#define NLSOL_EX_NL_H

#include "mp/nl-feeder.h"
#include "nlsol_ex_mdl.h"

class ExampleNLFeeder
    : public mp::NLFeeder<ExampleModel, void*> {
public:
  /// Construct
  ExampleNLFeeder(const ExampleModel& mdl, bool binary)
    : mdl_(mdl), binary_(binary) { }

  /** Provide NLHeader.
   *
   *	This method is called first.
   *
   *  NLHeader summarizes the model and provides some
   *  technical parameters,
   *  such as text/binary NL format. */
  mp::NLHeader Header() {
		auto hdr = Model().Header();

    hdr.format = binary_
        ? mp::NLHeader::BINARY : mp::NLHeader::TEXT;

    return hdr;
  }

  /// NL comments?
  bool WantNLComments() const { return true; }

  const char* ObjDescription(int i) { return "Obj_t2t3"; }

  int ObjType(int ) { return Model().obj_sense; }

  template <class ObjGradWriter>
  void FeedObjGradient(int i, ObjGradWriter& gw)
  { WriteSparseVec(gw, Model().obj_linpart); }

  template <class ObjExprWriter>
  void FeedObjExpression(int i, ObjExprWriter& ew)
  { Model().WriteObjExpr(ew); }


  template <class DefVarWriterFactory>
  void FeedDefinedVariables(int i, DefVarWriterFactory& dvw) {
    for (int dvar0: Model().DefVarsInItem(i)) {
      int lin_nnz = CountNNZ(Model().dvar_linpart[dvar0]);
      auto dv = dvw.StartDefVar(
            dvar0 + Model().n_var,  // index after direct vars
            lin_nnz, Model().dvar_name[dvar0]);
      /////////// Write the linear part:
      auto linw = dv.GetLinExprWriter();
      WriteDense2Sparse_Pure(linw,
                             Model().dvar_linpart[dvar0]);
      /////////// Write the expression tree:
      auto ew = dv.GetExprWriter();
      Model().WriteDVarExpr(dvar0, ew);
    }
  }

  template <class VarBoundsWriter>
  void FeedVarBounds(VarBoundsWriter& vbw) {
    for (int i = 0; i < Model().n_var; i++)
      vbw.WriteLbUb(Model().var_lb[i], Model().var_ub[i]);
  }

  template <class ConBoundsWriter>
  void FeedConBounds(ConBoundsWriter& cbw) {
    for (int j=0; j<Model().n_con; j++) {
      AlgConRange bnd;
      bnd.L = Model().con_lb[j];
      bnd.U = Model().con_ub[j];
      cbw.WriteAlgConRange(bnd);
    }
  }

  const char* ConDescription(int i)
  { return Model().con_name[i]; }

  template <class ConLinearExprWriter>
  void FeedLinearConExpr(int i, ConLinearExprWriter& clw)
  { WriteSparseVec(clw, Model().con_linpart[i]); }

  template <class ConExprWriter>
  void FeedConExpression(int i, ConExprWriter& ew)
  { Model().WriteConExpr(i, ew); }

  template <class ColSizeWriter>
  void FeedColumnSizes(ColSizeWriter& csw) {
    if (WantColumnSizes())
      for (int i=0; i < Model().n_var-1; ++i)
        csw.Write(Model().col_sizes[i]);
  }

  template <class IGWriter>
  void FeedInitialGuesses(IGWriter& igw)
  { WriteDense2Sparse(igw, Model().ini_x); }

 template <class SuffixWriterFactory>
 void FeedSuffixes(SuffixWriterFactory& swf) {
   for (int i=0; i<Model().n_suf; ++i) {
     if (Model().suf[i].kind_ & 4) {     // double's
       auto sw = swf.StartDblSuffix(
             Model().suf[i].name_,
             Model().suf[i].kind_,
             CountNNZ(Model().suf[i].values_));
       WriteDense2Sparse_Pure(sw, Model().suf[i].values_);
     } else {                            // int's
       auto sw = swf.StartIntSuffix(
             Model().suf[i].name_,
             Model().suf[i].kind_,
             CountNNZ(Model().suf[i].values_));
       WriteDense2Sparse_Pure(sw, Model().suf[i].values_);
     }
   }
 }

 template <class RowObjNameWriter>
 void FeedRowAndObjNames(RowObjNameWriter& wrt) {
   if (wrt) {     // && output_desired
     for (int i=0; i<Model().n_con; ++i)
       wrt << Model().con_name[i];
     wrt << Model().obj_name;
   }
 }

 /** Provide variable names. */
 template <class ColNameWriter>
 void FeedColNames(ColNameWriter& wrt) {
   if (wrt) {     // && output_desired
     for (int i=0; i<Model().n_var; ++i)
       wrt << Model().var_name[i];
   }
 }


protected:
  const ExampleModel& Model() const { return mdl_; }

  /// Write sparse vector.
  template <class GradWriterFactory, class SparseVec>
  void WriteSparseVec(GradWriterFactory& gw, const SparseVec& vec) {
    auto gvw = gw.MakeVectorWriter(vec.size());
    for (auto v: vec)
      gvw.Write(v.first, v.second);
  }

  /// Write dense vector as sparse,
  /// using provided writer factory.
  template <class GradWriterFactory, class DenseVec>
  void WriteDense2Sparse(GradWriterFactory& gw, const DenseVec& vec) {
    auto gvw = gw.MakeVectorWriter(CountNNZ(vec));
    WriteDense2Sparse_Pure(gvw, vec);
  }

  /// Count nnz
  template <class DenseVec>
  int CountNNZ(const DenseVec& vec) {
    int n=0;
    for (auto it=std::begin(vec);
         std::end(vec) != it; ++it)
      n += (*it != 0);
    return n;
  }

  /// Write dense vector as sparse,
  /// given pure vector writer.
  template <class GradWriter, class DenseVec>
  void WriteDense2Sparse_Pure(GradWriter& gw, const DenseVec& vec) {
    for (int i=0; i<vec.size(); ++i)
      if (vec[i])
        gw.Write(i, vec[i]);
  }

  /// Write dense vector as sparse,
  /// given pure vector writer.
  template <class GradWriter, size_t n>
  void WriteDense2Sparse_Pure(GradWriter& gw, const double (&vec) [n]) {
    int sz = sizeof(vec) / sizeof(vec[0]);
    for (int i = 0; i < sz; ++i)
      if (vec[i])
        gw.Write(i, vec[i]);
  }


private:
  const ExampleModel& mdl_;
  bool binary_ {true};
};

#endif // NLSOL_EX_NL_H
