/**
 NL Writer Python bindings.

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

#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "mp/nl-solver.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// Variables' data by pointers
struct NLWPY_ColData {
  /// Num vars
  int num_col_;
  /// lower bounds
  std::vector<double> lower_;
  /// upper bounds
  std::vector<double> upper_;
  /// type: NLW2_VarType...
  /// Set to NULL if all continuous.
  std::vector<int> type_;
};

/// Sparse matrix.
struct NLWPY_SparseMatrix {
  /// Size of the start_ array:
  /// N cols (for colwise) / N rows (for rowwise),
  /// depending on format_.
  int num_colrow_;
  /// Format (NLW2_MatrixFormat...).
  /// Only rowwise supported.
  int format_;
  /// Nonzeros
  size_t num_nz_;
  /// Row / col starts
  std::vector<size_t> start_;
  /// Entry index
  std::vector<int> index_;
  /// Entry value
  std::vector<double> value_;
};

/// NLWPY_NLModel.
/// @todo check array lengths etc.
///
/// @note In contrast to C/C++, NLWPY copies all data
///   so the provided arrays/views can be deleted straightaway.
class NLWPY_NLModel {
public:
  /// Construct
  NLWPY_NLModel(std::string nm={})
    : prob_name_(std::move(nm)),
      nlme_(prob_name_.c_str())
  { }

  /// Add variables (all at once.).
  /// @todo ty can be None.
  void SetCols(int n,
               std::vector<double> lb,
               std::vector<double> ub,
               std::vector<int> ty) {
    vars_.num_col_ = n;
    vars_.lower_ = std::move(lb);
    vars_.upper_ = std::move(ub);
    vars_.type_  = std::move(ty);
    nlme_.SetCols({n,
                   vars_.lower_.data(),
                   vars_.upper_.data(),
                   vars_.type_.data()
                  });
  }

  /// Add variable names
  void SetColNames(std::vector<std::string> nm) {
    var_names_=std::move(nm);
    var_names_c_.resize(var_names_.size());
    for (auto i=var_names_.size(); i--; )
      var_names_c_[i] = var_names_[i].c_str();
    nlme_.SetColNames(var_names_c_.data());
  }

  /// Add linear constraints (all at once).
  /// Only rowwise matrix supported.
  void SetRows(
      int nr,
      std::vector<double> rlb, std::vector<double> rub,
      int format,     // TODO enum
      size_t nnz,
      std::vector<size_t> st,
      /// Entry index
      std::vector<int> ind,
      /// Entry value
      std::vector<double> val
      ) {
    num_row_=nr; row_lb_=std::move(rlb); row_ub_=std::move(rub);
    A_={
      nr, format,
      nnz, std::move(st), std::move(ind), std::move(val)
    };
    nlme_.SetRows(nr, row_lb_.data(), row_ub_.data(),
                  {
                    nr, format, nnz,
                    A_.start_.data(), A_.index_.data(),
                    A_.value_.data()
                  });
  }

  /// Add constraint names
  void SetRowNames(std::vector<std::string> nm) {
    row_names_=std::move(nm);
    row_names_c_.resize(row_names_.size());
    for (auto i=row_names_.size(); i--; )
      row_names_c_[i] = row_names_[i].c_str();
    nlme_.SetRowNames(row_names_c_.data());
  }

  /// Add linear objective (only single objective supported.)
  /// Sense: NLW2_ObjSenseM....
  /// Coefficients: dense vector.
  void SetLinearObjective(int sense, double c0,
                          std::vector<double> c) {
    obj_sense_=sense; obj_c0_=c0; obj_c_=std::move(c);
    nlme_.SetLinearObjective(sense, c0, obj_c_.data());
  }

  /// Add Q for the objective quadratic part 0.5 @ x.T @ Q @ x.
  /// Format: NLW2_HessianFormat...
  void SetHessian(int nr,
                  int format,     // TODO enum
                  size_t nnz,
                  std::vector<size_t> st,
                  /// Entry index
                  std::vector<int> ind,
                  /// Entry value
                  std::vector<double> val
                  ) {
    Q_format_ = format;
    Q_={
      nr, 0,
      nnz, std::move(st), std::move(ind), std::move(val)
    };
    nlme_.SetHessian(format, {
      nr, 0, nnz,
      Q_.start_.data(), Q_.index_.data(),
      Q_.value_.data()
    });
  }

  /// Set obj name
  void SetObjName(std::string nm) {
    obj_name_=std::move(nm);
    nlme_.SetObjName(obj_name_.c_str());
  }

  /// Get the model
  const mp::NLModel& GetModel() const { return nlme_; }

private:
  /// Store the strings/arrays to keep the memory
  std::string prob_name_ {"NLWPY_Model"};
  mp::NLModel nlme_;
  NLWPY_ColData vars_ {};
  std::vector<std::string> var_names_ {};
  std::vector<const char*> var_names_c_ {};
  NLWPY_SparseMatrix A_ {};
  int num_row_ {};
  std::vector<double> row_lb_ {};
  std::vector<double> row_ub_ {};
  std::vector<std::string> row_names_ {};
  std::vector<const char*> row_names_c_ {};
  int obj_sense_ {};
  double obj_c0_ {};
  std::vector<double> obj_c_ {};
  int Q_format_ {};
  NLWPY_SparseMatrix Q_ {};
  std::string obj_name_ {"obj[1]"};
};

mp::NLSolution NLW2_Solve(mp::NLSolver& nls,
                          const NLWPY_NLModel& mdl,
                          const std::string& solver,
                          const std::string& solver_opts) {
  return nls.Solve(mdl.GetModel(), solver, solver_opts);
}

///////////////////////////////////////////////////////////////////////////////
PYBIND11_MODULE(nlwpy, m) {
    m.doc() = R"pbdoc(
        AMPL NL Writer library Python API
        ---------------------------------

        .. currentmodule:: nlwpy

        .. autosummary::
           :toctree: _generate

           NLW2_MakeNLOptionsBasic_Default
           add
           subtract
    )pbdoc";

    /// NLOptionsBasic
    py::class_<NLW2_NLOptionsBasic_C>(m, "NLW2_NLOptionsBasic")
        .def(py::init<>())
        .def_readwrite("n_text_mode_", &NLW2_NLOptionsBasic_C::n_text_mode_)
        .def_readwrite("want_nl_comments_", &NLW2_NLOptionsBasic_C::want_nl_comments_)
        .def_readwrite("flags_", &NLW2_NLOptionsBasic_C::flags_)
        ;

    m.def("NLW2_MakeNLOptionsBasic_Default", &NLW2_MakeNLOptionsBasic_C_Default, R"pbdoc(
        Use this to create default options for NLModel.
    )pbdoc");

    /// NLModel
    py::class_<NLWPY_NLModel>(m, "NLW2_NLModel")
        .def(py::init<const char*>())
        .def("SetCols", &NLWPY_NLModel::SetCols)
        .def("SetColNames", &NLWPY_NLModel::SetColNames)
        .def("SetRows", &NLWPY_NLModel::SetRows)
        .def("SetRowNames", &NLWPY_NLModel::SetRowNames)
        .def("SetLinearObjective", &NLWPY_NLModel::SetLinearObjective)
        .def("SetHessian", &NLWPY_NLModel::SetHessian)
        .def("SetObjName", &NLWPY_NLModel::SetObjName)
        ;

    /// NLSolution
    py::class_<mp::NLSolution>(m, "NLW2_NLSolution")
        .def_readwrite("solve_result_", &mp::NLSolution::solve_result_)
        .def_readwrite("nbs_", &mp::NLSolution::nbs_)
        .def_readwrite("solve_message_", &mp::NLSolution::solve_message_)
        .def_readwrite("obj_val_", &mp::NLSolution::obj_val_)
        .def_readwrite("x_", &mp::NLSolution::x_)
        .def_readwrite("y_", &mp::NLSolution::y_)
        ;

    /// NLSolver
    py::class_<mp::NLSolver>(m, "NLW2_NLSolver")
        .def(py::init<>())
        .def("SetFileStub", &mp::NLSolver::SetFileStub)
        .def("SetNLOptions", &mp::NLSolver::SetNLOptions)
        .def("GetErrorMessage", &mp::NLSolver::GetErrorMessage)
        .def("Solve", &NLW2_Solve)
        ;

    // -------------------------------------------------------------------
#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
