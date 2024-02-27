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
#include <pybind11/stl.h>

#include "mp/nl-solver.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// Variables' data
struct NLWPY_ColData
{
  /// lower bounds
  std::vector<double> lower_;
  /// upper bounds
  std::vector<double> upper_;
  /// type: NLW2_VarType...
  /// Set to NULL if all continuous.
  std::vector<int> type_;
};

/// Sparse matrix.
///
/// Size of the start_ array:
/// N cols (for colwise) / N rows (for rowwise),
/// depending on format_.
struct NLWPY_SparseMatrix
{
  /// Format (NLW2_MatrixFormat...).
  /// Only rowwise supported.
  NLW2_MatrixFormat format_;
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
class NLWPY_NLModel
{
public:
  /// Construct
  NLWPY_NLModel(std::string nm = {})
      : prob_name_(std::move(nm)),
        nlme_(prob_name_.c_str())
  {
  }

  /// Add variables (all at once.).
  /// @todo \a ty can be None.
  void SetCols(std::vector<double> lb,
               std::vector<double> ub,
               std::vector<int> ty)
  {
    assert(lb.size() == ub.size());
    assert(ty.empty() || ty.size() == lb.size());
    num_col_ = lb.size();
    vars_.lower_ = std::move(lb);
    vars_.upper_ = std::move(ub);
    vars_.type_ = std::move(ty);
    nlme_.SetCols({num_col_,
                   vars_.lower_.data(),
                   vars_.upper_.data(),
                   vars_.type_.data()});
  }

  /// Add variable names
  void SetColNames(std::vector<std::string> nm)
  {
    assert(nm.size() == num_col_);
    var_names_ = std::move(nm);
    var_names_c_.resize(var_names_.size());
    for (auto i = var_names_.size(); i--;)
      var_names_c_[i] = var_names_[i].c_str();
    nlme_.SetColNames(var_names_c_.data());
  }

  /// Add linear constraints (all at once).
  /// Only rowwise matrix supported.
  void SetRows(
      std::vector<double> rlb, std::vector<double> rub,
      NLW2_MatrixFormat format,
      std::vector<size_t> st,
      /// Entry index
      std::vector<int> ind,
      /// Entry value
      std::vector<double> val)
  {
    assert(rlb.size() == rub.size());
    num_row_ = rlb.size();
    row_lb_ = std::move(rlb);
    row_ub_ = std::move(rub);
    A_ = {
        format,
        std::move(st), std::move(ind), std::move(val)};
    nlme_.SetRows(num_row_,
                  row_lb_.data(), row_ub_.data(),
                  {num_row_, format, A_.value_.size(),
                   A_.start_.data(), A_.index_.data(),
                   A_.value_.data()});
  }

  /// Add constraint names
  void SetRowNames(std::vector<std::string> nm)
  {
    assert(nm.size() == num_row_);
    row_names_ = std::move(nm);
    row_names_c_.resize(row_names_.size());
    for (auto i = row_names_.size(); i--;)
      row_names_c_[i] = row_names_[i].c_str();
    nlme_.SetRowNames(row_names_c_.data());
  }

  /// Add linear objective (only single objective supported.)
  /// Sense: NLW2_ObjSenseM....
  /// Coefficients: dense vector.
  void SetLinearObjective(NLW2_ObjSense sense, double c0,
                          std::vector<double> c)
  {
    assert(c.size() == num_col_);
    obj_sense_ = sense;
    obj_c0_ = c0;
    obj_c_ = std::move(c);
    nlme_.SetLinearObjective(sense, c0, obj_c_.data());
  }

  /// Add Q for the objective quadratic part 0.5 @ x.T @ Q @ x.
  /// Format: NLW2_HessianFormat...
  void SetHessian(NLW2_HessianFormat format,
                  std::vector<size_t> st,
                  /// Entry index
                  std::vector<int> ind,
                  /// Entry value
                  std::vector<double> val)
  {
    assert(st.size() - 1 == num_col_);
    Q_format_ = format;
    Q_ = {
        NLW2_MatrixFormatIrrelevant,
        std::move(st), std::move(ind), std::move(val)};
    nlme_.SetHessian(
        format, {(int)Q_.start_.size() - 1, // don't need the last element
                 NLW2_MatrixFormatIrrelevant,
                 Q_.value_.size(),
                 Q_.start_.data(), Q_.index_.data(),
                 Q_.value_.data()});
  }

  /// Set obj name
  void SetObjName(std::string nm)
  {
    obj_name_ = std::move(nm);
    nlme_.SetObjName(obj_name_.c_str());
  }

  /// Set initial solution.
  void SetWarmstart(
      std::vector<int> i, std::vector<double> v)
  {
    assert(i.size() == v.size());
    ini_x_i_ = std::move(i);
    ini_x_v_ = std::move(v);
    nlme_.SetWarmstart(
        {(int)ini_x_i_.size(),
         ini_x_i_.data(), ini_x_v_.data()});
  }

  /// Set dual initial solution.
  void SetDualWarmstart(
      std::vector<int> i, std::vector<double> v)
  {
    assert(i.size() == v.size());
    ini_y_i_ = std::move(i);
    ini_y_v_ = std::move(v);
    nlme_.SetDualWarmstart(
        {(int)ini_y_i_.size(),
         ini_y_i_.data(), ini_y_v_.data()});
  }

  /// Add suffix.
  /// @return true iff new suffix added (vs replaced.)
  /// @note SOS constraints can be modeled as suffixes
  ///   for some AMPL solvers.
  bool AddSuffix(mp::NLSuffix suf)
  {
    return nlme_.AddSuffix(std::move(suf));
  }

  /// Get the model
  const mp::NLModel &GetModel() const { return nlme_; }

private:
  /// Store the strings/arrays to keep the memory
  std::string prob_name_{"NLWPY_Model"};
  mp::NLModel nlme_;
  int num_col_{};
  NLWPY_ColData vars_{};
  std::vector<std::string> var_names_{};
  std::vector<const char *> var_names_c_{};
  NLWPY_SparseMatrix A_{};
  int num_row_{};
  std::vector<double> row_lb_{};
  std::vector<double> row_ub_{};
  std::vector<std::string> row_names_{};
  std::vector<const char *> row_names_c_{};
  int obj_sense_{};
  double obj_c0_{};
  std::vector<double> obj_c_{};
  NLW2_HessianFormat Q_format_{};
  NLWPY_SparseMatrix Q_{};
  std::string obj_name_{"obj[1]"};

  std::vector<int> ini_x_i_;
  std::vector<double> ini_x_v_;
  std::vector<int> ini_y_i_;
  std::vector<double> ini_y_v_;
};

mp::NLSolution NLW2_Solve(mp::NLSolver &nls,
                          const NLWPY_NLModel &mdl,
                          const std::string &solver,
                          const std::string &solver_opts)
{
  return nls.Solve(mdl.GetModel(), solver, solver_opts);
}

///////////////////////////////////////////////////////////////////////////////
PYBIND11_MODULE(_nlwpy, m)
{
  m.doc() = R"pbdoc(
        AMPL NL Writer library Python API
        ---------------------------------

        .. currentmodule:: nlwpy

        .. autosummary::
           :toctree: _generate

           ObjSense
VarType
MatrixFormat
HessianFormat

NLOptionsBasic
           MakeNLOptionsBasic_Default
NLSuffix
NLModel
NLSolution
NLSolver
    )pbdoc";

  py::enum_<NLW2_ObjSense>(m, "ObjSense", py::arithmetic())
      .value("Minimize", NLW2_ObjSenseMinimize)
      .value("Maximize", NLW2_ObjSenseMaximize);
  // .export_values();     -- Leave them scoped

  py::enum_<NLW2_VarType>(m, "VarType", py::arithmetic())
      .value("Continuous", NLW2_VarTypeContinuous)
      .value("Integer", NLW2_VarTypeInteger);

  py::enum_<NLW2_MatrixFormat>(m, "MatrixFormat", py::arithmetic())
      .value("Rowwise", NLW2_MatrixFormatRowwise);

  py::enum_<NLW2_HessianFormat>(m, "HessianFormat", py::arithmetic())
      .value("Triangular", NLW2_HessianFormatTriangular)
      .value("Square", NLW2_HessianFormatSquare);

  /// NLOptionsBasic
  py::class_<NLW2_NLOptionsBasic_C>(m, "NLOptionsBasic")
      .def(py::init<>())
      .def_readwrite("n_text_mode_", &NLW2_NLOptionsBasic_C::n_text_mode_)
      .def_readwrite("want_nl_comments_", &NLW2_NLOptionsBasic_C::want_nl_comments_)
      .def_readwrite("flags_", &NLW2_NLOptionsBasic_C::flags_);

  m.def("MakeNLOptionsBasic_Default", &NLW2_MakeNLOptionsBasic_C_Default, R"pbdoc(
        Use this to create default options for NLModel.
    )pbdoc");

  /// NLSuffix
  py::class_<mp::NLSuffix>(m, "NLSuffix")
      .def(py::init<std::string, int, std::vector<double>>())
      .def(py::init<std::string, std::string, int, std::vector<double>>())
      .def_readwrite("name_", &mp::NLSuffix::name_)
      .def_readwrite("table_", &mp::NLSuffix::table_)
      .def_readwrite("kind_", &mp::NLSuffix::kind_)
      .def_readwrite("values_", &mp::NLSuffix::values_);

  /// NLSuffixSet
  py::class_<mp::NLSuffixSet>(m, "NLSuffixSet")
      .def("Find", // Find(): return None if not found
           [=](mp::NLSuffixSet const &ss,
               std::string const &name, int kind) -> py::object
           {
             auto pelem = ss.Find(name, kind);
             return py::cast(pelem);
           })
      .def("__len__", // &mp::NLSuffixSet::size - does not work)
           [](const mp::NLSuffixSet &ss)
           { return ss.size(); })
      .def(
          "__iter__", [](const mp::NLSuffixSet &ss)
          { return py::make_iterator(ss.begin(), ss.end()); },
          py::keep_alive<0, 1>()) /* Keep vector alive while iterator is used */
      .def("empty",
           [](const mp::NLSuffixSet &ss)
           { return ss.empty(); })
      .def("clear", [](mp::NLSuffixSet &ss)
           { ss.clear(); });

  /// NLModel
  py::class_<NLWPY_NLModel>(m, "NLModel")
      .def(py::init<const char *>())
      .def("SetCols", &NLWPY_NLModel::SetCols)
      .def("SetColNames", &NLWPY_NLModel::SetColNames)
      .def("SetRows", &NLWPY_NLModel::SetRows)
      .def("SetRowNames", &NLWPY_NLModel::SetRowNames)
      .def("SetLinearObjective", &NLWPY_NLModel::SetLinearObjective)
      .def("SetHessian", &NLWPY_NLModel::SetHessian)
      .def("SetObjName", &NLWPY_NLModel::SetObjName)
      .def("SetWarmstart", &NLWPY_NLModel::SetWarmstart)
      .def("SetDualWarmstart", &NLWPY_NLModel::SetDualWarmstart)
      .def("AddSuffix", &NLWPY_NLModel::AddSuffix);

  /// NLSolution
  py::class_<mp::NLSolution>(m, "NLSolution")
      .def_readwrite("solve_result_", &mp::NLSolution::solve_result_)
      .def_readwrite("nbs_", &mp::NLSolution::nbs_)
      .def_readwrite("solve_message_", &mp::NLSolution::solve_message_)
      .def_readwrite("obj_val_", &mp::NLSolution::obj_val_)
      .def_readwrite("x_", &mp::NLSolution::x_)
      .def_readwrite("y_", &mp::NLSolution::y_)
      .def_readwrite("suffixes_", &mp::NLSolution::suffixes_);

  /// NLSolver
  py::class_<mp::NLSolver>(m, "NLSolver")
      .def(py::init<>())
      .def("SetFileStub", &mp::NLSolver::SetFileStub)
      .def("SetNLOptions", &mp::NLSolver::SetNLOptions)
      .def("GetErrorMessage", &mp::NLSolver::GetErrorMessage)
      .def("Solve", &NLW2_Solve);

  // -------------------------------------------------------------------
#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
