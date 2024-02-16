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
  py::array_t<const double> lower_;
  /// upper bounds
  py::array_t<const double> upper_;
  /// type: NLW2_VarType...
  /// Set to NULL if all continuous.
  py::array_t<const int> type_;
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
  py::array_t<const size_t> start_;
  /// Entry index
  py::array_t<const int> index_;
  /// Entry value
  py::array_t<const double> value_;
};

/// NLWPY_NLModel.
/// TODO check array lengths etc.
class NLWPY_NLModel {
public:
  /// Construct
  NLWPY_NLModel(const char* nm=nullptr)
    : prob_name_(nm ? nm : "NLWPY_Model"),
      nlme_(prob_name_)
  { }

  /// Add variables (all at once.).
  /// TODO ty can be None.
  void SetCols(int n,
               py::array_t<const double> lb,
               py::array_t<const double> ub,
               py::array_t<const int> ty) {
    vars_.num_col_ = n;
    vars_.lower_ = lb;
    vars_.upper_ = ub;
    vars_.type_ = ty;
    nlme_.SetCols({n,
                   vars_.lower_.data(),
                   vars_.upper_.data(),
                   vars_.type_.data()
                  });
  }

  /// Add variable names
  void SetColNames(std::vector<const char *> nm) {
    var_names_=std::move(nm);
    nlme_.SetColNames(var_names_.data());
  }

  /// Add linear constraints (all at once).
  /// Only rowwise matrix supported.
  void SetRows(
      int nr,
      py::array_t<const double> rlb, py::array_t<const double> rub,
      int format,     // TODO enum
      size_t nnz,
      py::array_t<const size_t> st,
      /// Entry index
      py::array_t<const int> ind,
      /// Entry value
      py::array_t<const double> val
      ) {
    num_row_=nr; row_lb_=rlb; row_ub_=rub;
    A_={
      nr, format,
      nnz, st, ind, val
    };
    nlme_.SetRows(nr, row_lb_.data(), row_ub_.data(),
                  {
                    nr, format, nnz,
                    A_.start_.data(), A_.index_.data(),
                    A_.value_.data()
                  });
  }

  /// Add constraint names
  void SetRowNames(std::vector<const char *> nm) {
    row_names_=std::move(nm);
    nlme_.SetRowNames(row_names_.data());
  }

  /// Add linear objective (only single objective supported.)
  /// Sense: NLW2_ObjSenseM....
  /// Coefficients: dense vector.
  void SetLinearObjective(int sense, double c0,
                          py::array_t<const double> c) {
    obj_sense_=sense; obj_c0_=c0; obj_c_=c;
    nlme_.SetLinearObjective(sense, c0, obj_c_.data());
  }

  /// Add Q for the objective quadratic part 0.5 @ x.T @ Q @ x.
  /// Format: NLW2_HessianFormat...
  void SetHessian(int nr,
                  int format,     // TODO enum
                  size_t nnz,
                  py::array_t<const size_t> st,
                  /// Entry index
                  py::array_t<const int> ind,
                  /// Entry value
                  py::array_t<const double> val
                  ) {
    Q_format_ = format;
    Q_={
      nr, 0,
      nnz, st, ind, val
    };
    nlme_.SetHessian(format, {
      nr, 0, nnz,
      Q_.start_.data(), Q_.index_.data(),
      Q_.value_.data()
    });
  }

  /// Set obj name
  void SetObjName(const char* nm) {
    obj_name_=(nm ? nm : "");
    nlme_.SetObjName(obj_name_);
  }

  /// Get the model
  const mp::NLModel& GetModel() const { return nlme_; }

private:
  /// Store the strings/arrays to keep the memory
  const char* prob_name_ {"NLWPY_Model"};
  mp::NLModel nlme_;
  NLWPY_ColData vars_ {};
  std::vector<const char *> var_names_ {};
  NLWPY_SparseMatrix A_ {};
  int num_row_ {};
  py::array_t<const double> row_lb_ {};
  py::array_t<const double> row_ub_ {};
  std::vector<const char *> row_names_ {};
  int obj_sense_ {};
  double obj_c0_ {};
  py::array_t<const double> obj_c_ {};
  int Q_format_ {};
  NLWPY_SparseMatrix Q_ {};
  const char* obj_name_ {"obj[1]"};
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
