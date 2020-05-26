/*
 Abstract and basic model converters.

 Copyright (C) 2020 AMPL Optimization Inc

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

 Author: Gleb Belov <gleb.belov@monash.edu>
 */

#ifndef BASIC_CONVERTERS_H_
#define BASIC_CONVERTERS_H_

#include "mp/convert/model.h"
#include "mp/backend.h"
#include "mp/solver.h"

namespace mp {

/// An abstract MP converter - only complains, all conversions need to be redefined
/// in derived classes.
/// Responsible for model modification and solving, typical 'exported' solver API
/// Backend access is hidden (the backend itself is a parameter)
template <class Impl, class Backend,
          class Model = BasicProblem<std::allocator<char> > >
class BasicMPConverter {
protected:
  /// This is to wrap some old dependencies from MP
  using SolverAdapter = SolverImpl<Model>;
public:
  Model model_;
  Backend backend_;
public:
  using Converter = Impl;
  using ModelType = Model;
  using ProblemBuilder = Model;           // for old MP stuff
  using BackendType = Backend;
  const Model& GetModel() const { return model_; }    // Can be used for NL file input
  Model& GetModel() { return model_; }    // Can be used for NL file input
  const Backend& GetBackend() const { return backend_; }
  Backend& GetBackend() { return backend_; }
public:

  bool ParseOptions(char **argv, unsigned flags = 0) {
    return GetBackend().ParseOptions(argv, flags);
  }

  struct NLReadResult {
    std::unique_ptr< internal::SolverNLHandler<SolverAdapter> > handler_;
  };
  NLReadResult ReadNLFile(const std::string& nl_filename, int nl_reader_flags) {
    NLReadResult result;
    result.handler_.reset(
          new internal::SolverNLHandler<SolverAdapter>(GetModel(), GetBackend()));
    internal::NLFileReader<> reader;
    reader.Read(nl_filename, *result.handler_, nl_reader_flags);
    return result;
  }
  NLReadResult ReadNLFileAndUpdate(const std::string& nl_filename, int nl_reader_flags) {
    NLReadResult result = ReadNLFile(nl_filename, nl_reader_flags);
    ConvertModelAndUpdateBackend();
    return result;
  }

  /// These guys used from outside to feed a model to be converted
  /// and forwarded to a backend
  void InputVariables(int n, const double* lb, const double* ub, const var::Type* ty) {
    model_.AddVars(n, lb, ub, ty);
  }
  void InputObjective(obj::Type t,
                      int nnz, const double* c, const int* v, NumericExpr e=NumericExpr()) {
    mp::Problem::LinearObjBuilder lob = model_.AddObj(t, e);
    for (int i=0; i!=nnz; ++i) {
      lob.AddTerm(c[i], v[i]);
    }
  }
  void InputAlgebraicCon(int nnz, const double* c, const int* v,
                         double lb, double ub, NumericExpr e=NumericExpr()) {
    mp::Problem::MutAlgebraicCon mac = model_.AddCon(lb, ub);
    mp::Problem::LinearConBuilder lcb = mac.set_linear_expr(nnz);
    for (int i=0; i!=nnz; ++i)
      lcb.AddTerm(v[i], c[i]);
    mac.set_nonlinear_expr(e);
  }

  /// Says we finished problem modification,
  /// so we run model conversion and communicate the result into the backend
  void ConvertModelAndUpdateBackend() {
    MP_DISPATCH( ConvertModel() );
    MP_DISPATCH( PushWholeModelToBackend() );
  }

  void Solve(SolutionHandler &sh) {
    GetBackend().Solve(GetModel(), sh);   // TODO no model any more
  }


protected:
  /// Convert the whole model, e.g., after reading from NL
  void ConvertModel() {
    MP_DISPATCH( ConvertStandardItems() );
    MP_DISPATCH( ConvertExtraItems() );
  }

  void ConvertStandardItems() {
    int num_common_exprs = model_.num_common_exprs();
    for (int i = 0; i < num_common_exprs; ++i)
      MP_DISPATCH( Convert( model_.common_expr(i) ) );
    if (int num_objs = model_.num_objs())
      for (int i = 0; i < num_objs; ++i)
        MP_DISPATCH( Convert( model_.obj(i) ) );
    if (int n_cons = model_.num_algebraic_cons())
      for (int i = 0; i < n_cons; ++i)
        MP_DISPATCH( Convert( model_.algebraic_con(i) ) );
    if (int n_lcons = model_.num_logical_cons())
      for (int i = 0; i < n_lcons; ++i)
        MP_DISPATCH( Convert( model_.logical_con(i) ) );
  }

  void ConvertExtraItems() { }

  void Convert(typename Model::MutCommonExpr e) {
    throw std::runtime_error("No common exprs conversion implemented");
  }

  void Convert(typename Model::MutObjective obj) {
    throw std::runtime_error("No objectives conversion implemented");
  }

  void Convert(typename Model::MutAlgebraicCon con) {
    throw std::runtime_error("No algebraic constraints conversion implemented");
  }

  void Convert(typename Model::MutLogicalCon e) {
    throw std::runtime_error("No logical constraints conversion implemented");
  }

  void PushWholeModelToBackend() {
    GetModel().PushModelTo(GetBackend());
  }

public:
  /// These methods to be used by converter objects
  static constexpr double Infty() { return std::numeric_limits<double>::infinity(); }
  static constexpr double MinusInfty() { return -std::numeric_limits<double>::infinity(); }
  int AddVar(double lb=MinusInfty(), double ub=Infty(), var::Type type = var::CONTINUOUS) {
    auto var = GetModel().AddVar(lb, ub, type);
    return var.index();
  }
  std::vector<int> AddVars(int nvars, double lb, double ub, var::Type type = var::CONTINUOUS) {
    std::vector<int> newVars(nvars);                          // The binary flags
    for (int i=0; i<nvars; ++i)
      newVars[i] = AddVar(lb, ub, type);
    return newVars;
  }

};

/// A null converter - does not change anything
template <class Impl, class Backend,
          class Model = BasicProblem<std::allocator<char> > >
class NullMPConverter : public BasicMPConverter<Impl, Backend, Model> {
public:
  void ConvertModel() { }
};

/// A 'final' converter in a hierarchy, no static polymorphism
template <template <typename, typename, typename> class Converter,
          class Backend, class Model = BasicModel<std::allocator<char> > >
class ConverterImpl :
    public Converter<ConverterImpl<Converter, Backend, Model>, Backend, Model> { };

template <template <typename, typename, typename> class Converter,
          class Backend, class Model = BasicModel<std::allocator<char> > >
using Interface = ConverterImpl<Converter, Backend, Model>;

}  // namespace mp

#endif  // BASIC_CONVERTERS_H_
