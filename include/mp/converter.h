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

#ifndef CONVERTER_H_
#define CONVERTER_H_

#include <mp/problem.h>
#include <mp/backend.h>
#include <mp/solver.h>

namespace mp {

/// An abstract MP converter - does not change anything
/// Responsible for model modification and solving, typical 'exported' solver API
/// Backend access is hidden (the backend itself is a parameter)
template <class Impl, class Backend,
          class Model = BasicProblem<std::allocator<char> > >
class MPConverter {
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
  Model& GetModel() { return model_; }    // Can be used for NL file input
  Backend& GetBackend() { return backend_; }
public:

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
    ReadNLFile(nl_filename, nl_reader_flags);
    UpdateBackend();
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
  void UpdateBackend() {
    MP_DISPATCH( ConvertModel() );
    MP_DISPATCH( PushChangesToBackend() );
  }

  void Solve(SolutionHandler &sh) {
    GetBackend().Solve(GetModel(), sh);   // TODO no model any more
  }


protected:
  void PushChangesToBackend() {
    ModelToBackendFeeder<Model, Backend>
        feeder(GetModel(), GetBackend());
    feeder.PushWholeProblem();
  }
};

/// One of the converters requiring a "minimal" output interface
template <class Impl, class Backend,
          class Model = BasicProblem<std::allocator<char> > >
class BasicMPToMIPConverter : public MPConverter<Impl, Backend, Model> {
public:
  void ConvertModel() { /* DO NOTHING YET */ }

};

template <class Impl, template <typename, typename, typename> class Converter,
          class Backend, class Model = BasicProblem<std::allocator<char> > >
class BasicInterface : public Converter<Impl, Backend, Model> {
};

}  // namespace mp

#endif  // CONVERTER_H_
