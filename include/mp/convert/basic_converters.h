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

#include "mp/problem.h"
#include "mp/convert/model_adapter.h"
#include "mp/convert/backend.h"
#include "mp/solver.h"

namespace mp {

/// An abstract MP converter - only complains, all conversions need to be redefined
/// in derived classes.
/// Responsible for model modification and solving, typical 'exported' solver API
/// Backend access is hidden (the backend itself is a parameter)
template <class Impl, class Backend,
          class Model = BasicProblem< > >
class BasicMPConverter :
    public BasicConstraintConverter {

  ModelAdapter<Model> model_adapter_;
  Backend backend_;
  /// This is to wrap some dependencies from MP
  using SolverAdapter = SolverImpl< ModelAdapter<Model> >;

  std::unique_ptr<ConverterQuery> p_converter_query_;
  SolutionHandler* p_sol_h_;

public:
  static const char* GetConverterName() { return "BasicMPConverter"; }
  using Converter = Impl;
  using ModelType = Model;
  using OutputModelType = ModelAdapter<Model>;
  using BackendType = Backend;

  /// MP API requires a 'ProblemBuilder' type
  using ProblemBuilder = OutputModelType;
  using MPUtils = typename Backend::MPUtils;

  /// The working model
  const Model& GetModel() const { return GetOutputModel().GetModel(); }
  /// The working model
  Model& GetModel() { return GetOutputModel().GetModel(); }

  /// Dirty: returning the output model
  ///            // Can be used for NL file input
  const OutputModelType& GetInputModel() const { return GetOutputModel(); }
  OutputModelType& GetInputModel() { return GetOutputModel(); }

  const OutputModelType& GetOutputModel() const { return model_adapter_; }   // TODO
  OutputModelType& GetOutputModel() { return model_adapter_; }

  const Backend& GetBackend() const { return backend_; }
  Backend& GetBackend() { return backend_; }

  const MPUtils& GetMPUtils() const { return GetBackend().GetMPUtils(); }
  MPUtils& GetMPUtils() { return GetBackend().GetMPUtils(); }

  const ConverterQuery& GetCQ() const {
    assert(p_converter_query_);
    return *p_converter_query_;
  }
  ConverterQuery& GetCQ() {
    assert(p_converter_query_);
    return *p_converter_query_;
  }

  const SolutionHandler& GetSolH() const { assert(p_sol_h_); return *p_sol_h_; }
  SolutionHandler& GetSolH() { assert(p_sol_h_); return *p_sol_h_; }
  void SetSolHandler(SolutionHandler& sh) { assert(&sh); p_sol_h_ = &sh; }
  void RemoveSolHandler() { p_sol_h_=nullptr; }

public:

  BasicMPConverter() {
    InitConverterQueryObject();
    InitOptions();
    GetBackend().InitMetaInfoAndOptions();
  }

  void InitConverterQueryObject() {
    p_converter_query_ = MP_DISPATCH( MakeConverterQuery() );
    GetBackend().ProvideConverterQueryObject( &MP_DISPATCH( GetCQ() ) );
  }

  bool ParseOptions(char **argv, unsigned flags = 0) {
    return GetMPUtils().ParseOptions(argv, flags);
  }

  struct NLReadResult {
    std::unique_ptr< internal::SolverNLHandler<SolverAdapter> > handler_;
  };
  NLReadResult ReadNLFile(const std::string& nl_filename, int nl_reader_flags) {
    NLReadResult result;
    result.handler_.reset(
          new internal::SolverNLHandler<SolverAdapter>(GetInputModel(), GetMPUtils()));
    internal::NLFileReader<> reader;
    reader.Read(nl_filename, *result.handler_, nl_reader_flags);
    return result;
  }
  NLReadResult ReadNLFileAndUpdate(const std::string& nl_filename, int nl_reader_flags) {
    NLReadResult result = ReadNLFile(nl_filename, nl_reader_flags);
    ConvertModelAndUpdateBackend();
    return result;
  }

  /// INCREMENTAL INTERFACE
  /// These guys used from outside to feed a model to be converted
  /// and forwarded to a backend
  /// Currently only used for testing
  void InputVariables(int n, const double* lb, const double* ub, const var::Type* ty) {
    GetModel().AddVars(n, lb, ub, ty);
  }
  void InputObjective(obj::Type t,
                      int nnz, const double* c, const int* v, NumericExpr e=NumericExpr()) {
    typename Model::LinearObjBuilder lob = GetModel().AddObj(t, e);
    for (int i=0; i!=nnz; ++i) {
      lob.AddTerm(v[i], c[i]);
    }
  }
  void InputAlgebraicCon(int nnz, const double* c, const int* v,
                         double lb, double ub, NumericExpr e=NumericExpr()) {
    typename Model::MutAlgebraicCon mac = GetModel().AddCon(lb, ub);
    typename Model::LinearConBuilder lcb = mac.set_linear_expr(nnz);
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
    SetSolHandler(sh);
    GetBackend().SolveAndReport();
    RemoveSolHandler();
  }


protected:
  /// Convert the whole model, e.g., after reading from NL
  void ConvertModel() {
    MP_DISPATCH( PrepareConversion() );
    MP_DISPATCH( ConvertStandardItems() );
    MP_DISPATCH( ConvertExtraItems() );
  }

  void PrepareConversion() {
    MP_DISPATCH( MemorizeModelSize() );
  }

  void MemorizeModelSize() {
    GetOutputModel().set_num_vars(GetModel().num_vars());
    GetOutputModel().set_num_alg_cons(GetModel().num_algebraic_cons());
  }

  void ConvertStandardItems() {
    int num_common_exprs = GetModel().num_common_exprs();
    for (int i = 0; i < num_common_exprs; ++i)
      MP_DISPATCH( Convert( GetModel().common_expr(i) ) );
    if (int num_objs = GetModel().num_objs())
      for (int i = 0; i < num_objs; ++i)
        MP_DISPATCH( Convert( GetModel().obj(i) ) );
    if (int n_cons = GetModel().num_algebraic_cons())
      for (int i = 0; i < n_cons; ++i)
        MP_DISPATCH( Convert( GetModel().algebraic_con(i) ) );
    if (int n_lcons = GetModel().num_logical_cons())
      for (int i = 0; i < n_lcons; ++i)
        MP_DISPATCH( Convert( GetModel().logical_con(i) ) );

    ////////////////////////////////////////////////
    MP_DISPATCH( ConvertSOSConstraints() );
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

  void ConvertSOSConstraints() {
    if (sos()) {
      auto sosno = GetInputModel().
          ReadIntSuffix( { "sosno", suf::Kind::VAR } );
      auto ref = GetInputModel().
          ReadDblSuffix( { "ref", suf::Kind::VAR } );
      if (sosno && ref)
        MP_DISPATCH( ConvertSOSCollection(sosno, ref, false) );
    }
    if (sos2_ampl_pl()) {
      auto sos = GetInputModel().
          ReadIntSuffix( { "sos", suf::Kind::VAR } );
      auto sosref = GetInputModel().
          ReadDblSuffix( { "sosref", suf::Kind::VAR } );
      if (sos && sosref)
        MP_DISPATCH( ConvertSOSCollection(sos, sosref, true) );
    }
  }

  /// fAllSOS2: if false, only groups with sosno<0
  void ConvertSOSCollection(ArrayRef<int> sosno, ArrayRef<double> ref,
                            bool fAllSOS2) {
    assert(sosno.size() == ref.size());
    internal::Unused(fAllSOS2);
    UNSUPPORTED("ConvertSOSCollection");
  }

  void PushWholeModelToBackend() {
    GetModel().PushModelTo(GetBackend());
  }

public:
  /// These methods to be used by converter helper objects
  static constexpr double Infty() { return std::numeric_limits<double>::infinity(); }
  static constexpr double MinusInfty() { return -std::numeric_limits<double>::infinity(); }
  int AddVar(double lb=MinusInfty(), double ub=Infty(), var::Type type = var::CONTINUOUS) {
    auto var = GetModel().AddVar(lb, ub, type);
    return var.index();
  }
  std::vector<int> AddVars(std::size_t nvars, double lb, double ub, var::Type type = var::CONTINUOUS) {
    std::vector<int> newVars(nvars);
    for (std::size_t  i=0; i<nvars; ++i)
      newVars[i] = AddVar(lb, ub, type);
    return newVars;
  }

  ////////////////////////////// OPTIONS ////////////////////////////////
protected:
  // Add more text to be displayed before option descriptions.
  void add_to_long_name(fmt::StringRef name) { GetBackend().add_to_long_name(name); }
  void add_to_version(fmt::StringRef version) { GetBackend().add_to_version(version); }
  void add_to_option_header(const char* header_more) {
    GetBackend().add_to_option_header(header_more);
  }

  /// Simple stored option referencing a variable
  template <class Value>
  void AddOption(const char *name_list, const char *description,
                 Value& value, ValueArrayRef values = ValueArrayRef()) {
    GetBackend().AddStoredOption(name_list, description, value, values);
  }

  ////////////////////////////// UTILITIES //////////////////////////////
  template <typename... Args> \
  void Print(fmt::CStringRef format, const Args & ... args) {
    GetMPUtils().Print(format, args...);
  }

private:
  struct Options {
    int sos_ = 1;
    int sos2_ = 1;
  };
  Options options_;

protected:
  int sos() const { return options_.sos_; }
  int sos2_ampl_pl() const { return options_.sos2_; }

private:
  void InitOptions() {
    this->AddOption("cvt:sos sos",
        "0/1*: Whether to honor declared suffixes .sosno and .ref describing "
        "SOS sets. Each distinct nonzero .sosno "
        "value designates an SOS set, of type 1 for "
        "positive .sosno values and of type 2 for "
        "negative values.  The .ref suffix contains "
        "corresponding reference values used to order the variables.",
        options_.sos_);
    this->AddOption("cvt:sos2 sos2",
        "0/1*: Whether to honor SOS2 constraints for nonconvex "
        "piecewise-linear terms, using suffixes .sos and .sosref "
        "provided by AMPL.",
        options_.sos_);
  }

};

/// A null converter - does not change anything
template <class Impl, class Backend,
          class Model = BasicProblem< > >
class NullMPConverter : public BasicMPConverter<Impl, Backend, Model> {
public:
  void ConvertModel() { }
};

/// A 'final' converter in a hierarchy
template <template <typename, typename, typename> class Converter,
          class Backend, class Model = BasicModel< > >
class ConverterImpl :
    public Converter<ConverterImpl<Converter, Backend, Model>, Backend, Model> { };

template <template <typename, typename, typename> class Converter,
          class Backend, class Model = BasicModel< > >
using Interface = ConverterImpl<Converter, Backend, Model>;

/// Conversion failure helper
class ConstraintConversionFailure {
  const std::string msg_;
public:
  ConstraintConversionFailure(std::string&& msg) noexcept :
    msg_(std::move(msg)) { }
  const std::string& message() const { return msg_; }
};

}  // namespace mp

#endif  // BASIC_CONVERTERS_H_
