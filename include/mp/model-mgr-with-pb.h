/*
 Model Manager with Std Problem Representation (mp::Problem)

 Copyright (C) 2020-2022 AMPL Optimization Inc

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

#ifndef MODEL_MANAGER_STD_H_
#define MODEL_MANAGER_STD_H_

#include <cmath>

#include "mp/env.h"
#include "mp/clock.h"
#include "mp/model-mgr-base.h"
#include "mp/solver-io.h"

namespace mp {

/// Model Manager with a given ProblemBuilder as intermediate instance storage.
///
/// Reads a problem from an .nl file, writes solution(s).
///
/// ModelManagerWithProblemBuilder uses Converter to convert the instance,
/// e.g., walk the NL tree stored in a given ProblemBuilder which can be,
/// e.g., mp::Problem.
/// The Converter must implement the BasicConverter<> interface.
///
/// TODO consider walking the tree while reading
template <class Converter>
class ModelManagerWithProblemBuilder :
    public BasicModelManager,
    public EnvKeeper
{
protected:
  using ConverterType = Converter;

public:
  /// Class name
  static const char* GetName() { return "ModelManagerWithFlatBackend TODO"; }

  ModelManagerWithProblemBuilder(std::unique_ptr<Converter> pc)
    : EnvKeeper(pc->GetEnv()), pcvt_(std::move(pc)) { }


protected:
  void SetBasename(const std::string& filename_no_ext) override {
    MakeUpTemporarySolHandler(filename_no_ext); ///< So that Abort() works
  }

  /// Before reading the NL file, a generic solution handler
  /// for error reporting
  void MakeUpTemporarySolHandler(const std::string& filename_no_ext) {
    ArrayRef<int> options({1, 1, 0});        ///< 'Default' NL options
    SetSolHandler(new internal::AppSolutionHandlerImpl<SolverType, ModelType>(
                             filename_no_ext, GetEnv(), GetModel(),
                             options,
                             GetEnv().get_output_handler().has_output ?
                                 0 :
                                 GetEnv().get_output_handler().banner_size));
  }

  void ReadNLFileAndUpdate(const std::string& nl_filename,
                           const std::string& filename_no_ext) override {
    steady_clock::time_point start = steady_clock::now();

    ReadNLFile(nl_filename);

    double read_time = GetTimeAndReset(start);
    if (GetEnv().timing())
      GetEnv().Print("NL model read time = {:.6f}s\n", read_time);

    MakeProperSolutionHandler(filename_no_ext);
    ConvertModelAndUpdateBackend();

    double cvt_time = GetTimeAndReset(start);
    if (GetEnv().timing())
      GetEnv().Print("NL model conversion time = {:.6f}s\n", cvt_time);
  }

  void ReadNLFile(const std::string& nl_filename) {
    set_nl_read_result_handler(
          new SolverNLHandlerType(GetModelBuilder(), GetEnv()));
    internal::NLFileReader<> reader;
    reader.Read(nl_filename, *nl_read_result_.handler_, 0);
  }

  /// Once NL is read
  void MakeProperSolutionHandler(const std::string& filename_no_ext) {
    ArrayRef<int> options(get_nl_read_result_handler().options(),
                          get_nl_read_result_handler().num_options());
    SetSolHandler(
          new internal::AppSolutionHandlerImpl<SolverType, ModelType>(
          filename_no_ext, GetEnv(), GetModelBuilder(), options,
          GetEnv().get_output_handler().has_output ?
            0 :
            GetEnv().get_output_handler().banner_size));
  }

  /// Says we finished problem modification,
  /// so we run model conversion and communicate the result into the backend
  void ConvertModelAndUpdateBackend() {
    GetCvt().ConvertModel();
  }


  ////////////////////////////// OPTIONS ////////////////////////////////
protected:
  // Add more text to be displayed before option descriptions.
  void add_to_long_name(fmt::StringRef name) { GetEnv().add_to_long_name(name); }
  void add_to_version(fmt::StringRef version) { GetEnv().add_to_version(version); }
  void add_to_option_header(const char* header_more) {
    GetEnv().add_to_option_header(header_more);
  }

  /// Simple stored option referencing a variable
  template <class Value>
  void AddOption(const char *name_list, const char *description,
                 Value& value, ValueArrayRef values = ValueArrayRef()) {
    GetEnv().AddStoredOption(name_list, description, value, values);
  }

  ////////////////////////////// UTILITIES //////////////////////////////
  template <typename... Args> \
  void Print(fmt::CStringRef format, const Args & ... args) {
    GetEnv().Print(format, args...);
  }

public:
  /// Initialize solver options
  void InitOptions() override {
    GetCvt().InitOptions();
    InitOwnOptions();
  }

private:
  void InitOwnOptions() {
  }

public:
  /// SolverType, needed by SolutionHandler
  using SolverType = BasicSolver;

  /// ModelType, needed by NL readers and solution handler
  using ModelType = typename ConverterType::ModelType;

  const ModelType& GetModelBuilder() const { return GetModel(); }
  ModelType& GetModelBuilder() { return GetModel(); }

  const ModelType& GetModel() const
  { return GetCvt().GetModel(); }
  ModelType& GetModel() { return GetCvt().GetModel(); }

  const Converter& GetCvt() const { return *pcvt_; }
  Converter& GetCvt() { return *pcvt_; }


protected:
  using SolverNLHandlerType =
    internal::SolverNLHandlerImpl<BasicSolver, ModelType,
                                  internal::NLProblemBuilder<ModelType>>;

  struct NLReadResult {
    using HandlerType = SolverNLHandlerType;
    std::unique_ptr<HandlerType> handler_;
  };

  void set_nl_read_result_handler(typename NLReadResult::HandlerType* ph)
  { nl_read_result_.handler_.reset(ph); }
  const typename NLReadResult::HandlerType& get_nl_read_result_handler() const
  { return *nl_read_result_.handler_; }
  typename NLReadResult::HandlerType& get_nl_read_result_handler()
  { return *nl_read_result_.handler_; }

  bool HaveSolH() const { return (bool)p_sol_handler_; }
  const SolutionHandler& GetSolH() const
  { assert(p_sol_handler_); return *p_sol_handler_; }
  SolutionHandler& GetSolH()
  { assert(p_sol_handler_); return *p_sol_handler_; }
  void SetSolHandler(SolutionHandler* psh)
  { assert(psh); p_sol_handler_.reset(psh); }
  void RemoveSolHandler() { p_sol_handler_.reset(nullptr); }

  const BasicModelManager& GetCQ() const { return *this; }
  BasicModelManager& GetCQ() { return *this; }

  /// Define original model access from BasicOriginalModelProxy
  ArrayRef<double> InitialValues() override {
    return GetModel().InitialValues();
  }
  ArrayRef<double> InitialDualValues() override {
    return GetModel().InitialDualValues();
  }

  ArrayRef<int> ReadSuffix(const SuffixDef<int>& suf) override {
    return GetModel().ReadSuffix(suf);
  }
  ArrayRef<double> ReadSuffix(const SuffixDef<double>& suf) override {
    return GetModel().ReadSuffix(suf);
  }

  void ReportSuffix(const SuffixDef<int>& suf,
                    ArrayRef<int> values) override {
    GetModel().ReportSuffix(suf, values);
  }
  void ReportSuffix(const SuffixDef<double>& suf,
                    ArrayRef<double> values) override {
    GetModel().ReportSuffix(suf, values);
  }
  size_t GetSuffixSize(int kind) override {
    return GetModel().GetSuffixSize((suf::Kind)kind);
  }

  void HandleSolution(int status, fmt::CStringRef msg,
      const double *x, const double *y, double obj) override {
    if ( HaveSolH() )
      GetSolH().HandleSolution(status, msg, x, y, obj);
    else
      MP_RAISE_WITH_CODE(0, msg);
  }

  void HandleFeasibleSolution(fmt::CStringRef msg,
      const double *x, const double *y, double obj) override {
    GetSolH().HandleFeasibleSolution(msg, x, y, obj);
  }


  /// TODO Why here?
  const std::vector<bool>& IsVarInt() const override {
    return GetModel().IsVarInt();
  }


private:
  std::unique_ptr<Converter> pcvt_;
  NLReadResult nl_read_result_;
  std::unique_ptr<SolutionHandler> p_sol_handler_;
};

}  // namespace mp

#endif  // MODEL_MANAGER_STD_H_
