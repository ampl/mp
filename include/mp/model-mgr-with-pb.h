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

/// Model Manager with a ProblemBuilder.
///
/// Reads a problem from an .nl file, writes solution(s).
///
/// ProblemBuilder is defined by the Converter
/// and used as an intermediate instance storage.
///
/// ModelManagerWithProblemBuilder uses Converter to convert the instance,
/// e.g., walk the NL tree and convert other constraints
/// stored in a given ProblemBuilder which can be, e.g., mp::Problem.
///
/// The Converter must implement the BasicConverter<ProblemBuilder> interface.
///
/// TODO consider walking the tree while reading
template <class Converter>
class ModelManagerWithProblemBuilder :
    public BasicModelManager,
    public EnvKeeper
{
protected:
  /// Convenience typedef
  using ConverterType = Converter;

  /// SolverType, needed by SolutionHandler
  using SolverType = BasicSolver;

  /// ProblemBuilder, needed by NL readers and solution handler.
  /// Converter does not have to know this is a ProblemBuilder.
  using ProblemBuilder = typename ConverterType::ModelType;


public:
  /// Class name
  static const char* GetTypeName() { return "ModelManagerWithProblemBuilder"; }

  ModelManagerWithProblemBuilder(std::unique_ptr<Converter> pc)
    : EnvKeeper(pc->GetEnv()), pcvt_(std::move(pc)) { }


protected:
  /// Initialize solver options
  void InitOptions() override {
    GetCvt().InitOptions();
    InitOwnOptions();
  }

  void SetBasename(const std::string& filename_no_ext) override {
     /// So that Abort() + sol writing work
     /// before we have parsed the NL file
     MakeUpTemporarySolHandler(filename_no_ext);
  }

  void ReadNLModel(const std::string& nl_filename,
                   const std::string& filename_no_ext,
                   void (*cb_checkmodel)(size_t, size_t, size_t))
  override {
    steady_clock::time_point start = steady_clock::now();

    ReadNLFile(nl_filename);

    if (cb_checkmodel)
      cb_checkmodel(GetPB().num_vars(),
                    GetPB().num_algebraic_cons(),
                    GetPB().num_logical_cons());

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
          new SolverNLHandlerType(GetPB(), GetEnv()));
    internal::NLFileReader<> reader;
    reader.Read(nl_filename, *nl_read_result_.handler_, 0);
  }

  /// Before reading the NL file, a generic solution handler
  /// for error reporting
  void MakeUpTemporarySolHandler(const std::string& filename_no_ext) {
    static int opt_static[] = {1, 1, 0};
    ArrayRef<int> options(opt_static);             ///< 'Default' NL options
    SetSolHandler(new internal::AppSolutionHandlerImpl<SolverType, ProblemBuilder>(
                             filename_no_ext, GetEnv(), GetModel(),
                             options,
                             GetEnv().get_output_handler().has_output ?
                                 0 :
                                 GetEnv().get_output_handler().banner_size));
  }

  /// Once NL is read
  void MakeProperSolutionHandler(const std::string& filename_no_ext) {
    ArrayRef<int> options(get_nl_read_result_handler().options(),
                          get_nl_read_result_handler().num_options());
    SetSolHandler(
          new internal::AppSolutionHandlerImpl<SolverType, ProblemBuilder>(
          filename_no_ext, GetEnv(), GetPB(), options,
          GetEnv().get_output_handler().has_output ?
            0 :
            GetEnv().get_output_handler().banner_size));
  }

  /// Says we finished problem modification,
  /// so we run model conversion and communicate the result into the backend
  void ConvertModelAndUpdateBackend() {
    GetCvt().ConvertModel();
  }


protected:
  const ProblemBuilder& GetPB() const { return GetModel(); }
  ProblemBuilder& GetPB() { return GetModel(); }

  const ProblemBuilder& GetModel() const
  { return GetCvt().GetModel(); }
  ProblemBuilder& GetModel() { return GetCvt().GetModel(); }

  const Converter& GetCvt() const { return *pcvt_; }
  Converter& GetCvt() { return *pcvt_; }

  using SolverNLHandlerType =
    internal::SolverNLHandlerImpl< BasicSolver, ProblemBuilder,
                                   internal::NLProblemBuilder<ProblemBuilder> >;

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
  void InitOwnOptions() {
  }


private:
  std::unique_ptr<Converter> pcvt_;
  NLReadResult nl_read_result_;
  std::unique_ptr<SolutionHandler> p_sol_handler_;
};

}  // namespace mp

#endif  // MODEL_MANAGER_STD_H_
