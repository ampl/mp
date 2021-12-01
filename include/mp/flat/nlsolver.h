/*
 NL Solver control class

 Copyright (C) 2020-2021 AMPL Optimization Inc

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

#include <cmath>

#include "mp/flat/expr_flattener.h"
#include "mp/flat/backend_base.h"
#include "mp/solver.h"

namespace mp {

/// NL solver control class.
///
/// Reads a problem from an .nl file, parses solver options
/// from argv and environment variables, solves the problem and writes
/// solution(s)
///
/// BasicNLSolverWithFlatBackend uses ExprFlattener to flatten the NL tree
/// TODO consider walking the tree while reading
template <class ExprFlattener>
class BasicNLSolverWithFlatBackend :
    public NLSolverProxy
{
protected:
  using ExprFlattenerType = ExprFlattener;

public:
  static const char* GetName() { return "BasicNLSolverWithFlatBackend"; }

  BasicNLSolverWithFlatBackend() { }

  /// Parse solver options such as "outlev=1".
  /// @param flags: 0 or \a Solver::NO_OPTION_ECHO
  bool ParseOptions(const char* filename_no_ext,
                    char **argv, unsigned flags = 0) {
    /// Chance e.g. for the Backend to init solver environment, etc
    InitOptionParsing(filename_no_ext);
    if (GetMPUtils().ParseOptions(argv, flags)) {
      /// Chance to consider options immediately (open cloud, etc)
      GetExprFlattener().FinishOptionParsing();
      return true;
    }
    return false;
  }

  /// Runs Solver given the NL file name.
  void RunFromNLFile(const std::string& nl_filename,
                     const std::string& filename_no_ext,
                     int nl_flags) {
    ReadNLFileAndUpdate(nl_filename, filename_no_ext, nl_flags);
    Solve();
  }

protected:
  void InitOptionParsing(const std::string& filename_no_ext) {
    InitConverterQueryObject();                 ///< For the Backend
    MakeUpTemporarySolHandler(filename_no_ext); ///< So that Abort() works
    GetExprFlattener().InitOptionParsing();
  }

  void InitConverterQueryObject() {
    GetBasicBackend().ProvideNLSolverProxyObject( &GetCQ() );
  }

  /// Before reading the NL file, a generic solution handler
  /// for error reporting
  void MakeUpTemporarySolHandler(const std::string& filename_no_ext) {
    ArrayRef<int> options({1, 1, 0});        ///< 'Default' NL options
    SetSolHandler(new internal::AppSolutionHandler<MPUtils>(
                             filename_no_ext, GetMPUtils(), GetModelBuilder(),
                             options,
                             GetMPUtils().get_output_handler().has_output ?
                                 0 :
                                 GetMPUtils().get_output_handler().banner_size));
  }

  void ReadNLFileAndUpdate(const std::string& nl_filename,
                           const std::string& filename_no_ext,
                           int nl_reader_flags) {
    steady_clock::time_point start = steady_clock::now();

    ReadNLFile(nl_filename, nl_reader_flags);

    double read_time = GetTimeAndReset(start);
    if (GetMPUtils().timing())
      GetMPUtils().Print("NL model read time = {:.6f}s\n", read_time);

    MakeProperSolutionHandler(filename_no_ext);
    ConvertModelAndUpdateBackend();

    double cvt_time = GetTimeAndReset(start);
    if (GetMPUtils().timing())
      GetMPUtils().Print("NL model conversion time = {:.6f}s\n", cvt_time);
  }

  void ReadNLFile(const std::string& nl_filename, int nl_reader_flags) {
    set_nl_read_result_handler(
          new internal::SolverNLHandler<
            SolverAdapter>(GetModelBuilder(), GetMPUtils()));
    internal::NLFileReader<> reader;
    reader.Read(nl_filename, *nl_read_result_.handler_, nl_reader_flags);
  }

  /// Once NL is read
  void MakeProperSolutionHandler(const std::string& filename_no_ext) {
    ArrayRef<int> options(get_nl_read_result_handler().options(),
                          get_nl_read_result_handler().num_options());
    SetSolHandler(new internal::AppSolutionHandler<MPUtils>(
          filename_no_ext, GetMPUtils(), GetModelBuilder(), options,
          GetMPUtils().get_output_handler().has_output ?
            0 :
            GetMPUtils().get_output_handler().banner_size));
  }

  /// Says we finished problem modification,
  /// so we run model conversion and communicate the result into the backend
  void ConvertModelAndUpdateBackend() {
    GetExprFlattener().ConvertModel();
  }

  void Solve() {
    GetBasicBackend().SolveAndReport();
  }


  ////////////////////////////// OPTIONS ////////////////////////////////
protected:
  // Add more text to be displayed before option descriptions.
  void add_to_long_name(fmt::StringRef name) { GetMPUtils().add_to_long_name(name); }
  void add_to_version(fmt::StringRef version) { GetMPUtils().add_to_version(version); }
  void add_to_option_header(const char* header_more) {
    GetMPUtils().add_to_option_header(header_more);
  }

  /// Simple stored option referencing a variable
  template <class Value>
  void AddOption(const char *name_list, const char *description,
                 Value& value, ValueArrayRef values = ValueArrayRef()) {
    GetMPUtils().AddStoredOption(name_list, description, value, values);
  }

  ////////////////////////////// UTILITIES //////////////////////////////
  template <typename... Args> \
  void Print(fmt::CStringRef format, const Args & ... args) {
    GetMPUtils().Print(format, args...);
  }

public:
  void InitOptions() {
    GetExprFlattener().InitOptions();
    InitOwnOptions();
  }

private:
  void InitOwnOptions() {
  }

public:
  /// TODO use universal Env instead
  /// (which can well use these "MPUtils",
  /// but ideally an appr new base class of Solver)
  using MPUtils = typename ExprFlattenerType::MPUtils;
  const MPUtils& GetMPUtils() const { return GetExprFlattener().GetMPUtils(); }
  MPUtils& GetMPUtils() { return GetExprFlattener().GetMPUtils(); }

protected:
  using ModelType = typename ExprFlattenerType::ModelType;

  const ModelType& GetModelBuilder() const { return GetModel(); }
  ModelType& GetModelBuilder() { return GetModel(); }

  const ModelType& GetModel() const
  { return GetExprFlattener().GetModel(); }
  ModelType& GetModel() { return GetExprFlattener().GetModel(); }

  const ExprFlattener& GetExprFlattener() const { return expr_flt_; }
  ExprFlattener& GetExprFlattener() { return expr_flt_; }

  /// Expose abstract Backend from FlatCvt
  const BasicBackend& GetBasicBackend() const
  { return GetExprFlattener().GetBasicBackend(); }
  BasicBackend& GetBasicBackend()
  { return GetExprFlattener().GetBasicBackend(); }

  /// This is to wrap some dependencies from MP
  /// TODO hide
  using SolverAdapter = SolverImpl< ModelType >;

  struct NLReadResult {     // TODO decouple from SolverAdapter
    using HandlerType = internal::SolverNLHandler<SolverAdapter>;
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

  const NLSolverProxy& GetCQ() const { return *this; }
  NLSolverProxy& GetCQ() { return *this; }

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
  ExprFlattener expr_flt_;
  NLReadResult nl_read_result_;
  std::unique_ptr<SolutionHandler> p_sol_handler_;
};

template <class Backend,
          template <typename, typename, typename> class Converter,
          class Model = Problem >
using NLSolverWithFlatBackend = BasicNLSolverWithFlatBackend<
                                  ExprFlattenerImpl<ExprFlattener, Model,
                                    Interface<Converter, Backend> > >;

}  // namespace mp

#endif  // BASIC_CONVERTERS_H_
