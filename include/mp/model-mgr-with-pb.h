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

  void ReadNLModel(const std::string& nl_filename,
                   const std::string& filename_no_ext,
                   Checker_AMPLS_ModeltTraits cb_checkmodel,
                   std::function<void()> after_header)
  override {
    steady_clock::time_point start = steady_clock::now();

    ReadNLFile(nl_filename,
               [this, &filename_no_ext, after_header](){
      MakeProperSolutionHandler(filename_no_ext);
      if (after_header)
        after_header();                   // parse options
    });
    ReadNames(filename_no_ext);

    double read_time = GetTimeAndReset(start);
    if (GetEnv().timing())
      GetEnv().Print("NL model read time = {:.6f}s\n", read_time);

    ConvertModelAndUpdateBackend();

    if (cb_checkmodel) {
      AMPLS_ModelTraits mtraits;
      GetCvt().FillModelTraits(mtraits);
      cb_checkmodel(&mtraits);
    }

    double cvt_time = GetTimeAndReset(start);
    if (GetEnv().timing())
      GetEnv().Print("NL model conversion time = {:.6f}s\n", cvt_time);
  }

  /// Read the NL file
  void ReadNLFile(
      const std::string& nl_filename,
      std::function<void()> after_header) {
    set_nl_read_result_handler(
          new SolverNLHandlerType(GetPB(), GetEnv(), after_header));
    internal::NLFileReader<> reader;
    reader.Read(nl_filename, *nl_read_result_.handler_, 0);
  }

  /// Read var / con / obj names.
  /// The .row file has cons + objs.
  void ReadNames(const std::string& namebase) {
    if (WantNames()) {
      NameProvider npv("_svar");
      NameProvider npc("_scon");
      if (WantNames()<=2) {
        npv.ReadNames(namebase + ".col",
                      GetModel().num_vars());
        npc.ReadNames(namebase + ".row",
                      GetModel().num_cons()
                      + GetModel().num_objs());
      }
      if (WantNames()>=2
          || npv.number_read()+npc.number_read()) {
        GetModel().SetVarNames(
              npv.get_names(GetModel().num_vars()));
        GetModel().SetConNames(
              npc.get_names(GetModel().num_cons()) );
        SetObjNames(npc);
      }
    }
  }

  /// Obj names.
  /// @param npco: NameProvider of row+obj names,
  ///   read from .row or generated
  ///
  /// We have to consider that obj:no=n
  /// selects objective n.
  void SetObjNames(NameProvider& npco) {
    if (GetModel().num_objs()) {
      auto num_c = GetModel().num_cons();
      auto o1 = GetEnv().objno_used()-1;
      assert(o1>=0);
      auto o2 = o1+1;
      if (GetEnv().multiobj()) {
        o1 = 0;
        o2 = GetModel().num_objs();
      }
      std::vector<std::string> names_o;
      for (auto io=num_c+o1; io<num_c+o2; ++io) {
        if (npco.number_read()>(size_t)io)
          names_o.push_back(npco.name(io));
        else
          names_o.push_back("_sobj["
                            + std::to_string(io-num_c+1)
                            + ']');
      }
      GetModel().SetObjNames( std::move( names_o ) );
    }
  }

  /// Once NL header is read
  void MakeProperSolutionHandler(const std::string& filename_no_ext) {
    ArrayRef<int> options(get_nl_read_result_handler().options(),
                          get_nl_read_result_handler().num_options());
    SetSolHandler(
          new internal::AppSolutionHandlerImpl<SolverType, ProblemBuilder>(
          filename_no_ext, GetEnv(), GetPB(), options,
            GetEnv().get_output_handler().has_output
            ? 0 : GetEnv().get_output_handler().banner_size));
  }

  /// Says we finished problem modification,
  /// so we run model conversion and communicate the result into the ModelAPI.
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
  ArrayRef<int> InitialValuesSparsity() override {
    return GetModel().InitialValuesSparsity();
  }
  ArrayRef<double> InitialDualValues() override {
    return GetModel().InitialDualValues();
  }
  ArrayRef<int> InitialDualValuesSparsity() override {
    return GetModel().InitialDualValuesSparsity();
  }

  ArrayRef<int> ReadSuffix(const SuffixDef<int>& suf) override {
    return GetModel().ReadIntSuffix(suf);
  }
  ArrayRef<double> ReadSuffix(const SuffixDef<double>& suf) override {
    return GetModel().ReadDblSuffix(suf);
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
  void SetSolutionFileName(const std::string& fileName) override {
    if (HaveSolH())
      GetSolH().OverrideSolutionFileName(fileName);
  }
  void HandleSolution(int status, fmt::CStringRef msg,
      const double *x, const double *y, double obj) override {
    if ( HaveSolH() )
      GetSolH().HandleSolution(status, msg, x, y, obj);
    else         // throw pure exception
      throw std::runtime_error(msg.c_str());
  }

  void HandleFeasibleSolution(
      int solve_code, fmt::CStringRef msg,
      const double *x, const double *y, double obj) override {
    GetSolH().HandleFeasibleSolution(solve_code, msg, x, y, obj);
  }

  const std::vector<bool>& IsVarInt() const override {
    return GetModel().IsVarInt();
  }


protected:
  /// Whether and what names
  int WantNames() const { return options_.nNames_; }


private:
  const mp::OptionValueInfo values_want_names_[4] = {
    { "0", "No names", 0},
    { "1", "(Default) Only provide names if at least one of "
           ".col / .row name files was written by AMPL "
           "(AMPL: `option [<solver>_]auxfiles rc;`) ", 1},
    { "2", "Read names from AMPL, but create generic "
           "names if not provided", 2},
    { "3", "Create generic names.", 3}
  };

  struct Options {
    int nNames_ = 1;
  };
  Options options_;

  void InitOwnOptions() {
    GetEnv().AddStoredOption("cvt:names names modelnames",
      "Whether to read or generate variable / constraint / "
      "objective names:\n"
      "\n.. value-table::\n",
      options_.nNames_, values_want_names_);
  }


private:
  std::unique_ptr<Converter> pcvt_;
  NLReadResult nl_read_result_;
  std::unique_ptr<SolutionHandler> p_sol_handler_;
};

}  // namespace mp

#endif  // MODEL_MANAGER_STD_H_
