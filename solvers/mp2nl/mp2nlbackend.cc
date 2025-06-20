#include <vector>
#include <climits>
#include <cfloat>
#include <filesystem>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "mp2nlbackend.h"

extern "C" {
  #include "mp2nl-ampls-c-api.h"    // MP2NL AMPLS C API
}
#include "mp/ampls-cpp-api.h"


std::unique_ptr<mp::BasicBackend> CreateMP2NLBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::MP2NLBackend()};
}


namespace mp {

/// Create MP2NL Model Manager
/// @param gc: the MP2NL common handle
/// @param e: environment
/// @param pre: presolver to be returned,
/// need it to convert solution data
/// @return MP2NLModelMgr
std::unique_ptr<BasicModelManager>
CreateMP2NLModelMgr(MP2NLCommon&, Env&, pre::BasicValuePresolver*&);


// We don't provde bigM, it's user's risk.
// 1e5 is usually the maximum numerically stable value
MP2NLBackend::ConfigMap MP2NLBackend::config_map_ {
    { "baron",
     "acc:alldiff=0 acc:numberofconst=0 acc:numberofvar=0 "
     "acc:count=0 "
     "acc:indle=0 acc:indge=0 acc:indeq=0 acc:not=0 "
     "acc:condlinle=0 acc:condlineq=0 acc:condlinge=0 "
     "acc:condlinlt=0 acc:condlingt=0 acc:ifthen=0 "
     "acc:and=0 acc:or=0 acc:abs=0 acc:max=0 acc:min=0 "
     "cvt:multoutcard=0 "
     "acc:sin=0 acc:cos=0 acc:tan=0 acc:asin=0 acc:acos=0 acc:atan=0 "
     "acc:sinh=0 acc:cosh=0 acc:tanh=0 acc:asinh=0 acc:acosh=0 acc:atanh=0 "
     // "acc:log=0 "  // Actually PLApprox performs better
     "acc:sos1=0 acc:sos2=0 acc:compl=0 acc:impl=0" },
    { "baronmp",              // Don't need MP2NL actually
     "cvt:multoutcard=0 " },
    { "cbc", "cvt:multoutcard=0 " },
    { "gurobi", "cvt:multoutcard=0 " },
    { "gcg", "cvt:multoutcard=0 " },
    { "highs", "cvt:multoutcard=0 " },
    { "cplex", "cvt:multoutcard=0 " },
    { "mosek", "cvt:multoutcard=0 " },
    { "copt", "cvt:multoutcard=0 " },
    { "xpress", "cvt:multoutcard=0 " },
    { "scip", "cvt:multoutcard=0 " },
    { "gurobiasl",
     "acc:count=0 "
     "acc:alldiff=0 acc:numberofconst=0 acc:numberofvar=0 "
     "acc:indle=0 acc:indge=0 acc:indeq=0 acc:not=0 "
     "acc:condlinle=0 acc:condlineq=0 acc:condlinge=0 "
     "acc:condlinlt=0 acc:condlingt=0 acc:ifthen=0 "
     "acc:and=0 acc:or=0 acc:abs=0 acc:max=0 acc:min=0 "
     "cvt:multoutcard=0 cvt:socp=1 "
     "acc:pow=0 acc:powconstexp=0 acc:div=0 "
     "acc:log=0 acc:logA=0 acc:exp=0 acc:expA=0 "
     "acc:sin=0 acc:cos=0 acc:tan=0 acc:asin=0 acc:acos=0 acc:atan=0 "
     "acc:sinh=0 acc:cosh=0 acc:tanh=0 acc:asinh=0 acc:acosh=0 acc:atanh=0 "
     "acc:sos1=0 acc:sos2=0 acc:compl=0 acc:impl=0" },
    { "cplexasl",
     "acc:count=0 "
     "acc:alldiff=0 acc:numberofconst=0 acc:numberofvar=0 "
     "acc:indle=0 acc:indge=0 acc:indeq=0 acc:not=0 "
     "acc:condlinle=0 acc:condlineq=0 acc:condlinge=0 "
     "acc:condlinlt=0 acc:condlingt=0 acc:ifthen=0 "
     "acc:and=0 acc:or=0 acc:abs=0 acc:max=0 acc:min=0 "
     "cvt:multoutcard=0 cvt:socp=1 cvt:prod=7 "
     "acc:pow=0 acc:powconstexp=0 acc:div=0 "
     "acc:log=0 acc:logA=0 acc:exp=0 acc:expA=0 "
     "acc:sin=0 acc:cos=0 acc:tan=0 acc:asin=0 acc:acos=0 acc:atan=0 "
     "acc:sinh=0 acc:cosh=0 acc:tanh=0 acc:asinh=0 acc:acosh=0 acc:atanh=0 "
     "acc:sos1=0 acc:sos2=0 acc:compl=0 acc:impl=0" },
    { "xpressasl",
     "acc:count=0 "
     "acc:alldiff=0 acc:numberofconst=0 acc:numberofvar=0 "
     "acc:indle=0 acc:indge=0 acc:indeq=0 acc:not=0 "
     "acc:condlinle=0 acc:condlineq=0 acc:condlinge=0 "
     "acc:condlinlt=0 acc:condlingt=0 acc:ifthen=0 "
     "acc:and=0 acc:or=0 acc:abs=0 acc:max=0 acc:min=0 "
     "cvt:multoutcard=0 cvt:socp=1 "
     "acc:pow=0 acc:powconstexp=0 acc:div=0 "
     "acc:log=0 acc:logA=0 acc:exp=0 acc:expA=0 "
     "acc:sin=0 acc:cos=0 acc:tan=0 acc:asin=0 acc:acos=0 acc:atan=0 "
     "acc:sinh=0 acc:cosh=0 acc:tanh=0 acc:asinh=0 acc:acosh=0 acc:atanh=0 "
     "acc:sos1=0 acc:sos2=0 acc:compl=0 acc:impl=0" },
    { "knitro",
     "acc:count=0 "
     "acc:alldiff=0 acc:numberofconst=0 acc:numberofvar=0 "
     "acc:indle=0 acc:indge=0 acc:indeq=0 acc:not=0 "
     "acc:condlinle=0 acc:condlineq=0 acc:condlinge=0 "
     "acc:condlinlt=0 acc:condlingt=0 acc:ifthen=0 "
     "acc:and=0 acc:or=0 acc:abs=0 acc:max=0 acc:min=0 "
     "cvt:multoutcard=0 "
     // minmaxabs: indeed better linearize as Bob said
     "acc:sos1=0 acc:sos2=0 acc:compl=0 acc:impl=0" },  // @todo compl
    { "lindoglobal",
     "acc:count=0 "
     "acc:alldiff=0 acc:numberofconst=0 acc:numberofvar=0 "
     "acc:indle=0 acc:indge=0 acc:indeq=0 acc:not=0 "
     "acc:condlinle=0 acc:condlineq=0 acc:condlinge=0 "
     "acc:condlinlt=0 acc:condlingt=0 acc:ifthen=0 "
     "acc:and=0 acc:or=0 "   // acc:abs=0 acc:max=0 acc:min=0 "
     "cvt:multoutcard=0 "
     "acc:sos1=0 acc:sos2=0 acc:compl=0 acc:impl=0" },   // @todo compl
    { "couenne",
     "acc:count=0 "
     "acc:alldiff=0 acc:numberofconst=0 acc:numberofvar=0 "
     "acc:indle=0 acc:indge=0 acc:indeq=0 acc:not=0 "
     "acc:condlinle=0 acc:condlineq=0 acc:condlinge=0 "
     "acc:condlinlt=0 acc:condlingt=0 acc:ifthen=0 "
     "acc:and=0 acc:or=0 acc:abs=0 acc:max=0 acc:min=0 "
     "cvt:multoutcard=0 cvt:socp=1 "
     "acc:asin=0 acc:acos=0 acc:atan=0 "
     "acc:asinh=0 acc:acosh=0 acc:atanh=0 "
     "acc:sos1=0 acc:sos2=0 acc:compl=0 acc:impl=0" },  // @todo compl
    { "bonmin",
     "acc:count=0 "
     "acc:alldiff=0 acc:numberofconst=0 acc:numberofvar=0 "
     "acc:indle=0 acc:indge=0 acc:indeq=0 acc:not=0 "
     "acc:condlinle=0 acc:condlineq=0 acc:condlinge=0 "
     "acc:condlinlt=0 acc:condlingt=0 acc:ifthen=0 "
     "acc:and=0 acc:or=0 acc:abs=0 acc:max=0 acc:min=0 "
     "cvt:multoutcard=0 cvt:socp=1 "
     "acc:sos1=0 acc:sos2=0 acc:compl=0 acc:impl=0" },  // @todo compl
    { "conopt",
     "acc:count=0 "
     "acc:alldiff=0 acc:numberofconst=0 acc:numberofvar=0 "
     "acc:indle=0 acc:indge=0 acc:indeq=0 acc:not=0 "
     "acc:condlinle=0 acc:condlineq=0 acc:condlinge=0 "
     "acc:condlinlt=0 acc:condlingt=0 acc:ifthen=0 "
     "acc:and=0 acc:or=0 acc:abs=0 acc:max=0 acc:min=0 "
     "cvt:multoutcard=0 cvt:socp=1 "
     "acc:sos1=0 acc:sos2=0 acc:compl=0 acc:impl=0" },  // @todo compl
    { "conopt4",
     "acc:count=0 "
     "acc:alldiff=0 acc:numberofconst=0 acc:numberofvar=0 "
     "acc:indle=0 acc:indge=0 acc:indeq=0 acc:not=0 "
     "acc:condlinle=0 acc:condlineq=0 acc:condlinge=0 "
     "acc:condlinlt=0 acc:condlingt=0 acc:ifthen=0 "
     "acc:and=0 acc:or=0 acc:abs=0 acc:max=0 acc:min=0 "
     "cvt:multoutcard=0 cvt:socp=1 "
     "acc:sos1=0 acc:sos2=0 acc:compl=0 acc:impl=0" },  // @todo compl
    { "snopt",
     "acc:count=0 "
     "acc:alldiff=0 acc:numberofconst=0 acc:numberofvar=0 "
     "acc:indle=0 acc:indge=0 acc:indeq=0 acc:not=0 "
     "acc:condlinle=0 acc:condlineq=0 acc:condlinge=0 "
     "acc:condlinlt=0 acc:condlingt=0 acc:ifthen=0 "
     "acc:and=0 acc:or=0 acc:abs=0 acc:max=0 acc:min=0 "
     "cvt:multoutcard=0 cvt:socp=1 "
     // minmaxabs: rely on ASL, or better linearize?
     "acc:sos1=0 acc:sos2=0 acc:compl=0 acc:impl=0" },  // @todo compl
    { "ipopt",
     "acc:count=0 "
     "acc:alldiff=0 acc:numberofconst=0 acc:numberofvar=0 "
     "acc:indle=0 acc:indge=0 acc:indeq=0 acc:not=0 "
     "acc:condlinle=0 acc:condlineq=0 acc:condlinge=0 "
     "acc:condlinlt=0 acc:condlingt=0 acc:ifthen=0 "
     "acc:and=0 acc:or=0 acc:abs=0 acc:max=0 acc:min=0 "
     "cvt:multoutcard=0 cvt:socp=1 "
     // minmaxabs: indeed better linearize as Bob said?
     "acc:sos1=0 acc:sos2=0 acc:compl=0 acc:impl=0" },  // @todo compl
    { "minos",
     "acc:count=0 "
     "acc:alldiff=0 acc:numberofconst=0 acc:numberofvar=0 "
     "acc:indle=0 acc:indge=0 acc:indeq=0 acc:not=0 "
     "acc:condlinle=0 acc:condlineq=0 acc:condlinge=0 "
     "acc:condlinlt=0 acc:condlingt=0 acc:ifthen=0 "
     "acc:and=0 acc:or=0 acc:abs=0 acc:max=0 acc:min=0 "
     "cvt:multoutcard=0 cvt:socp=1 "
     // minmaxabs: indeed better linearize as Bob said?
     "acc:sos1=0 acc:sos2=0 acc:compl=0 acc:impl=0" },  // @todo compl
    { "lgo",
     "acc:count=0 "
     "acc:alldiff=0 acc:numberofconst=0 acc:numberofvar=0 "
     "acc:indle=0 acc:indge=0 acc:indeq=0 acc:not=0 "
     "acc:condlinle=0 acc:condlineq=0 acc:condlinge=0 "
     "acc:condlinlt=0 acc:condlingt=0 acc:ifthen=0 "
     "acc:and=0 acc:or=0 acc:abs=0 acc:max=0 acc:min=0 "
     "cvt:multoutcard=0 cvt:socp=1 "
     "acc:sos1=0 acc:sos2=0 acc:compl=0 acc:impl=0" },  // @todo compl
    { "loqo",
     "acc:count=0 "
     "acc:alldiff=0 acc:numberofconst=0 acc:numberofvar=0 "
     "acc:indle=0 acc:indge=0 acc:indeq=0 acc:not=0 "
     "acc:condlinle=0 acc:condlineq=0 acc:condlinge=0 "
     "acc:condlinlt=0 acc:condlingt=0 acc:ifthen=0 "
     "acc:and=0 acc:or=0 acc:abs=0 acc:max=0 acc:min=0 "
     "cvt:multoutcard=0 cvt:socp=1 "
     "acc:sos1=0 acc:sos2=0 acc:compl=0 acc:impl=0" },  // @todo compl
    { "raposa",
     "acc:count=0 "
     "acc:alldiff=0 acc:numberofconst=0 acc:numberofvar=0 "
     "acc:indle=0 acc:indge=0 acc:indeq=0 acc:not=0 "
     "acc:condlinle=0 acc:condlineq=0 acc:condlinge=0 "
     "acc:condlinlt=0 acc:condlingt=0 acc:ifthen=0 "
     "acc:and=0 acc:or=0 acc:abs=0 acc:max=0 acc:min=0 "
     "cvt:multoutcard=0 cvt:socp=1 "
     "acc:sos1=0 acc:sos2=0 acc:compl=0 acc:impl=0" },  // @todo compl
    { "octeract",
     "acc:count=0 "
     "acc:alldiff=0 acc:numberofconst=0 acc:numberofvar=0 "
     "acc:indle=0 acc:indge=0 acc:indeq=0 acc:not=0 "
     "acc:condlinle=0 acc:condlineq=0 acc:condlinge=0 "
     "acc:condlinlt=0 acc:condlingt=0 acc:ifthen=0 "
     "acc:and=0 acc:or=0 acc:abs=0 acc:max=0 acc:min=0 "
     "cvt:multoutcard=0 cvt:socp=1 "
     "acc:sos1=0 acc:sos2=0 acc:compl=0 acc:impl=0" },  // @todo compl
    { "ilogcp",
     "acc:compl=0 cvt:multoutcard=0 " },
    { "gecode",
     "acc:compl=0 cvt:multoutcard=0 " },
};


/// Implement MP2NLSolverQueryCallbacks
class MP2NLSolverQueryCallbacksImpl
    : public MP2NLSolverQueryCallbacks {
public:
  /// Construct
  MP2NLSolverQueryCallbacksImpl(MP2NLBackend& be) : be_(be) { }


  /// Get initial X
  ArrayRef< std::pair<int, double> > GetInitialGuesses() override {
    return be_.GetInitialGuess();
  }

  /// Get initial Y
  ArrayRef<double> GetInitialDualGuesses() override {
    return be_.GetInitialDualGuess();
  }

  /// Suffix names.
  std::set<std::string> GetSuffixNames() override {
    std::set<std::string> suf_skip {
      "sosno", "ref", "sos", "sosref"          // AMPL SOS constraints
    };
    auto result = be_.GetSuffixNames();
    for (const auto& k2: suf_skip)
      result.erase(k2);
    return result;
  }

  /// Get model suffix with given name
  MP2NLModelSuffix GetModelSuffix(
      const std::string& name) override {
    MP2NLModelSuffix result;
    assert(internal::NUM_SUFFIX_KINDS == result.values_.size());
    result.name_ = name;
    int ni=0, nd=0;
    for (int kind=0; kind<internal::NUM_SUFFIX_KINDS; ++kind) {
      // ReadDblSuffix also reads int if that's the case
      int fint=0;
      result.values_[kind] = be_.ReadDblSuffix({name, kind}, &fint);
      ni += fint;
      nd += (result.values_[kind].size() != 0);
    }
    if (nd>ni)
      result.flags_ |= suf::FLOAT;            // it's really real-valued
    else {
      assert(nd==ni);
    }
    // Treating all suffixes generically,
    // with range constraints this should work for .status
    auto suf_pre = be_.GetValuePresolver().PresolveGenericDbl({
        result.values_[0],
        result.values_[1],
        result.values_[2]
    });
    result.values_[0] = suf_pre.GetVarValues()();
    result.values_[1] = suf_pre.GetConValues()(CG_Algebraic);
    const auto& suf_log = suf_pre.GetConValues()(CG_Logical);
    result.values_[1].insert(result.values_[1].end(),
                                    suf_log.begin(), suf_log.end());
    if (be_.multiobj())   // Skip because we compress objectives
      result.values_[2].clear(); // @todo native multiobj
    else
      result.values_[2] = suf_pre.GetObjValues()();
    return result;
  }
private:
  MP2NLBackend& be_;
};


MP2NLBackend::MP2NLBackend() {
  OpenSolver();

  /// Create a ModelManager
  pre::BasicValuePresolver* pPre;
  auto data = CreateMP2NLModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);

  p_nlsi_ = get_other()->p_nlsi_;    // before copying
  /// Copy env/lp to ModelAPI
  copy_common_info_to_other();

  p_qc_ = std::make_unique<MP2NLSolverQueryCallbacksImpl>(*this);
  p_nlsi_->ProvideQueryCallbacks(*p_qc_.get());
}

MP2NLBackend::~MP2NLBackend() {
  CloseSolver();
}


const char* MP2NLBackend::GetBackendName()
  { return "MP2NLBackend"; }

std::string MP2NLBackend::GetSolverVersion() {
  return "0.1";
}

std::string MP2NLBackend::set_external_libs() {
  return "NLWriter2";
}

ArrayRef<double> MP2NLBackend::PrimalSolution() {
  return GetNLSolver().GetX();
}

pre::ValueMapDbl MP2NLBackend::DualSolution() {
  return {{ { CG_Algebraic, DualSolution_LP() } }};
}

ArrayRef<double> MP2NLBackend::DualSolution_LP() {
  return GetNLSolver().GetY();
}

double MP2NLBackend::ObjectiveValue() const {
  return 0.0;                        // SOL does not provide one
}

double MP2NLBackend::NodeCount() const {
  return 0.0;
}

double MP2NLBackend::SimplexIterations() const {
  return 0.0;
}

int MP2NLBackend::BarrierIterations() const {
  return 0;
}


void MP2NLBackend::SetInterrupter(mp::Interrupter *inter) {
  // TODO
}

void MP2NLBackend::Solve() {
  GetNLSolver().Solve(
      storedOptions_.solver_.c_str(),
      storedOptions_.solver_options_.c_str());
  WindupMP2NLSolve();
}

void MP2NLBackend::WindupMP2NLSolve() { }

void MP2NLBackend::ReportResults() {
  ReportMP2NLResults();
  BaseBackend::ReportResults();
}

void MP2NLBackend::ReportSuffixes() {
  auto sufnames = GetNLSolver().GetSuffixNames();
  // We pass all suffixes currently.
  // To control them, now the subsolver's options should be used.
  // @todo adapt them for each target subsolver (basis, etc.)
  for (const auto& sufname: sufnames) {
    ReportModelSuffix( GetNLSolver().GetModelSuffix(sufname) );
  }
}

void MP2NLBackend::ReportModelSuffix(
    const MP2NLModelSuffix& modelsuf) {
  assert(internal::NUM_SUFFIX_KINDS == modelsuf.values_.size());
  auto n_alg_cons = GetNLSolver().GetNumAlgCons();
  assert(modelsuf.values_[1].empty()
         || n_alg_cons <= modelsuf.values_[1].size());
  // Treating all suffixes generically,
  // with range constraints this should work for .status.
  // @todo adapted handling of basis, etc.
  auto suf_post = GetValuePresolver().PostsolveGenericDbl(
      {
          modelsuf.values_[0],
          modelsuf.values_[1].empty() ? pre::ValueMapDbl{} :
              pre::ValueMapDbl{{
                  { CG_Algebraic,
                   { modelsuf.values_[1].begin(), modelsuf.values_[1].begin()+n_alg_cons } },
                  { CG_Logical,
                   { modelsuf.values_[1].begin()+n_alg_cons, modelsuf.values_[1].end() } }
              }},
          modelsuf.values_[2]
      });
  if (modelsuf.values_[0].size())           // Only if presolved values were supplied
    ReportSuffix(modelsuf.name_,
                 0 | modelsuf.flags_, suf_post.GetVarValues()(),
                 modelsuf.table_);
  if (modelsuf.values_[1].size())
    ReportSuffix(modelsuf.name_,
                 1 | modelsuf.flags_, suf_post.GetConValues()(),
                 modelsuf.table_);
  if (modelsuf.values_[2].size())
    ReportSuffix(modelsuf.name_,
                 2 | modelsuf.flags_, suf_post.GetObjValues()(),
                 modelsuf.table_);
  if (modelsuf.values_[3].size())
    ReportSuffix(modelsuf.name_,
                 3 | modelsuf.flags_, modelsuf.values_[3],
                 modelsuf.table_);
}

void MP2NLBackend::ReportMP2NLResults() {
  SetStatus( GetSolveResult() );
  AddMP2NLMessages();
  if (need_multiple_solutions())
    ReportMP2NLPool();
}
std::vector<double> MP2NLBackend::getPoolSolution(int i)
{
  int num_vars = NumVars();
  std::vector<double> vars(num_vars);
  return vars;
}
double MP2NLBackend::getPoolObjective(int i)
{
  double obj {0.0};
  return obj;
}
void MP2NLBackend::ReportMP2NLPool() {
  if (!IsMIP())
    return;
  int iPoolSolution = -1;
  int nsolutions = 0;
  
  while (++iPoolSolution < nsolutions) {
    ReportIntermediateSolution(
      { getPoolSolution(iPoolSolution),
        {}, { getPoolObjective(iPoolSolution) } });
  }
}


void MP2NLBackend::AddMP2NLMessages() {
}

std::pair<int, std::string> MP2NLBackend::GetSolveResult() {
  return {
          GetNLSolver().GetSolveResult(),
      GetNLSolver().GetSolveMessage() };
}


void MP2NLBackend::FinishOptionParsing() {
	if (storedOptions_.config_attempted_.size()
			&& !storedOptions_.config_.size())
		AddWarning("SolverConfig",
							 fmt::format("MP2NL: configuration '{}', \n"
													 "assumed for solver '{}',\n"
													 "is unknown. Defaults applied\n"
													 "(use nl:printconfigs and/or nl:config)",
													 storedOptions_.config_attempted_,
													 storedOptions_.solver_) );
}


////////////////////////////// OPTIONS /////////////////////////////////

static const mp::OptionValueInfo childsel[] = {
  { "d", "down", 0},
  { "u", "up", 1},
  { "p", "pseudo costs", 2},
  { "i", "inference", 3},
  { "l", "lp value", 4},
  { "r", "root LP value difference", 5},
  { "h", "hybrid inference/root LP value difference (default)", 6}
};


void MP2NLBackend::InitCustomOptions() {

  set_option_header(
    "MP2NL Optimizer Options for AMPL\n"
    "--------------------------------------------\n"
    "\n"
		"To configure MP2NL, assign a string specifying option values to the "
		"AMPL option ``mp2nl_options``, as well as the subsolver option values "
		"to its own AMPL option. For example::\n"
    "\n"
    "  ampl: option mp2nl_options 'solver=baron';\n"
    "  ampl: option baron_options 'outlev=1 iisfind=1';"
      );

	AddIntOption("tech:outlev outlev",
    "0*/1: Verbosity for the MP2NL driver. "
                  "For the underlying solver, use the <subsolver>_options "
                  "environment variable.",
		&MP2NLBackend::GetOutlev, &MP2NLBackend::SetOutlev);

  // AddStoredOption("tech:logfile logfile",
  //                 "Log file name.",
  //                 storedOptions_.logFile_);

	AddStrOption("nl:solver solver nlsolver",
									"Subsolver (underlying AMPL solver.)\n"
							 "\n"
							 "Automatically sets solver configuration, if known.",
							 &MP2NLBackend::GetSolver,
									 &MP2NLBackend::SetSolver);

	AddFlagOption(
		"nl:printconfigs printconfigs",
		"Print the configurations for various values of the "
		"solver/config options.",
		&MP2NLBackend::GetDummyFlagOption,
				&MP2NLBackend::SetConfigPrintFlag);

	AddStrOption(
		"nl:config config",
		"Choose subsolver configuration (see nl:printconfigs).",
		&MP2NLBackend::GetConfig,
				&MP2NLBackend::SetConfig);

	AddStoredOption("nl:solver_options solver_options slv_opts",
                  "Subsolver options.\n\n"
                  "This way is for convenience; the preferred and "
                  "dominating way is to use the <subsolver>_options "
                  "environment variable (which can be set in AMPL "
                  "in the usual way as 'option <subsolver>_options \"...\";')",
                  storedOptions_.solver_options_);

  /// Enforce time limit?
  AddStoredOption("lim:time timelim timelimit time_limit",
    "Enforce limit on solve time (in seconds; default: 1e+20). DOES NOTHING CURRENTLY.",
    storedOptions_.tilim_);


  /// Custom infinity value?
  // AddSolverOption("num:infinity infinity",
  //   "Values larger than this are considered infinity (default: 1e+20)",
  //   "numerics/infinity", 1e+10, DBL_MAX);

}


double MP2NLBackend::MIPGap() {
  return 0.0;
}
double MP2NLBackend::BestDualBound() {
  return 0.0;
}

double MP2NLBackend::MIPGapAbs() {
  return 0.0;
}


ArrayRef<int> MP2NLBackend::VarStatii() {
  
  std::vector<int> vars(NumVars());
  /*
  MP2NL_GetBasis(lp(), vars.data(), NULL);
  for (auto& s : vars) {
    switch (s) {
    case MP2NL_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case MP2NL_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case MP2NL_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case MP2NL_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case MP2NL_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown MP2NL VBasis value: {}", s));
    }
  }
  */
  return vars;
}

ArrayRef<int> MP2NLBackend::ConStatii() {

  std::vector<int> cons(NumLinCons());
  /*
  MP2NL_GetBasis(lp(), NULL, cons.data());
  for (auto& s : cons) {
    switch (s) {
    case MP2NL_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case MP2NL_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case MP2NL_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case MP2NL_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case MP2NL_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown MP2NL VBasis value: {}", s));
    }
  }*/
  return cons;
}

void MP2NLBackend::VarStatii(ArrayRef<int> vst) {
  int index[1];
  std::vector<int> stt(vst.data(), vst.data() + vst.size());
  /*
  for (auto j = stt.size(); j--; ) {
    auto& s = stt[j];
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = MP2NL_BASIS_BASIC;
      break;
    case BasicStatus::low:
      s = MP2NL_BASIS_LOWER;
      break;
    case BasicStatus::equ:
      s = MP2NL_BASIS_FIXED;
      break;
    case BasicStatus::upp:
      s = MP2NL_BASIS_UPPER;
      break;
    case BasicStatus::sup:
    case BasicStatus::btw:
      s = MP2NL_BASIS_SUPERBASIC;
      break;
    case BasicStatus::none:
      /// 'none' is assigned to new variables. Compute low/upp/sup:
      /// Depending on where 0.0 is between bounds
      double lb, ub;
      index[0] = (int)j;
      if(!MP2NL_GetColInfo(lp(), MP2NL_DBLINFO_LB, 1, index, &lb) &&
        !MP2NL_GetColInfo(lp(), MP2NL_DBLINFO_UB, 1, index, &ub))
      { 
        if (lb >= -1e-6)
          s = -1;
        else if (ub <= 1e-6)
          s = -2;
        else
          s = -3;  // or, leave at 0?
      }
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL var status value: {}", s));
    }
  }
  MP2NL_SetBasis(lp(), stt.data(), NULL);
  */
}

void MP2NLBackend::ConStatii(ArrayRef<int> cst) {
  /*
  std::vector<int> stt(cst.data(), cst.data() + cst.size());
  for (auto& s : stt) {
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = MP2NL_BASIS_BASIC;
      break;
    case BasicStatus::none:   // for 'none', which is the status
    case BasicStatus::upp:    // assigned to new rows, it seems good to guess
    case BasicStatus::sup:    // a valid status.
    case BasicStatus::low:    // 
    case BasicStatus::equ:    // For active constraints, it is usually 'sup'.
    case BasicStatus::btw:    // We could compute slack to decide though.
      s = MP2NL_BASIS_SUPERBASIC;
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL con status value: {}", s));
    }
  }
  MP2NL_SetBasis(lp(), NULL, stt.data());
  */
}

SolutionBasis MP2NLBackend::GetBasis() {
  std::vector<int> varstt = VarStatii();
  std::vector<int> constt = ConStatii();
  if (varstt.size() && constt.size()) {
    auto mv = GetValuePresolver().PostsolveBasis(
      { std::move(varstt),
        {{{ CG_Linear, std::move(constt) }}} });
    varstt = mv.GetVarValues()();
    constt = mv.GetConValues()();
    assert(varstt.size());
  }
  return { std::move(varstt), std::move(constt) };
}

void MP2NLBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetValuePresolver().PresolveBasis(
    { basis.varstt, basis.constt });
  auto varstt = mv.GetVarValues()();
  auto constt = mv.GetConValues()(CG_Linear);
  assert(varstt.size());
  assert(constt.size());
  VarStatii(varstt);
  ConStatii(constt);
}

void MP2NLBackend::AddPrimalDualStart(Solution sol)
{
  auto mv = GetValuePresolver().PresolveSolution(
        { sol.primal, sol.dual } );
  auto x0 = mv.GetVarValues()();
  y0_ = mv.GetConValues()(CG_Algebraic);
  auto ms = GetValuePresolver().PresolveGenericInt(
      { sol.spars_primal } );
  auto s0 = ms.GetVarValues()();
  x0_.clear();
  x0_.reserve(x0.size());
  for (int i=0; i<(int)x0.size(); ++i) {
    if (s0[i]) {
      x0_.push_back( {i, x0[i]} );
    }
  }
}

void MP2NLBackend::AddMIPStart(
    ArrayRef<double> x0_unpres, ArrayRef<int> s0_unpres) {
  auto mv = GetValuePresolver().PresolveSolution( { x0_unpres } );
  auto ms = GetValuePresolver().PresolveGenericInt( { s0_unpres } );
  auto x0 = mv.GetVarValues()();
  auto s0 = ms.GetVarValues()();
  x0_.clear();
  x0_.reserve(x0.size());
  for (int i=0; i<(int)x0.size(); ++i) {
    if (s0[i]) {
      x0_.push_back( {i, x0[i]} );
    }
  }
}


void MP2NLBackend::SetOutlev(const SolverOption& , int val) {
	storedOptions_.outlev_ = val>0;
	set_verbose_mode(val>0);
}

void MP2NLBackend::SetSolver(const SolverOption& , fmt::StringRef val) {
	storedOptions_.solver_ = val;
	auto cfg = ExtractSolverConfigName(val);
	if (DoSetConfig(cfg)) {
		storedOptions_.config_ = cfg;
		if (verbose_mode())
			this->Print("Set solver configuration '{}'\n", cfg);
	} else
		storedOptions_.config_attempted_ = cfg;
}

std::string MP2NLBackend::ExtractSolverConfigName(
		fmt::StringRef solver){
	std::filesystem::path p {solver.to_string()};
	if (p.has_stem())
		return p.stem().string();
	return "";
}

void MP2NLBackend::SetConfigPrintFlag(const SolverOption& , bool ) {
	fmt::print("Solver configurations (options solver, config):\n");
	fmt::print("-----------------------------------------------\n");
	for (const auto& [key, val]: config_map_) {
		fmt::MemoryWriter writer;
		internal::FormatRST(writer, val, 8);
		this->Print("\"{}\":\n\t{}\n", key, writer.c_str());
	}
	fmt::print("-----------------------------------------------\n\n");
}

void MP2NLBackend::SetConfig(const SolverOption& , fmt::StringRef val) {
	MP_ASSERT_ALWAYS( DoSetConfig(val),
										fmt::format("MP2NL: solver configuration '{}' unknown\n"
																"(see nl:printconfigs)", val) );
	storedOptions_.config_ = val;
}

bool MP2NLBackend::DoSetConfig(fmt::StringRef val) {
  auto it = config_map_.find(val);
  if (config_map_.end() == it)
    return false;
  ParseOptionString(it->second.c_str(), NO_OPTION_ECHO);
  return true;
}

void MP2NLBackend::DoWriteProblem(const std::string& name) {
}

} // namespace mp


// AMPLs
AMPLS_MP_Solver* Open_MP2NL(CCallbacks cb = {}) {
  return AMPLS__internal__Open(std::unique_ptr<mp::BasicBackend>{new mp::MP2NLBackend()},
    cb);
}

void AMPLSClose_MP2NL(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* AMPLSGetModel_MP2NL(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::MP2NLBackend*>(AMPLSGetBackend(slv));
}
