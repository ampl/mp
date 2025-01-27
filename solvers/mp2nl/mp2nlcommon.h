#ifndef MP2NLCOMMON_H
#define MP2NLCOMMON_H

#include <array>
#include <vector>
#include <set>
#include <string>

#include "mp/arrayref.h"
#include "mp/backend-to-model-api.h"

#include "mp/format.h"


namespace mp {

/// A (MP2NL)ModelSuffix has
/// all suffixes with given name,
/// for all kinds provided.
/// The reason we send them together is that
/// FlatConverter pre-/postsolves them together,
/// see FlatBackend::ReadModelSuffix().
///
/// We store all suffixes as double's,
/// but transmit as int's if originally so.
struct MP2NLModelSuffix {
  std::string name_;
  int flags_ {};       // suf::FLOAT for float-valued
  std::string table_;
  /// Values for each kind.
  std::array< std::vector<double>, 4 > values_;
};


/// Interface for query callbacks
/// from ModelAPI into Backend
class MP2NLSolverQueryCallbacks {
public:
  /// Virtual destruct
  virtual ~MP2NLSolverQueryCallbacks() { }

  /// Get initial X
  virtual ArrayRef< std::pair<int, double> > GetInitialGuesses() = 0;

  /// Get initial Y
  virtual ArrayRef<double> GetInitialDualGuesses() = 0;

  /// Suffix names.
  virtual std::set<std::string> GetSuffixNames() = 0;

  /// Get model suffix with given name
  virtual MP2NLModelSuffix GetModelSuffix(
      const std::string& name) = 0;
};


/// MP2NLSolverNLParams
struct MP2NLSolverNLParams {
  /// NL format text mode?
  int if_nl_text_ {};
  /// In text mode, comments?
  int if_nl_comments_ {};
  /// Custom stub file?
  const char* stub_file_ {};
};


/// Interface for NLSolver
/// to be used from Backend
class MP2NLSolverIntf {
public:
  /// Virtual destruct
  virtual ~MP2NLSolverIntf() { }

  /// Provide query callback (e.g., for suffixes)
  void ProvideQueryCallbacks(MP2NLSolverQueryCallbacks& qc)
  { p_nlsq_ = &qc; }

  void ProvideNLParams(MP2NLSolverNLParams prm) { prm_ = prm; }

  /// Solve
  virtual void Solve(const char* solver, const char* sopts) = 0;

  /// AMPL solve result code
  virtual int GetSolveResult() const = 0;

  /// Solve result message or error message
  virtual const char* GetSolveMessage() const = 0;

  /// Number of backspaces to print
  /// if printing the solve message right here,
  /// or skip so many symbols first.
  virtual int GetSolveMessageNbs() const = 0;

  /// Stub file used
  virtual const char* GetFileStub() const = 0;

  /// Objno used.
  /// @note AMPL .sol file does not provide
  ///   objective values, only in the solve message if at all.
  virtual int GetObjnoUsed() const = 0;

  /// Primal solution
  virtual ArrayRef<double> GetX() const = 0;

  /// Dual solution
  virtual ArrayRef<double> GetY() const = 0;

  /// Suffix names
  virtual std::set<std::string> GetSuffixNames() = 0;

  /// Get model suffix with given name
  virtual const MP2NLModelSuffix& GetModelSuffix(
      const std::string& name) = 0;

  /// Number of algebraic cons
  virtual int GetNumAlgCons() const = 0;
public:
  MP2NLSolverQueryCallbacks* GetCallbacks()
  { assert(p_nlsq_); return p_nlsq_; }
  const MP2NLSolverNLParams& GetParams() const { return prm_; }
  MP2NLSolverNLParams& GetParams() { return prm_; }
private:
  MP2NLSolverQueryCallbacks* p_nlsq_ {};
  MP2NLSolverNLParams prm_ {};
};


/// Information shared by both
/// `MP2NLBackend` and `MP2NLModelAPI`
struct MP2NLCommonInfo {
  /// MP2NLSolver
  MP2NLSolverIntf* p_nlsi_ {};
private:
};


/// Common API for MP2NL classes
class MP2NLCommon :
		public Backend2ModelAPIConnector<MP2NLCommonInfo> {
public:
  static constexpr double Infinity() { return INFINITY; }
  static constexpr double MinusInfinity() { return -INFINITY; }

protected:
  void OpenSolver();
  void CloseSolver();

  const MP2NLSolverIntf& GetNLSolver() const
  { assert(p_nlsi_); return *p_nlsi_; }
  MP2NLSolverIntf& GetNLSolver()
  { assert(p_nlsi_); return *p_nlsi_; }

  int NumLinCons() const;
  int NumVars() const;
  int NumObjs() const;
  int NumQPCons() const;
  int NumSOSCons() const;
  int NumIndicatorCons() const;
};

} // namespace mp

#endif // MP2NLCOMMON_H
