#ifndef CONVERTERMIPTEST_H
#define CONVERTERMIPTEST_H

#include "mp/arrayref.h"

#include "mp/flat/problem_flattener.h"
#include "mp/flat/redef/MIP/converter_mip.h"
#include "mp/flat/model_api_base.h"
#include "mp/flat/nl_expr/model_api_base.h"

namespace mip_converter_test {

enum Sense {
  minimize_ = mp::obj::MIN,
  maximize_ = mp::obj::MAX
};
constexpr auto I_ = mp::var::INTEGER;
constexpr auto F_ = mp::var::CONTINUOUS;
constexpr double infty_ = 1e100;

struct MIPInstance {
  struct SparseVec {       // TODO Use LinearExpr for this
    std::vector<double> c_;
    std::vector<int>    v_;
    SparseVec(const std::vector<double>& c, const std::vector<int>& v) :
      c_(c), v_(v) { }
    SparseVec(int nnz, const double* c, const int* v) :
      c_(c, c+nnz), v_(v, v+nnz) { }
    bool operator==(const SparseVec& sv) const { return c_==sv.c_ && v_==sv.v_; }
    size_t size() const { return c_.size(); }
  };
  struct Objective {
    Sense sense_;
    SparseVec le_;        // lin expression
    bool operator==(const Objective& obj) const {
      bool eq1 = sense_==obj.sense_;
      bool eq2 = le_==obj.le_;
      if (!eq1)
        printf("Diff sense\n");
      if (!eq2) {
        printf("Diff coefs\n");
        for (size_t i=0; i<le_.size(); ++i) {
          printf(" %d:%f", le_.v_[i], le_.c_[i]);
        }
        printf("\n");
        for (size_t i=0; i<obj.le_.size(); ++i) {
          printf(" %d:%f", obj.le_.v_[i], obj.le_.c_[i]);
        }
        printf("\n");
      }
      return eq1 && eq2;
    }
  };
  using ObjectivesContainer = std::vector<Objective>;
  ObjectivesContainer objs_;
  /// Variables
  mp::VarArrayDef vars_;

  struct Constraint {
    SparseVec le_;
    double lb_;
    double ub_;
  };
  using ConstraintsContainer = std::vector<Constraint>;
  ConstraintsContainer cons_;

  bool ObjsEqual(const MIPInstance& mip) const {
    return
        objs_ == mip.objs_;
  }
  bool VarBoundsEqual(const MIPInstance& mip) const {
    return
        vars_.size() == mip.vars_.size() &&
        std::equal(vars_.plb(), vars_.plb()+vars_.size(),
                   mip.vars_.plb()) &&
        std::equal(vars_.pub(), vars_.pub()+vars_.size(),
                   mip.vars_.pub());
  }
  bool NConstrEqual(const MIPInstance& mip) const {
    return
        cons_.size() == mip.cons_.size();
  }
};

template <class Interface>
void feedInstance( Interface& interface, const MIPInstance& mip ) {
  interface.InputVariables(mip.vars_.size(),
                           mip.vars_.plb(), mip.vars_.pub(),
                           mip.vars_.ptype());
  for (MIPInstance::ObjectivesContainer::const_iterator it = mip.objs_.begin();
       it!=mip.objs_.end(); ++it) {
    interface.InputObjective((mp::obj::Type)it->sense_,
                             it->le_.size(), it->le_.c_.data(), it->le_.v_.data());
  }
  for (MIPInstance::ConstraintsContainer::const_iterator it = mip.cons_.begin();
       it!=mip.cons_.end(); ++it) {
    interface.InputAlgebraicCon(it->le_.size(), it->le_.c_.data(), it->le_.v_.data(),
                             it->lb_, it->ub_);
  }
}

/// A toy backend using struct MIPInstance
class MIPInstanceBackend :
    public mp::BasicFlatModelAPI,
    public mp::EnvKeeper
{
  MIPInstance instance_;
public:
  MIPInstanceBackend(mp::Env& e) : mp::EnvKeeper(e) { }

  static constexpr const char* GetTypeName() { return "tester"; }

  MIPInstance& GetInstance() { return instance_; }

  /// These things the concrete interface currently has to define
  void AddVariables(const mp::VarArrayDef& v) {
    instance_.vars_ = v;
  }
  void SetLinearObjective(int , const mp::LinearObjective& lo) {
    mip_converter_test::MIPInstance::SparseVec lin_part {lo.coefs(), lo.vars()};
    instance_.objs_.push_back({(Sense)lo.obj_sense(),
                              std::move(lin_part)});
  }

  /// Allow all constraint types to be compiled
  USE_BASE_CONSTRAINT_HANDLERS(mp::BasicFlatModelAPI)

  ACCEPT_CONSTRAINT(mp::LinConLE, Recommended, mp::CG_Default)

  /// Specialize for LinearConstraint
  void AddConstraint(const mp::LinConLE& lc) {
    instance_.cons_.push_back({ { (int)lc.size(), lc.pcoefs(), lc.pvars() },
                                lc.lb(), lc.ub() });
  }
};

/// Testing the default MIP interface layer
class MIPConverterTester :
    public mp::ProblemFltImpl<mp::ProblemFlattener, mp::Problem,
      mp::FlatCvtImpl<mp::MIPFlatConverter, MIPInstanceBackend> >
//    public mp::MPToMIPConverter<MIPConverterTester, MIPInstanceBackend>
{
  using Base = mp::ProblemFltImpl<mp::ProblemFlattener, mp::Problem,
    mp::FlatCvtImpl<mp::MIPFlatConverter, MIPInstanceBackend> >;
public:
  /// Construct
  MIPConverterTester(mp::Env& e) : Base(e) { InitOptions(); }
  /// This is testing API
  bool ObjsEqual(const MIPInstance& mip) {
    return GetBackend().GetInstance().ObjsEqual( mip );
  }
  bool VarBoundsEqual(const MIPInstance& mip) {
    return GetBackend().GetInstance().VarBoundsEqual( mip );
  }
  bool NConstrEqual(const MIPInstance& mip) {
    return GetBackend().GetInstance().NConstrEqual( mip );
  }
protected:
  MIPInstanceBackend& GetBackend() {
    return GetFlatCvt().GetModelAPI();
  }
};


/// Record constraints of \a Type if \a fAccept, otherwise
/// leave them for conversion
#define D2I_RECORD_CONSTR(Type, fAccept) \
public: \
  static mp::ConstraintAcceptanceLevel \
  AcceptanceLevel(const mp::Type*) { \
    return (fAccept) ? mp::ConstraintAcceptanceLevel::Recommended \
                     : mp::ConstraintAcceptanceLevel::NotAccepted; \
  } \
  static constexpr int GroupNumber(const mp::Type*) \
  { return mp::CG_Default; } \
  void AddConstraint(const mp::Type& con) \
  { cons_ ## Type ## _.push_back(con); } \
  const std::vector<mp::Type>& Get ## Type() const \
  { return cons_ ## Type ## _; } \
private: \
  std::vector<mp::Type> cons_ ## Type ## _;


/// A toy backend for disjunction absorption, auxiliary-variable,
/// bounds and piecewise-linear tests.
/// Accepts and records linear constraints, indicators, OR,
/// conditional <=, and SOS2. Tests switch native acceptance off with
/// acc:... options (see kAcc*), so that a single instantiation of the
/// converter serves all variants: each distinct backend type
/// instantiates the whole MIP converter, which makes this file
/// slow and memory-hungry to compile.
class Disj2IndBackend :
    public mp::BasicFlatModelAPI,
    public mp::EnvKeeper
{
public:
  Disj2IndBackend(mp::Env& e) : mp::EnvKeeper(e) { }

  static constexpr const char* GetTypeName() { return "d2i_tester"; }

  void AddVariables(const mp::VarArrayDef& v) { vars_ = v; }
  void SetLinearObjective(int , const mp::LinearObjective& lo)
  { obj_vars_.assign(lo.vars().begin(), lo.vars().end()); }

  /// Number of variables
  int NumVars() const { return (int)vars_.size(); }
  /// Variable bounds
  double lb(int v) const { return vars_.plb()[v]; }
  double ub(int v) const { return vars_.pub()[v]; }
  /// Submitted variable type
  mp::var::Type type(int v) const { return vars_.ptype()[v]; }
  /// Variables of the linear objective
  const std::vector<int>& ObjVars() const { return obj_vars_; }
  /// Number of submitted integer variables
  int NumIntVars() const {
    int n=0;
    for (int v=NumVars(); v--; )
      n += (mp::var::INTEGER == type(v));
    return n;
  }

  USE_BASE_CONSTRAINT_HANDLERS(mp::BasicFlatModelAPI)

  D2I_RECORD_CONSTR(LinConLE, true)
  D2I_RECORD_CONSTR(LinConGE, true)
  D2I_RECORD_CONSTR(LinConEQ, true)
  D2I_RECORD_CONSTR(LinConRange, true)
  D2I_RECORD_CONSTR(IndicatorConstraintLinLE, true)
  D2I_RECORD_CONSTR(IndicatorConstraintLinGE, true)
  D2I_RECORD_CONSTR(IndicatorConstraintLinEQ, true)
  D2I_RECORD_CONSTR(OrConstraint, true)
  D2I_RECORD_CONSTR(CondLinConLE, true)
  D2I_RECORD_CONSTR(SOS2Constraint, true)

public:
  /// Number of linear constraints of all kinds
  int NumLinCons() const {
    return int(GetLinConLE().size() + GetLinConGE().size()
               + GetLinConEQ().size() + GetLinConRange().size());
  }
  /// Number of indicators of all kinds
  int NumIndicators() const {
    return int(GetIndicatorConstraintLinLE().size()
               + GetIndicatorConstraintLinGE().size()
               + GetIndicatorConstraintLinEQ().size());
  }

private:
  mp::VarArrayDef vars_;
  std::vector<int> obj_vars_;
};


/// Acceptance options selecting what Disj2IndBackend accepts natively;
/// the rest is converted by MP.
/// Linear constraints only (indicators become big-M)
constexpr const char* kAccLinear =
    "acc:indle=0 acc:indge=0 acc:indeq=0 acc:or=0 acc:condlinle=0 "
    "acc:sos2=0";
/// Linear constraints and indicators
constexpr const char* kAccInd = "acc:or=0 acc:condlinle=0 acc:sos2=0";
/// Linear constraints, indicators and OR
constexpr const char* kAccIndOr = "acc:condlinle=0 acc:sos2=0";
/// Linear constraints, indicators, OR and conditional <=
constexpr const char* kAccIndOrCond = "acc:sos2=0";
/// Linear constraints and SOS2
constexpr const char* kAccLinearSOS2 =
    "acc:indle=0 acc:indge=0 acc:indeq=0 acc:or=0 acc:condlinle=0";


/// Interface (flattener + MIP converter + toy backend)
/// for disjunction absorption and related tests
class Disj2IndTester :
    public mp::ProblemFltImpl<mp::ProblemFlattener, mp::Problem,
      mp::FlatCvtImpl<mp::MIPFlatConverter, Disj2IndBackend> >
{
  using Base = mp::ProblemFltImpl<mp::ProblemFlattener, mp::Problem,
    mp::FlatCvtImpl<mp::MIPFlatConverter, Disj2IndBackend> >;
public:
  /// Construct with acceptance options \a acc (kAcc...)
  /// and further solver options \a opts
  Disj2IndTester(mp::Env& e, const char* acc, const char* opts="")
    : Base(e) {
    this->InitOptions();
    e.ParseOptionString((std::string(acc) + ' ' + opts).c_str());
  }

  /// The backend
  const Disj2IndBackend& Backend() {
    return this->GetFlatCvt().GetModelAPI();
  }

  /// Number of used (active or converted) items of type \a Con
  template <class Con>
  int NumUsed() {
    return this->GetFlatCvt().GetConstraintKeeper((Con*)nullptr)
        .GetNumberOfUsed();
  }

  /// Add variable
  int Var(double lb, double ub, mp::var::Type ty=mp::var::CONTINUOUS)
  { return this->AddVars(1, lb, ub, ty)[0]; }

  /// Make  (var (kind) value)
  mp::LogicalExpr Cmp(mp::expr::Kind kind, int var, double value) {
    auto& m = this->GetModel();
    return m.MakeRelational(kind, m.MakeVariable(var),
                            m.MakeNumericConstant(value));
  }

  /// Make  (a (kind) b)  for logical a, b
  mp::LogicalExpr Bin(mp::expr::Kind kind,
                      mp::LogicalExpr a, mp::LogicalExpr b) {
    return this->GetModel().MakeBinaryLogical(kind, a, b);
  }

  /// Add static logical constraint
  void Add(mp::LogicalExpr e) { this->GetModel().AddCon(e); }

  /// Add linear constraint  c0*v0 + c1*v1 (in) [lb, ub]
  void Lin2(double c0, int v0, double c1, int v1, double lb, double ub) {
    const double coefs[] = {c0, c1};
    const int vars[] = {v0, v1};
    this->InputAlgebraicCon(2, coefs, vars, lb, ub);
  }

  /// Make  (a1 (k1) r1) && (a2 (k2) r2)
  mp::LogicalExpr And2(mp::LogicalExpr a, mp::LogicalExpr b)
  { return Bin(mp::expr::AND, a, b); }

  /// Variable type as known to MP (vs the submitted one)
  mp::var::Type InternalType(int v)
  { return this->GetFlatCvt().GetModel().var_type(v); }
  /// Make  if (cond) then v1 else v2
  mp::NumericExpr If(mp::LogicalExpr cond, int v1, int v2) {
    auto& m = this->GetModel();
    return m.MakeIf(cond, m.MakeVariable(v1), m.MakeVariable(v2));
  }
  /// Make  if (cond) then c1 else c2  for constants
  mp::NumericExpr IfC(mp::LogicalExpr cond, double c1, double c2) {
    auto& m = this->GetModel();
    return m.MakeIf(cond, m.MakeNumericConstant(c1),
                    m.MakeNumericConstant(c2));
  }
  /// Make  (e (kind) value)  for a numeric expression
  mp::LogicalExpr CmpE(mp::expr::Kind kind, mp::NumericExpr e, double value) {
    auto& m = this->GetModel();
    return m.MakeRelational(kind, e, m.MakeNumericConstant(value));
  }
  /// Make  count(a, b)
  mp::NumericExpr Count2(mp::LogicalExpr a, mp::LogicalExpr b) {
    auto& m = this->GetModel();
    auto bld = m.BeginCount(2);
    bld.AddArg(a);
    bld.AddArg(b);
    return m.EndCount(bld);
  }
  /// Make  kind(var)  , e.g., floor(x)
  mp::NumericExpr Un(mp::expr::Kind kind, int var) {
    auto& m = this->GetModel();
    return m.MakeUnary(kind, m.MakeVariable(var));
  }
  /// Make  (v1 kind v2)  , e.g., x div y
  mp::NumericExpr BinN(mp::expr::Kind kind, int v1, int v2) {
    auto& m = this->GetModel();
    return m.MakeBinary(kind, m.MakeVariable(v1), m.MakeVariable(v2));
  }
  /// Make  round(var, n)
  mp::NumericExpr Round(int var, double n) {
    auto& m = this->GetModel();
    return m.MakeBinary(mp::expr::ROUND, m.MakeVariable(var),
                        m.MakeNumericConstant(n));
  }
  /// Minimize expression \a e
  void Minimize(mp::NumericExpr e)
  { this->GetModel().AddObj(mp::obj::MIN, e); }
  /// Number of auxiliary variables with MP type \a tint
  /// submitted with type \a tsub
  int NumAuxVars(mp::var::Type tint, mp::var::Type tsub) {
    int n=0;
    for (int v=this->GetFlatCvt().GetModel().num_vars(); v--; )
      if (!this->GetFlatCvt().GetModel().is_var_original(v)
          && tint == InternalType(v) && tsub == Backend().type(v))
        ++n;
    return n;
  }
};

/// Interface for piecewise-linear conversion tests.
/// Reuses Disj2IndTester's converter instantiation.
class PLTester : public Disj2IndTester {
public:
  /// Construct with acceptance options \a acc (kAcc...)
  /// and further solver options \a opts
  PLTester(mp::Env& e, const char* acc, const char* opts="")
    : Disj2IndTester(e, acc, opts) { }

  /// Add linear constraint  v == value
  void Fix(int v, double value) {
    const double coef = 1.0;
    this->InputAlgebraicCon(1, &coef, &v, value, value);
  }

  /// Make  <<0; -1, 1>> v,  i.e., abs(v)
  mp::NumericExpr AbsPL(int v) {
    auto& m = this->GetModel();
    auto bld = m.BeginPLTerm(1);
    bld.AddSlope(-1.0);
    bld.AddBreakpoint(0.0);
    bld.AddSlope(1.0);
    return m.EndPLTerm(bld, m.MakeVariable(v));
  }

  /// Maximize expression \a e
  void Maximize(mp::NumericExpr e)
  { this->GetModel().AddObj(mp::obj::MAX, e); }
};


using mp::ExpressionAcceptanceLevel;
using mp::ExprWrapper;

/// A toy expression backend: accepts linear rows, NL constraints and
/// objectives, explicifiers (NLAssign*) and exp/affine/quadratic
/// expressions; records the submitted variables and the explicified ones.
class ExprBackend :
    public mp::BasicExprModelAPI<ExprBackend, int>,
    public mp::EnvKeeper
{
public:
  using BaseModelAPI = mp::BasicExprModelAPI<ExprBackend, int>;
  ExprBackend(mp::Env& e) : mp::EnvKeeper(e) { }
  static constexpr const char* GetTypeName() { return "expr_tester"; }
  void AddVariables(const mp::VarArrayDef& v) { vars_ = v; }
  void SetLinearObjective(int , const mp::LinearObjective& ) { }
  void SetNLObjective(int , const mp::NLObjective& ) { }
  int NumVars() const { return (int)vars_.size(); }
  double lb(int v) const { return vars_.plb()[v]; }
  double ub(int v) const { return vars_.pub()[v]; }
  /// Variables assigned an explicified expression
  const std::vector<int>& Explicified() const { return explicified_; }

  USE_BASE_CONSTRAINT_HANDLERS(BaseModelAPI)
  USE_BASE_EXPRESSION_HANDLERS(BaseModelAPI)
  ACCEPT_EXPRESSION_INTERFACE(Recommended)
  Expr GetVarExpression(int i) { return i; }
  Expr GetZeroExpression() { return -1; }
  D2I_RECORD_CONSTR(LinConLE, true)
  D2I_RECORD_CONSTR(LinConGE, true)
  D2I_RECORD_CONSTR(LinConEQ, true)
  D2I_RECORD_CONSTR(LinConRange, true)
public:
  ACCEPT_CONSTRAINT(mp::NLConstraint, Recommended, mp::CG_Algebraic)
  void AddConstraint(const mp::NLConstraint& ) { }
  ACCEPT_CONSTRAINT(mp::NLAssignEQ, Recommended, mp::CG_General)
  void AddConstraint(const mp::NLAssignEQ& nle)
  { explicified_.push_back(GetVariable(nle)); }
  ACCEPT_CONSTRAINT(mp::NLAssignLE, Recommended, mp::CG_General)
  void AddConstraint(const mp::NLAssignLE& nle)
  { explicified_.push_back(GetVariable(nle)); }
  ACCEPT_CONSTRAINT(mp::NLAssignGE, Recommended, mp::CG_General)
  void AddConstraint(const mp::NLAssignGE& nle)
  { explicified_.push_back(GetVariable(nle)); }
  ACCEPT_EXPRESSION(mp::NLAffineExpression, Recommended)
  Expr AddExpression(const mp::NLAffineExpression& ) { return -2; }
  ACCEPT_EXPRESSION(mp::NLQuadExpression, Recommended)
  Expr AddExpression(const mp::NLQuadExpression& ) { return -3; }
  ACCEPT_EXPRESSION(mp::ExpExpression, Recommended)
  Expr AddExpression(const mp::ExpExpression& ) { return -4; }
private:
  mp::VarArrayDef vars_;
  std::vector<int> explicified_;
};

/// Interface (flattener + MIP converter + expression backend)
class ExprTester :
    public mp::ProblemFltImpl<mp::ProblemFlattener, mp::Problem,
      mp::FlatCvtImpl<mp::MIPFlatConverter, ExprBackend> >
{
  using Base = mp::ProblemFltImpl<mp::ProblemFlattener, mp::Problem,
    mp::FlatCvtImpl<mp::MIPFlatConverter, ExprBackend> >;
public:
  ExprTester(mp::Env& e, const char* opts) : Base(e) {
    this->InitOptions();
    e.ParseOptionString(opts);
  }
  const ExprBackend& Backend() { return this->GetFlatCvt().GetModelAPI(); }
  int Var(double lb, double ub, mp::var::Type ty=mp::var::CONTINUOUS)
  { return this->AddVars(1, lb, ub, ty)[0]; }
  mp::var::Type InternalType(int v)
  { return this->GetFlatCvt().GetModel().var_type(v); }
  mp::LogicalExpr Cmp(mp::expr::Kind kind, int var, double value) {
    auto& m = this->GetModel();
    return m.MakeRelational(kind, m.MakeVariable(var),
                            m.MakeNumericConstant(value));
  }
  mp::LogicalExpr CmpE(mp::expr::Kind kind, mp::NumericExpr e, double value) {
    auto& m = this->GetModel();
    return m.MakeRelational(kind, e, m.MakeNumericConstant(value));
  }
  void Add(mp::LogicalExpr e) { this->GetModel().AddCon(e); }
  void Minimize(mp::NumericExpr e)
  { this->GetModel().AddObj(mp::obj::MIN, e); }
};

} // namespace mip_converter_test


#endif // CONVERTERMIPTEST_H
