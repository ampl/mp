#include "mosekmodelapi.h"


namespace mp {

void MosekModelAPI::InitProblemModificationPhase(const FlatModelInfo* info) {
	/// Preallocate algebraic constraints.
	/// MOSEK 10 seems to handle all algebraic constraints as 1 group.
	/// CG_Algebraic, etc. are the constraint group indexes
  /// provided in ACCEPT_CONSTRAINT macros.
  MOSEK_CCALL(MSK_appendcons(lp(),
														 info->GetNumberOfConstraintsOfGroup(CG_Algebraic)));
}
std::string& myreplace(std::string& s, const std::string& from, const std::string& to)
{
  for (size_t pos = 0; (pos = s.find(from, pos)) != std::string::npos; pos += to.size())
    s.replace(pos, from.size(), to);
  return s;
}

std::string sanitizeName(std::string n) {
  // Mosek does not like square brackets or spaces in variable names
  std::replace(n.begin(), n.end(), '[', '(');
  std::replace(n.begin(), n.end(), ']', ')');
  std::replace(n.begin(), n.end(), ' ', '_');
  return n;
}
void MosekModelAPI::AddVariables(const VarArrayDef& v) {
  MOSEK_CCALL(MSK_appendvars(lp(), v.size()));
  for (size_t i = v.size(); i--; ) {
    MSKboundkeye k;
    bool freelb = v.plb()[i] == -std::numeric_limits<double>::infinity();
    bool freeub = v.pub()[i] == std::numeric_limits<double>::infinity();
    if (freelb && freeub)
      k = MSK_BK_FR;
    else if (freelb == freeub)
    {
      if (v.plb()[i] == v.pub()[i])
        k = MSK_BK_FX;
      else
        k = MSK_BK_RA;
    } else {
      if (freelb)
        k = MSK_BK_UP;
      else
        k = MSK_BK_LO;
    }
    MOSEK_CCALL(MSK_putvarbound(lp(), i, k, v.plb()[i], v.pub()[i]));
    MOSEK_CCALL(MSK_putvartype(lp(), i,
      var::Type::CONTINUOUS == v.ptype()[i] ? MSK_VAR_TYPE_CONT :
      MSK_VAR_TYPE_INT));
   }
  if (v.pnames()) 
    for (size_t i = v.size(); i--; )
      MOSEK_CCALL(MSK_putvarname(lp(), i, sanitizeName(v.pnames()[i]).c_str()));
}

void MosekModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (iobj < 1) {
    MOSEK_CCALL(MSK_putobjsense(lp(),
      obj::Type::MAX == lo.obj_sense() ? MSK_OBJECTIVE_SENSE_MAXIMIZE : MSK_OBJECTIVE_SENSE_MINIMIZE));
    for (int i = 0; i < lo.num_terms(); i++) {
      MOSEK_CCALL(MSK_putcj(lp(), lo.vars()[i], lo.coefs()[i]));
    }
  }
  else {
    MP_RAISE("Multiple objectives not supported");
  }
}

/// Add QP terms to objective
void AddObjQuadraticPart(MSKtask_t lp,
                         const QuadTerms& qt) {
  auto coefs = qt.coefs();
  auto vars1 = qt.vars1();
  auto vars2 = qt.vars2();
  for (auto i=coefs.size(); i--; ) {
    if (vars1[i]==vars2[i])
      coefs[i] *= 2;                 // MOSEK represents lower submatrix
    else if (vars1[i] < vars2[i])
      std::swap(vars1[i], vars2[i]);
  }
  MOSEK_CCALL( MSK_putqobj(lp,
                           coefs.size(),
                           vars1.data(),
                           vars2.data(),
                           coefs.data()) );
}

void MosekModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  if (1 > iobj) {
    SetLinearObjective(iobj, qo);                         // add the linear part
    const auto& qt = qo.GetQPTerms();
    AddObjQuadraticPart(lp(), qt);
  }
  else {
    MP_RAISE("Multiple objectives not supported");
  }
}

void MosekModelAPI::AddLinearConstraint(
    MSKtask_t lp, size_t size, MSKboundkey_enum key,
    double lb, double ub,
    const int* vindex, const double* values, const char* name) {
  /* Linear + quadratic constraints are preallocated in
     InitProblemModificationPhase() */
  MOSEK_CCALL(MSK_putarow(lp, n_alg_cons_, size, vindex, values));
  MOSEK_CCALL(MSK_putconbound(lp, n_alg_cons_, key, lb, ub));
  MOSEK_CCALL(MSK_putconname(lp, n_alg_cons_, name));

  ++n_alg_cons_;
}

// Some linear constraints are processed as Range, so we set the bound keys here.
// Sometimes the separate AddConstraint functions are also called (e.g. logical constraints), see below.
void MosekModelAPI::AddConstraint(const LinConRange& lc)
{
  double lb = lc.lb();
  double ub = lc.ub();
  MSKboundkey_enum key;

  if (lb == -std::numeric_limits<double>::infinity())
  {
    if (ub == std::numeric_limits<double>::infinity())
      key = MSK_BK_FR;
    else
      key = MSK_BK_UP;
  }
  else
  {
    if (ub == std::numeric_limits<double>::infinity())
      key = MSK_BK_LO;
    else if (lb == ub)
      key = MSK_BK_FX;
    else
      key = MSK_BK_RA;
  }

  AddLinearConstraint(lp(), lc.size(), key, lb, ub, lc.pvars(), lc.pcoefs(), lc.name());
}

void MosekModelAPI::AddConstraint(const LinConLE& lc)
{
  AddLinearConstraint(lp(), lc.size(), MSK_BK_UP, lc.lb(), lc.ub(), lc.pvars(), lc.pcoefs(), lc.name());
}
void MosekModelAPI::AddConstraint(const LinConEQ& lc)
{
  AddLinearConstraint(lp(), lc.size(), MSK_BK_FX, lc.lb(), lc.ub(), lc.pvars(), lc.pcoefs(), lc.name());
}
void MosekModelAPI::AddConstraint(const LinConGE& lc)
{
  AddLinearConstraint(lp(), lc.size(), MSK_BK_LO, lc.lb(), lc.ub(), lc.pvars(), lc.pcoefs(), lc.name());
}

/// Add QP terms to constraint
void AddConQuadraticPart(MSKtask_t lp,
                         int i_con,
                         const QuadTerms& qt) {
  auto coefs = qt.coefs();
  auto vars1 = qt.vars1();
  auto vars2 = qt.vars2();
  for (auto i=coefs.size(); i--; ) {
    if (vars1[i]==vars2[i])
      coefs[i] *= 2;                 // MOSEK represents lower submatrix
    else if (vars1[i] < vars2[i])
      std::swap(vars1[i], vars2[i]);
  }
  MOSEK_CCALL( MSK_putqconk(lp,
                           i_con,
                           coefs.size(),
                           vars1.data(),
                           vars2.data(),
                           coefs.data()) );
}

void MosekModelAPI::AddConstraint(const QuadConRange& qc) {
  double lb = qc.lb();
  double ub = qc.ub();
  MSKboundkey_enum key;

  if (lb == -std::numeric_limits<double>::infinity())
  {
    if (ub == std::numeric_limits<double>::infinity())
      key = MSK_BK_FR;
    else
      key = MSK_BK_UP;
  }
  else
  {
    if (ub == std::numeric_limits<double>::infinity())
      key = MSK_BK_LO;
    else if (lb == ub)
      key = MSK_BK_FX;
    else
      key = MSK_BK_RA;
  }

  const auto& lt = qc.GetLinTerms();
  AddLinearConstraint(lp(), lt.size(), key, lb, ub, lt.pvars(), lt.pcoefs(), qc.name());
  AddConQuadraticPart(lp(), n_alg_cons_-1, qc.GetQPTerms());
}

void MosekModelAPI::AddConstraint( const QuadConLE& qc ) {
  const auto& lt = qc.GetLinTerms();
  AddLinearConstraint(lp(), lt.size(), MSK_BK_UP, qc.lb(), qc.ub(),
                      lt.pvars(), lt.pcoefs(), qc.name());
  AddConQuadraticPart(lp(), n_alg_cons_-1, qc.GetQPTerms());
}

void MosekModelAPI::AddConstraint( const QuadConEQ& qc ) {
  const auto& lt = qc.GetLinTerms();
  AddLinearConstraint(lp(), lt.size(), MSK_BK_FX, qc.lb(), qc.ub(),
                      lt.pvars(), lt.pcoefs(), qc.name());
  AddConQuadraticPart(lp(), n_alg_cons_-1, qc.GetQPTerms());
}

void MosekModelAPI::AddConstraint( const QuadConGE& qc ) {
  const auto& lt = qc.GetLinTerms();
  AddLinearConstraint(lp(), lt.size(), MSK_BK_LO, qc.lb(), qc.ub(),
                      lt.pvars(), lt.pcoefs(), qc.name());
  AddConQuadraticPart(lp(), n_alg_cons_-1, qc.GetQPTerms());
}

void MosekModelAPI::AddConstraint(
		const QuadraticConeConstraint& qc) {
	MSKint64t numafe_prev=0;
	MOSEK_CCALL( MSK_getnumafe(lp(), &numafe_prev) );
	auto nnz = qc.GetArguments().size();
	// Is it too slow to add cones 1 by 1?
	MOSEK_CCALL( MSK_appendafes(lp(), nnz) );
	std::vector<MSKint64t> afeidx(nnz);
	for (auto idx=nnz; idx--; )        // Fill new AFE indexes
		afeidx[idx] = numafe_prev+idx;
	MOSEK_CCALL(
				MSK_putafefentrylist(lp(), nnz,
														 afeidx.data(),
														 qc.GetArguments().data(),  // assumes it's int32
														 qc.GetParameters().data()) );
	MSKint64t domidx=-1;
	MOSEK_CCALL( MSK_appendquadraticconedomain(lp(), nnz, &domidx) );
	MOSEK_CCALL( MSK_appendaccseq(lp(), domidx, nnz, numafe_prev, NULL) );
}

void MosekModelAPI::AddConstraint(
		const RotatedQuadraticConeConstraint& qc) {
	MSKint64t numafe_prev=0;
	MOSEK_CCALL( MSK_getnumafe(lp(), &numafe_prev) );
	auto nnz = qc.GetArguments().size();
	// Is it too slow to add cones 1 by 1?
	MOSEK_CCALL( MSK_appendafes(lp(), nnz) );
	std::vector<MSKint64t> afeidx(nnz);
	for (auto idx=nnz; idx--; )        // Fill new AFE indexes
		afeidx[idx] = numafe_prev+idx;
	MOSEK_CCALL(
				MSK_putafefentrylist(lp(), nnz,
														 afeidx.data(),
														 qc.GetArguments().data(),  // assumes it's int32
														 qc.GetParameters().data()) );
	MSKint64t domidx=-1;
	MOSEK_CCALL( MSK_appendrquadraticconedomain(lp(), nnz, &domidx) );
	MOSEK_CCALL( MSK_appendaccseq(lp(), domidx, nnz, numafe_prev, NULL) );
}

/// Add indicator as disjunction
/// (binvar!=binval) \/ (lincon).
/// 1st afe: binvar+binval-1=0
/// 2nd afe: (conbody)-(conrhs) (<=>) 0
template <int sens>
void AddIndicator(MSKtask_t lp,
             const IndicatorConstraint< LinConRhs<sens> >& ic) {
  const auto& lincon = ic.get_constraint();
  /////////////////////////////// 1. AFEs /////////////////////////
  MSKint64t numafe_prev=0;
  MOSEK_CCALL( MSK_getnumafe(lp, &numafe_prev) );
  // Is it too slow to add indicators 1 by 1?
  MOSEK_CCALL( MSK_appendafes(lp, 2) );
  const MSKint32t fvaridx[] = {ic.get_binary_var()};
  const MSKrealt  fval[]    = {1.0};
  const MSKrealt  g[]       = {ic.get_binary_value() - 1.0,
                               -lincon.rhs()};
  MOSEK_CCALL(                          // 1st afe
        MSK_putafefentrylist(lp, 1,
                             &numafe_prev, fvaridx, fval) );
  MOSEK_CCALL(
        MSK_putafeg(lp, numafe_prev, g[0]) );
  auto nnz = lincon.size();
  std::vector<MSKint64t> afeidx(nnz, numafe_prev+1);
  MOSEK_CCALL(                          // 2nd afe
        MSK_putafefentrylist(lp, nnz,
                             afeidx.data(),
                             lincon.pvars(),      // assumes it's int32
                             lincon.pcoefs()) );
  MOSEK_CCALL(
        MSK_putafeg(lp, numafe_prev+1, g[1]) );
  //////////////////////////////// 2. Domains /////////////////////
  MSKint64t dom1=-1, dom2=-1;
  MOSEK_CCALL( MSK_appendrzerodomain(lp, 1, &dom1) );
  if (sens<0)
    MOSEK_CCALL( MSK_appendrminusdomain(lp, 1, &dom2) );
  else if (sens>0)
    MOSEK_CCALL( MSK_appendrplusdomain(lp, 1, &dom2) );
  else
    MOSEK_CCALL( MSK_appendrzerodomain(lp, 1, &dom2) );
  //////////////////////////////// 3. DJC /////////////////////////
  MSKint64t numdjcs_prev=0;
  MOSEK_CCALL( MSK_getnumdjc(lp, &numdjcs_prev) );
  MOSEK_CCALL( MSK_appenddjcs(lp, 1) );
  const MSKint64t domidxlist[] = {dom1, dom2};
  const MSKint64t afeidxlist[] = {numafe_prev, numafe_prev+1};
  const MSKint64t termsizelist[] = {1, 1};

  MOSEK_CCALL( MSK_putdjc(lp,
                 numdjcs_prev,           // DJC index
                 2, domidxlist,
                 2, afeidxlist,
                 NULL,                   // Unused
                 2, termsizelist) );
}

void MosekModelAPI::AddConstraint(const IndicatorConstraintLinLE &ic) {
  AddIndicator(lp(), ic);
}

void MosekModelAPI::AddConstraint(const IndicatorConstraintLinEQ &ic) {
  AddIndicator(lp(), ic);
}

void MosekModelAPI::AddConstraint(const IndicatorConstraintLinGE &ic) {
  AddIndicator(lp(), ic);
}


void MosekModelAPI::AddConstraint(const SOS1Constraint& sos) {
  fmt::print("Adding SOS1 constraint {}\n", sos.GetTypeName());
  // TODO
/*  int type = MOSEK_SOS_TYPE1;
  int beg = 0;
  const int size = sos.size();
  MOSEK_CCALL(MOSEK_AddSOSs(lp(), 1, &type, &beg,
    &size, (int*)sos.get_vars().data(),
    (double*)sos.get_weights().data())); */
}

void MosekModelAPI::AddConstraint(const SOS2Constraint& sos) {
  fmt::print("Adding SOS1 constraint {}\n", sos.GetTypeName());
  // TODO
  /*int type = MOSEK_SOS_TYPE2;
  int beg = 0;
  const int size = sos.size();
  MOSEK_CCALL(MOSEK_AddSOSs(lp(), 1, &type, &beg,
    &size, (int*)sos.get_vars().data(),
    (double*)sos.get_weights().data()));*/
}


void MosekModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp
