/*-------------------------------------------------------------------------*/
/* AMPL/Concert driver                                       Robert Fourer */
/*                                                                         */
/* Name           : concert.cpp                                            */
/* Title          : AMPL/ILOG Concert driver                               */
/* By             : Robert Fourer                                          */
/* Date           : October 2000                                           */
/*                                                                         */
/* A driver to link AMPL linear integer programs with ILOG Concert 1.0     */
/* October 2000: Linear/Nonlinear version                                  */
/*-------------------------------------------------------------------------*/

#include "concert.h"

#include <algorithm>
#include <iostream>

#include <ilcplex/ilocplex.h>
#include <ilcp/cp.h>

#include "solvers/asl.h"
#include "solvers/nlp.h"
#include "solvers/getstub.h"
#include "solvers/opcode.hd"

using namespace std;

// for suppressing "String literal to char*" warnings
#define CSTR(s) const_cast<char*>(s)

namespace {
struct DriverOptionInfo : Option_Info {
  Driver *driver;
};

// Returns the constant term in the first objective.
real objconst0(ASL_fg *a) {
  expr *e = a->I.obj_de_->e;
  return reinterpret_cast<size_t>(e->op) == OPNUM ?
      reinterpret_cast<expr_n*>(e)->v : 0;
}

char *skip_space(char *s) {
  while (*s && !isspace(*s))
    ++s;
  return s;
}

struct CPOptionInfo {
  IloCP::IntParam param;
  int start;
};

const CPOptionInfo LogVerbosity = {IloCP::LogVerbosity, IloCP::Quiet};
const CPOptionInfo DefaultInferenceLevel =
  {IloCP::DefaultInferenceLevel, IloCP::Default};

class CPLEXOptimizer : public Optimizer {
 private:
  IloCplex cplex_;

 public:
  CPLEXOptimizer(IloEnv env) : cplex_(env) {}

  IloCplex cplex() const { return cplex_; }
  IloAlgorithm algorithm() const { return cplex_; }

  void set_option(const void *key, int value);
};

void CPLEXOptimizer::set_option(const void *key, int value) {
  cplex_.setParam(
      static_cast<IloCplex::IntParam>(reinterpret_cast<size_t>(key)), value);
}

class CPOptimizer : public Optimizer {
 private:
  IloSolver solver_;

 public:
  CPOptimizer(IloEnv env) : solver_(env) {}

  IloSolver solver() const { return solver_; }
  IloAlgorithm algorithm() const { return solver_; }

  void set_option(const void *key, int value);
};

void CPOptimizer::set_option(const void *key, int value) {
  const CPOptionInfo *info = static_cast<const CPOptionInfo*>(key);
  solver_.setParameter(info->param, info->start + value);
}
}

Optimizer::~Optimizer() {}

keyword Driver::keywords_[] = { /* must be in alphabetical order */
   KW(CSTR("debugexpr"), Driver::set_int_option, Driver::DEBUGEXPR,
      CSTR("print debugging information for expression trees")),
   KW(CSTR("defaultinferencelevel"), Driver::set_cp_int_option,
      &DefaultInferenceLevel, CSTR("default inference level for constraints")),
   KW(CSTR("ilogcplex"), Driver::use_cplex, Driver::ILOGOPTTYPE,
      CSTR("use ILOG CPLEX optimizer")),
   KW(CSTR("ilogsolver"), Driver::use_cpoptimizer, Driver::ILOGOPTTYPE,
      CSTR("use ILOG Constraint Programming optimizer")),
   KW(CSTR("logverbosity"), Driver::set_cp_int_option,
      &LogVerbosity, CSTR("verbosity of the search log")),
   KW(CSTR("timing"), Driver::set_bool_option, Driver::TIMING,
      CSTR("display timings for the run")),
   KW(CSTR("usenumberof"), Driver::set_bool_option, Driver::USENUMBEROF,
      CSTR("consolidate 'numberof' expressions"))
};

Driver::Driver() :
   mod_(env_), asl(reinterpret_cast<ASL_fg*>(ASL_alloc(ASL_read_fg))),
   gotopttype(false), n_badvals(0) {
   options_[DEBUGEXPR] = 0;
   options_[ILOGOPTTYPE] = DEFAULT_OPT;
   options_[TIMING] = 0;
   options_[USENUMBEROF] = 1;

   version_.resize(strlen(IloConcertVersion::_ILO_NAME) + 100);
   snprintf(&version_[0], version_.size() - 1,
       "%s %d.%d.%d", IloConcertVersion::_ILO_NAME,
       IloConcertVersion::_ILO_MAJOR_VERSION,
       IloConcertVersion::_ILO_MINOR_VERSION,
       IloConcertVersion::_ILO_TECH_VERSION);
   DriverOptionInfo *doi = 0;
   oinfo_.reset(doi = new DriverOptionInfo());
   oinfo_->sname = CSTR("concert");
   oinfo_->bsname = &version_[0];
   oinfo_->opname = CSTR("concert_options");
   oinfo_->keywds = keywords_;
   oinfo_->n_keywds = sizeof(keywords_) / sizeof(*keywords_);
   oinfo_->version = &version_[0];
   oinfo_->driver_date = 20120521;
   doi->driver = this;
}

Driver::~Driver() {
   env_.end();
}

char *Driver::use_cplex(Option_Info *oi, keyword *, char *value) {
   Driver *d = static_cast<DriverOptionInfo*>(oi)->driver;
   if (!d->gotopttype)
      d->options_[ILOGOPTTYPE] = CPLEX;
   return value;
}

char *Driver::use_cpoptimizer(Option_Info *oi, keyword *, char *value) {
   Driver *d = static_cast<DriverOptionInfo*>(oi)->driver;
   if (!d->gotopttype)
      d->options_[ILOGOPTTYPE] = CPOPTIMIZER;
   return value;
}

char *Driver::set_int_option(Option_Info *oi, keyword *kw, char *value) {
   Driver *d = static_cast<DriverOptionInfo*>(oi)->driver;
   if (!d->gotopttype)
      return skip_space(value);
   keyword thiskw(*kw);
   thiskw.info = d->options_ + reinterpret_cast<size_t>(kw->info);
   return I_val(oi, &thiskw, value);
}

char *Driver::set_bool_option(Option_Info *oi, keyword *kw, char *value) {
   Driver *d = static_cast<DriverOptionInfo*>(oi)->driver;
   if (!d->gotopttype)
      return skip_space(value);
   keyword thiskw(*kw);
   int intval = 0;
   thiskw.info = &intval;
   char *result = I_val(oi, &thiskw, value);
   if (intval != 0 && intval != 1) {
     ++d->n_badvals;
     cerr << "Invalid value " << value
        << " for directive " << kw->name << endl;
   } else d->options_[reinterpret_cast<size_t>(kw->info)] = intval;
   return result;
}

char *Driver::set_cp_int_option(Option_Info *oi, keyword *kw, char *value) {
   Driver *d = static_cast<DriverOptionInfo*>(oi)->driver;
   if (!d->gotopttype)
      return skip_space(value);
   keyword thiskw(*kw);
   int intval = 0;
   thiskw.info = &intval;
   char *result = I_val(oi, &thiskw, value);
   d->optimizer_->set_option(kw->info, intval);
   // TODO: check
   return result;
}

bool Driver::parse_options(char **argv) {
   // Get optimizer type.
   gotopttype = false;
   if (getopts(argv, oinfo_.get()))
      return false;

   int& ilogopttype = options_[ILOGOPTTYPE];
   if (ilogopttype == DEFAULT_OPT)
      ilogopttype = nlo + nlc + n_lcon == 0 ? CPLEX : CPOPTIMIZER;
   if (ilogopttype == CPLEX)
      optimizer_.reset(new CPLEXOptimizer(env_));
   else optimizer_.reset(new CPOptimizer(env_));

   // Parse remaining options.
   gotopttype = true;
   n_badvals = 0;
   if (getopts(argv, oinfo_.get()) || n_badvals != 0)
      return false;
   return true;
}

/*----------------------------------------------------------------------

  Main Program

----------------------------------------------------------------------*/

int Driver::run(char **argv) {
   /*** Initialize timers ***/

   IloTimer timer(env_);
   timer.start();

   IloNum Times[4];
   Times[0] = timer.getTime();

   /*** Get name of .nl file; read problem sizes ***/

   char *stub = getstub(&argv, oinfo_.get());
   if (!stub)
     usage_ASL(oinfo_.get(), 1);
   FILE *nl = jac0dim(stub, strlen(stub));

   /*** Read coefficients & bounds & expression tree from .nl file ***/

   Uvx = static_cast<real*>(Malloc(n_var * sizeof(real)));
   Urhsx = static_cast<real*>(Malloc(n_con * sizeof(real)));

   efunc *r_ops_int[N_OPS];
   for (int i = 0; i < N_OPS; i++)
      r_ops_int[i] = reinterpret_cast<efunc*>(i);
   asl->I.r_ops_ = r_ops_int;
   want_derivs = 0;
   fg_read(nl, ASL_allow_CLP);
   asl->I.r_ops_ = 0;

   if (!parse_options(argv))
      return 1;

   /*-------------------------------------------------------------------

     Set up optimization problem in ILOG Concert

   -------------------------------------------------------------------*/

   vars_ = IloNumVarArray(env_,n_var);

   int n_var_int = nbv + niv + nlvbi + nlvci + nlvoi;
   for (int j = 0; j < n_var - n_var_int; j++)
      vars_[j] = IloNumVar(env_, LUv[j], Uvx[j], ILOFLOAT);
   for (int j = n_var - n_var_int; j < n_var; j++)
      vars_[j] = IloNumVar(env_, LUv[j], Uvx[j], ILOINT);

   IloObjective MinOrMax(env_);

   if (n_obj > 0) {
      IloExpr objExpr(env_, objconst0(asl));
      if (0 < nlo)
         objExpr += build_expr (obj_de[0].e);
      for (ograd *og = Ograd[0]; og; og = og->next)
         objExpr += (og -> coef) * vars_[og -> varno];
      MinOrMax = IloObjective (env_, objExpr,
         objtype[0] == 0 ? IloObjective::Minimize : IloObjective::Maximize);
      IloAdd (mod_, MinOrMax);
   }

   IloRangeArray Con(env_,n_con);

   for (int i = 0; i < n_con; i++) {
      IloExpr conExpr(env_);
      for (cgrad *cg = Cgrad[i]; cg; cg = cg->next)
         conExpr += (cg -> coef) * vars_[cg -> varno];
      if (i < nlc) 
         conExpr += build_expr (con_de[i].e);
      Con[i] = (LUrhs[i] <= conExpr <= Urhsx[i]);
   }

   IloConstraintArray LCon(env_,n_lcon);

   for (int i = 0; i < n_lcon; i++)
      LCon[i] = build_constr (lcon_de[i].e);

   if (n_con > 0) mod_.add (Con);
   if (n_lcon > 0) mod_.add (LCon);

   finish_building_numberof ();

   int timing = get_option(TIMING);
   if (timing) Times[1] = timer.getTime();

   /*-------------------------------------------------------------------

     Solve integer/linear program in CPLEX

   -------------------------------------------------------------------*/

   IloAlgorithm alg(optimizer_->algorithm());
   alg.extract (mod_);
   if (timing) Times[2] = timer.getTime();
   IloBool successful = alg.solve();
   if (timing) Times[3] = timer.getTime();

   CPLEXOptimizer *cplex_opt = dynamic_cast<CPLEXOptimizer*>(optimizer_.get());
   if (cplex_opt) {
      IloCplex cplex(cplex_opt->cplex());
      IloNum objValue = cplex.getObjValue();

      int sSoFar = 0;
      char sMsg[256];

      sSoFar += Sprintf(sMsg,
         "\n%s: optimal solution found\n", oinfo_->bsname);

      if (cplex.isMIP()) {
         sSoFar += g_fmtop(sMsg+sSoFar,cplex.getNnodes());
         sSoFar += Sprintf(sMsg+sSoFar, " nodes, ");
         sSoFar += g_fmtop(sMsg+sSoFar,cplex.getNiterations());
         sSoFar += Sprintf(sMsg+sSoFar, " iterations, objective ");
         g_fmtop(sMsg+sSoFar, objValue);

         vector<real> Xopt(n_var);
         for (int j = 0; j < n_var; j++) Xopt[j] = cplex.getValue(vars_[j]);
         write_sol(sMsg, &Xopt[0], 0, oinfo_.get());
      }
      else {
         sSoFar += g_fmtop(sMsg+sSoFar,cplex.getNiterations());
         sSoFar += Sprintf(sMsg+sSoFar, " iterations, objective ");
         g_fmtop(sMsg+sSoFar, objValue);

         vector<real> Xopt(n_var);
         vector<real> Piopt(n_con);
         for (int j = 0; j < n_var; j++) Xopt[j] = cplex.getValue(vars_[j]);
         for (int i = 0; i < n_con; i++) Piopt[i] = cplex.getDual(Con[i]);
         write_sol(sMsg, &Xopt[0], &Piopt[0], oinfo_.get());
      }
   }

   /*-------------------------------------------------------------------

     Solve problem in ILOG Solver

   -------------------------------------------------------------------*/

   else {
      IloSolver solver(dynamic_cast<CPOptimizer&>(*optimizer_).solver());

      int sSoFar = 0;
      char sMsg[256];

      if (successful) {
         sSoFar += Sprintf(sMsg,
            "\n%s: solution found\n", oinfo_->bsname);

         sSoFar += g_fmtop(sMsg+sSoFar,solver.getNumberOfChoicePoints());
         sSoFar += Sprintf(sMsg+sSoFar, " choice points, ");
         sSoFar += g_fmtop(sMsg+sSoFar,solver.getNumberOfFails());
         sSoFar += Sprintf(sMsg+sSoFar, " fails");
         if (n_obj > 0) {
            sSoFar += Sprintf(sMsg+sSoFar, ", objective ");
            g_fmtop(sMsg+sSoFar, solver.getValue(MinOrMax));
         }

         real *Xopt = new real [n_var];

         for (int j = 0; j < n_var; j++) Xopt[j] = solver.getValue(vars_[j]);
         write_sol(sMsg, Xopt, 0, oinfo_.get());
         delete [] Xopt;
      }
      else {
         sSoFar += Sprintf(sMsg,
            "\n%s: no solution found!\n", oinfo_->bsname);
         write_sol(sMsg, 0, 0, oinfo_.get());
      }
   }

   if (timing) {
      Times[4] = timer.getTime();
      cerr << endl
           << "Define = " << Times[1] - Times[0] << endl
           << "Setup =  " << Times[2] - Times[1] << endl
           << "Solve =  " << Times[3] - Times[2] << endl
           << "Output = " << Times[4] - Times[3] << endl;
   }
   return 0;
}
