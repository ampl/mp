// ===============================================================
// Generic Ampl interface to facilitate usage from other languages
// Dominique Orban
// Vancouver, April 2014.
// Montreal, February 2015.
// ===============================================================

#include "aslinterface.h"

// Module functions.

void *asl_init(char *stub) {
  ASL_pfgh *asl = (ASL_pfgh*)ASL_alloc(ASL_read_pfgh);
  if (!asl) return NULL;

  FILE *ampl_file = jac0dim(stub, (fint)strlen(stub));

  // Allocate room to store problem data
  if (! (asl->i.X0_    = (real *)M1alloc(asl->i.n_var_ * sizeof(real)))) return NULL;
  if (! (asl->i.LUv_   = (real *)M1alloc(asl->i.n_var_ * sizeof(real)))) return NULL;
  if (! (asl->i.Uvx_   = (real *)M1alloc(asl->i.n_var_ * sizeof(real)))) return NULL;
  if (! (asl->i.pi0_   = (real *)M1alloc(asl->i.n_con_ * sizeof(real)))) return NULL;
  if (! (asl->i.LUrhs_ = (real *)M1alloc(asl->i.n_con_ * sizeof(real)))) return NULL;
  if (! (asl->i.Urhsx_ = (real* )M1alloc(asl->i.n_con_ * sizeof(real)))) return NULL;

  // Read in ASL structure
  asl->i.want_xpi0_ = 3;        // Read primal and dual estimates
  pfgh_read(ampl_file , 0);     // pfgh_read closes the file.

  return (void *)asl;
}

void asl_write_sol(void *asl, const char *msg, double *x, double *y) {
  ASL *this_asl = (ASL *)asl;
  write_sol_ASL(asl, msg, x, y, 0); // Do not handle Option_Info for now.
  return;
}

void asl_finalize(void *asl) {
  ASL *this_asl = (ASL *)asl;
  ASL_free((ASL **)(&this_asl));
  return;
}

// Problem setup.

int asl_objtype(void *asl) {
  return ((ASL *)asl)->i.objtype_[0];  // 0 means minimization problem.
}

int asl_nvar(void *asl) {
  return ((ASL *)asl)->i.n_var_;
}

int asl_ncon(void *asl) {
  return ((ASL *)asl)->i.n_con_;
}

int asl_nlc(void *asl) {
  return ((ASL *)asl)->i.nlc_;
}

int asl_nlnc(void *asl) {
  return ((ASL *)asl)->i.nlnc_;
}

int asl_nnzj(void *asl) {
  return ((ASL *)asl)->i.nzc_;
}

int asl_nnzh(void *asl) {
  ASL *this_asl = (ASL *)asl;
  return (int)((*(this_asl->p.Sphset))(this_asl, 0, -1, 1, 1, 1));
}

int asl_islp(void *asl) {
  ASL *this_asl = (ASL *)asl;
  return ((this_asl->i.nlo_ + this_asl->i.nlc_ + this_asl->i.nlnc_) > 0 ? 0 : 1);
}

double *asl_x0(void *asl) {
  return ((ASL *)asl)->i.X0_;
}

double *asl_y0(void *asl) {
  return ((ASL *)asl)->i.pi0_;
}

double *asl_lvar(void *asl) {
  return ((ASL *)asl)->i.LUv_;
}

double *asl_uvar(void *asl) {
  return ((ASL *)asl)->i.Uvx_;
}

double *asl_lcon(void *asl) {
  return ((ASL *)asl)->i.LUrhs_;
}

double *asl_ucon(void *asl) {
  return ((ASL *)asl)->i.Urhsx_;
}

// Objective.

void asl_varscale(void *asl, double *s) {
  fint ne;
  int this_nvar = ((ASL *)asl)->i.n_var_;

  for (int i = 0; i < this_nvar; i++)
    varscale_ASL((ASL *)asl, i, s[i], &ne);
  return;
}

double asl_obj(void *asl, double *x) {
  fint ne;
  return (*((ASL *)asl)->p.Objval)((ASL *)asl, 0, x, &ne);
}

void asl_grad(void *asl, double *x, double *g) {
  fint ne;

  (*((ASL *)asl)->p.Objgrd)((ASL *)asl, 0, x, g, &ne);
}

// Lagrangian.

void asl_lagscale(void *asl, double s) {
  fint ne;
  lagscale_ASL((ASL *)asl, s, &ne);
  return;
}

// Constraints and Jacobian.

void asl_conscale(void *asl, double *s) {
  fint ne;
  int this_ncon = ((ASL *)asl)->i.n_con_;

  for (int j = 0; j < this_ncon; j++)
    conscale_ASL((ASL *)asl, j, s[j], &ne);
  return;
}

void asl_cons(void *asl, double *x, double *c) {
  fint ne;

  (*((ASL *)asl)->p.Conval)((ASL *)asl, x, c, &ne);
}

double asl_jcon(void *asl, double *x, int j) {
  fint ne;
  return (*((ASL *)asl)->p.Conival)((ASL *)asl, j, x, &ne);
}

void asl_jcongrad(void *asl, double *x, double *g, int j) {
  ASL *this_asl = (ASL *)asl;
  fint ne;
  (*(this_asl->p.Congrd))(this_asl, j, x, g, &ne);
}

size_t asl_sparse_congrad_nnz(void *asl, int j) {
  ASL *this_asl = (ASL *)asl;
  size_t nzgj = 0;
  cgrad *cg;
  for (cg = this_asl->i.Cgrad_[j]; cg; cg = cg->next) nzgj++;

  return nzgj;
}

void asl_sparse_congrad(void *asl, double *x, int j, int64_t *inds, double *vals) {
  ASL *this_asl = (ASL *)asl;
  cgrad *cg;


  int congrd_mode_bkup = this_asl->i.congrd_mode;
  this_asl->i.congrd_mode = 1;  // Sparse gradient mode.

  fint ne;
  (*(this_asl->p.Congrd))(this_asl, j, x, vals, &ne);

  int k = 0;
  for (cg = this_asl->i.Cgrad_[j]; cg; cg = cg->next)
      inds[k++] = (long)(cg->varno) + 1;

  this_asl->i.congrd_mode = congrd_mode_bkup;  // Restore gradient mode.

}

// Evaluate Jacobian at x in triplet form (rows, vals, cols).
void asl_jac(void *asl, double *x, int64_t *rows, int64_t *cols, double *vals) {
  ASL *this_asl = (ASL *)asl;
  int this_nzc = this_asl->i.nzc_, this_ncon = this_asl->i.n_con_;

  fint ne;
  (*(this_asl->p.Jacval))(this_asl, x, vals, &ne);

  // Fill in sparsity pattern. Account for 1-based indexing.
  for (int j = 0; j < this_ncon; j++)
    for (cgrad *cg = this_asl->i.Cgrad_[j]; cg; cg = cg->next) {
        rows[cg->goff] = (long)j + 1;
        cols[cg->goff] = (long)(cg->varno) + 1;
    }
}

// Hessian.

void asl_hprod(void *asl, double *y, double *v, double *hv, double w) {
  ASL *this_asl = (ASL *)asl;
  double ow[1];  // Objective weight.

  ow[0]  = this_asl->i.objtype_[0] ? -w : w;
  hvpinit_ASL(this_asl, this_asl->p.ihd_limit_, 0, NULL, y);
  (*(this_asl->p.Hvcomp))(this_asl, hv, v, -1, ow, y); // nobj=-1 so ow takes precendence.
}

void asl_hvcompd(void *asl, double *v, double *hv, int nobj) {
  ASL *this_asl = (ASL *)asl;
  (*(this_asl->p.Hvcompd))(this_asl, hv, v, nobj);
}

void asl_ghjvprod(void *asl, double *g, double *v, double *ghjv) {
  ASL *this_asl = (ASL *)asl;
  int this_ncon = this_asl->i.n_con_;
  int this_nvar = this_asl->i.n_var_;
  int this_nlc  = this_asl->i.nlc_;
  double *hv    = (double *)Malloc(this_nvar * sizeof(real));

  fint ne;
  double prod;
  int i, j;

  // Process nonlinear constraints.
  for (j = 0 ; j < this_nlc ; j++) {
    (*(this_asl->p.Hvcompd))(this_asl, hv, v, j);

    // Compute dot product g'Hi*v. Should use BLAS.
    for (i = 0, prod = 0 ; i < this_nvar ; i++)
      prod += (hv[i] * g[i]);
    ghjv[j] = prod;
  }
  free(hv);

  // All terms corresponding to linear constraints are zero.
  for (j = this_nlc ; j < this_ncon ; j++) ghjv[j] = 0.;
}

// Return Hessian at (x,y) in triplet form (rows, vals, cols).
void asl_hess(void *asl, double *y, double w, int64_t *rows, int64_t *cols, double *vals) {
  ASL *this_asl = (ASL *)asl;
  double ow[1];  // Objective weight.
  size_t nnzh = (size_t)((*(this_asl->p.Sphset))(this_asl, 0, -1, 1, 1, 1)); // nobj=-1 so ow takes precendence.

  ow[0]  = this_asl->i.objtype_[0] ? -w : w;
  (*(this_asl->p.Sphes))(this_asl, 0, vals, -1, ow, y);

  // Fill in sparsity pattern. Account for 1-based indexing.
  int k = 0;
  for (int i = 0; i < this_asl->i.n_var_; i++)
    for (int j = this_asl->i.sputinfo_->hcolstarts[i]; j < this_asl->i.sputinfo_->hcolstarts[i+1]; j++) {
      rows[k] = this_asl->i.sputinfo_->hrownos[j] + 1;
      cols[k] = i + 1;
      k++;
    }
}
