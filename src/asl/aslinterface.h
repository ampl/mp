// ===============================================================
// Generic Ampl interface to facilitate usage from other languages
// Dominique Orban
// Vancouver, April 2014.
// Montreal, February 2015.
// ===============================================================

#include <sys/types.h>                 // For ssize_t
#include "asl_pfgh.h"

// ==========================================================================

//
//        P r o t o t y p e s   f o r   m o d u l e   f u n c t i o n s

// ==========================================================================

void *asl_init(char *stub);
void asl_finalize(void *asl);
void asl_write_sol(void *asl, const char *msg, double *x, double *y);

int asl_objtype(void *asl);
int asl_nvar(   void *asl);
int asl_ncon(   void *asl);
int asl_nlc(    void *asl);
int asl_nlnc(   void *asl);
int asl_nnzj(   void *asl);
int asl_nnzh(   void *asl);
int asl_islp(   void *asl);

double *asl_x0(  void *asl);
double *asl_y0(  void *asl);
double *asl_lvar(void *asl);
double *asl_uvar(void *asl);
double *asl_lcon(void *asl);
double *asl_ucon(void *asl);

void asl_varscale(void *asl, double *s);
void asl_lagscale(void *asl, double  s);
void asl_conscale(void *asl, double *s);

double  asl_obj(     void *asl, double *x);
void    asl_grad(    void *asl, double *x, double *g);
void    asl_cons(    void *asl, double *x, double *c);
double  asl_jcon(    void *asl, double *x, int j);
void    asl_jcongrad(void *asl, double *x, double *g, int j);
void    asl_hprod(   void *asl, double *y, double *v, double *hv, double w);
void    asl_hvcompd( void *asl, double *v, double *hv, int nobj);
void    asl_ghjvprod(void *asl, double *g, double *v, double *ghjv);

size_t asl_sparse_congrad_nnz(void *asl, int j);
void asl_sparse_congrad(void *asl, double *x, int j, int64_t *inds, double *vals);
void asl_jac( void *asl, double *x, int64_t *rows, int64_t *cols, double *vals);
void asl_hess(void *asl, double *y, double w, int64_t *rows, int64_t *cols, double *vals);
