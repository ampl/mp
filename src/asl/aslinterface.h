// ===============================================================
// Generic Ampl interface to facilitate usage from other languages
// Dominique Orban
// Vancouver, April 2014.
// Montreal, February 2015.
// ===============================================================

#include <stdint.h>
#include <sys/types.h>                 // For ssize_t
#include "asl_pfgh.h"

// ==========================================================================

//
//        P r o t o t y p e s   f o r   m o d u l e   f u n c t i o n s

// ==========================================================================

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ASL ASL;

ASL *asl_init(char *stub);
void asl_finalize(ASL *asl);
void asl_write_sol(ASL *asl, const char *msg, double *x, double *y);

int asl_objtype(ASL *asl);
int asl_nvar(   ASL *asl);
int asl_ncon(   ASL *asl);
int asl_nlc(    ASL *asl);
int asl_nlnc(   ASL *asl);
int asl_nnzj(   ASL *asl);
int asl_nnzh(   ASL *asl);
int asl_islp(   ASL *asl);

double *asl_x0(  ASL *asl);
double *asl_y0(  ASL *asl);
double *asl_lvar(ASL *asl);
double *asl_uvar(ASL *asl);
double *asl_lcon(ASL *asl);
double *asl_ucon(ASL *asl);

void asl_varscale(ASL *asl, double *s);
void asl_lagscale(ASL *asl, double  s);
void asl_conscale(ASL *asl, double *s);

double  asl_obj(     ASL *asl, double *x);
void    asl_grad(    ASL *asl, double *x, double *g);
void    asl_cons(    ASL *asl, double *x, double *c);
double  asl_jcon(    ASL *asl, double *x, int j);
void    asl_jcongrad(ASL *asl, double *x, double *g, int j);
void    asl_hprod(   ASL *asl, double *y, double *v, double *hv, double w);
void    asl_hvcompd( ASL *asl, double *v, double *hv, int nobj);
void    asl_ghjvprod(ASL *asl, double *g, double *v, double *ghjv);

size_t asl_sparse_congrad_nnz(ASL *asl, int j);
void asl_sparse_congrad(ASL *asl, double *x, int j, int64_t *inds, double *vals);
void asl_jac( ASL *asl, double *x, int64_t *rows, int64_t *cols, double *vals);
void asl_hess(ASL *asl, double *y, double w, int64_t *rows, int64_t *cols, double *vals);

#ifdef __cplusplus
}  // extern "C"
#endif
