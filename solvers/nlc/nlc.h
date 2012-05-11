/****************************************************************
Copyright (C) 1992, 1993, 1994 Lucent Technologies
All Rights Reserved

Permission to use, copy, modify, and distribute this software and
its documentation for any purpose and without fee is hereby
granted, provided that the above copyright notice appear in all
copies and that both that the copyright notice and this
permission notice and warranty disclaimer appear in supporting
documentation, and that the name of Lucent or any of its entities
not be used in advertising or publicity pertaining to
distribution of the software without specific, written prior
permission.

LUCENT DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
IN NO EVENT SHALL LUCENT OR ANY OF ITS ENTITIES BE LIABLE FOR ANY
SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.
****************************************************************/

#ifdef KR_headers
#define ANSI(x) ()
#else
#define ANSI(x) x
#endif

 typedef struct
v_i {
	union {
		real v;
		struct v_i *next;
		} u;
	int i;
	} v_i;

 typedef union
vpi {
	int	i;
	real	*vp;
	cgrad	*cg;
	ograd	*og;
	v_i	*vi;
	} vpi;

 typedef struct
Adjoint {
			vpi o;	/* offset into pd or dv, if appropriate */
			unsigned storage:4;
#define STOR_UNUSED	0
#define STOR_IMPLICIT	1
#define STOR_PD		2
#define STOR_DV		3
#define STOR_VP		4
#define STOR_VARVAL	5
#define STOR_GRAD	6
#define STOR_JAC	7
#define STOR_VI		8
#define STOR_DEFV	9
			unsigned neg:1;
			unsigned ifset:1;
			unsigned stored:1;
			unsigned seen:1;
		} Adjoint;

 typedef struct
dLR {
		int kind;
		union {
			double *vp;
			expr_if *eif;
			expr_va *eva;
			expr *ep;
			int i;
			} o;
		} dLR;

#define dLR_UNUSED	0
#define dLR_PD		1
#define dLR_one		2
#define dLR_negone	3
#define dLR_VP		4
#define dLR_VARVAL	5
#define dLR_VARARG	6
#define dLR_IF		7

#ifdef X64_bit_pointers
#define Adjp(x) (*(Adjoint **)x)
#define dLRp(x) (*(dLR **)&x)
#define Make_dLR(x) *(dLR **)x = (dLR *)mem(sizeof(dLR))
#define Make_dLRp(x) (*(dLR ***)&x = (dLR **)mem(sizeof(dLR*)), (dLR **)x)
#else
#define Adjp(x) ((Adjoint *)x)
#define dLRp(x) ((dLR *)&x)
#define Make_dLR(x) (dLR *)x
#define Make_dLRp(x) (dLR **)&x
#endif

 typedef char *efuncb(expr *, char *);
#define callb(a,b) (*(efuncb *)a->op)(a, b)

 typedef struct
expr_nx {	/* for numbers */
	efuncb *op;
	struct expr_nx *next;
	real v;
	} expr_nx;

extern char *e_val ANSI((expr*, char*));
extern char *f_OPNUM1 ANSI((expr*, char*));
extern void ifstart ANSI((char*, char*, char*));
extern void assign ANSI((char*, char*));
extern void binop ANSI((char*, char*, char*, char*));
extern void Goto ANSI((int));
extern void ifgo ANSI((char*, char*, char*, int));
extern void label ANSI((int));
extern void ifstart ANSI((char*, char*, char*));
extern void elsestart ANSI((void));
extern void elseif ANSI((char*, char*, char*));
extern void endif ANSI((void));
extern int  Switch ANSI((char*, int));
extern void Case ANSI((int));
extern void Break ANSI((int));
extern void endswitch ANSI((int));
extern void domain ANSI((char*, char*));
extern void zerdiv ANSI((char*));
extern void call ANSI((char*, char*));
extern char *call1 ANSI((char*, char*));
extern char *call2 ANSI((char*, char*, char*));
extern void introuble ANSI((char*, char*));
extern void introuble2 ANSI((char*, char*, char*));
extern char *num ANSI((int));

extern char *Half, *Negone, *One, *T, *T1, *Zero;
extern char *opEQ, *opGE, *opGT, *opLE, *opLT, *opNE, *pd_fmt;
extern int branches;
extern char *cond_fmt, *offlfmt1, *offlfmt2, *progname;
