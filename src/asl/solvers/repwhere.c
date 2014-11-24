/* The following Lucent copyright notice applies to report_where()
   and the original introuble* routines (moved here in 2011 from
   rops.c and rops2.c). */

/****************************************************************
Copyright (C) 1997 Lucent Technologies
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

/* The following AMPL copyright notice applies to material added in 2011,
   identified by #ifndef ASL_OLD_DERIV_CHECK . */

/****************************************************************
Copyright (C) 2011 AMPL Optimization LLC; written by David M. Gay.

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that the copyright notice and this permission notice and warranty
disclaimer appear in supporting documentation.

The author and AMPL Optimization LLC disclaim all warranties with
regard to this software, including all implied warranties of
merchantability and fitness.  In no event shall the author be liable
for any special, indirect or consequential damages or any damages
whatsoever resulting from loss of use, data or profits, whether in an
action of contract, negligence or other tortious action, arising out
of or in connection with the use or performance of this software.
****************************************************************/

#include "asl.h"
#include "errchk.h"

 void
#ifdef KR_headers
report_where(asl) ASL *asl;
#else
report_where(ASL *asl)
#endif
{
	int i, j, k, k1;
	static const char *what[2] = { "constraint", "objective" };
	static const char *nfmt[2] = { "%d: ", "function: " };
	char *b, buf[512];
	FILE *f;

	fflush(stdout);
	need_nl = 0;
	fprintf(Stderr, "Error evaluating ");

#define next_line fgets(buf,sizeof(buf),f)

	if ((i = cv_index)) {
		strcpy(stub_end, ".fix");
		j = 0;
		if ((f = fopen(filename, "r"))) {
			for(;;) {
				if (!next_line)
					goto eof;
				for(b = buf; *b; b++)
					if (*b == '=') {
						while(++j < i)
							if (!next_line)
								goto eof;
						b = buf;
						while(*b && *b != '=')
							b++;
						if (*b != '=' || b < buf + 2)
							j = 0;
						else
							b[-1] = 0;
						goto eof;
						}
				}
 eof:
			fclose(f);
			}
		if (j == i)
			fprintf(Stderr, "var %s: ", buf);
		else
			fprintf(Stderr, "\"var =\" definition %d: ", i);
		goto ret;
		}

	k = k1 = 0;
	if ((i = co_index) < 0) {
		k = 1;
		i = n_con -i - 1;
		if (n_obj <= 1)
			k1 = 1;
		}
	fprintf(Stderr, "%s ", what[k]);
	if (maxrownamelen) {
		strcpy(stub_end, ".row");
		if ((f = fopen(filename, "r"))) {
			for(j = 0; j <= i; j++)
				if (!next_line)
					break;
			fclose(f);
			if (j > i) {
				for(b = buf; *b; b++)
					if (*b == '\n') {
						*b = 0;
						break;
						}
				fprintf(Stderr, "%s: ", buf);
				goto ret;
				}
			}
		}
	fprintf(Stderr, nfmt[k1], i + 1);
 ret:
	errno = 0;	/* in case it was set by fopen */
	fflush(Stderr);
	}

 static void
jmp_check(Jmp_buf *J, int jv)
{
	if (J)
		longjmp(J->jb, jv);
	}

 static void
Errprint(const char *fmt, ...)
{
	va_list ap;

	va_start(ap, fmt);
#ifndef NO_PERROR
	if (errno)
		fprintf(Stderr, "\n%s: ", strerror(errno));
#endif
	vfprintf(Stderr, fmt, ap);
	va_end(ap);
	fflush(Stderr);
	}

#ifdef ASL_OLD_DERIV_CHECK /*{{*/
 void deriv_errclear_ASL(Edaginfo *I) {}
 void deriv_errchk_ASL(ASL *asl, fint *nerror, int coi, int n) {}
#else /*}{*/
 typedef struct DerrRecord DerrRecord;
 typedef void (*DerrPrint)(ASL*, DerrRecord*);

 struct
DerrRecord {
	DerrPrint errprint;
	const char *fmt, *who;
	real a;
	union { const char *s; real b; } u;
	int jv;
	};

 static void
derrprint1(ASL *asl, DerrRecord *R)
{
	fprintf(Stderr, R->fmt, R->who, R->a);
	}

 static void
derrprint2(ASL *asl, DerrRecord *R)
{
	fprintf(Stderr, R->fmt, R->who, R->a, R->u.b);
	}

 static void
derrprintf(ASL *asl, DerrRecord *R)
{
	fprintf(Stderr, R->fmt, R->who, R->u.s);
	}

 typedef struct DerrMblock DerrMblock;
 struct
DerrMblock {
	DerrMblock *next;
	size_t len;
	real align[1]; /* would prefer align[0], but some older compilers would complain */
	};

 struct
DerivErrInfo {
	DerrMblock *curmb, *freemb;
	char *mbnext, *mblast;
	DerrRecord **R;
	int *busy;
	int nbusy;
	};

 void
deriv_errchk_ASL(ASL *asl, fint *nerror, int coi, int n)
{
	DerivErrInfo *D;
	DerrRecord *R, **Rp, **Rpe;
	int k;

	D = asl->i.Derrs;
	if ((k = coi) < 0) {
		k = -(k + 1);
		if (k >= nlo)
			return;
		}
	else if (k >= nlc)
		return;
	for(Rp = D->R + k, Rpe = Rp + n; Rp < Rpe; ++Rp, ++coi)
		if ((R = *Rp)) {
			jmp_check(err_jmp, R->jv);
			co_index = coi;
			report_where(asl);
			R->errprint(asl,R);
			fflush(Stderr);
			jmp_check(err_jmp1, R->jv);
			exit(1);
			}
	}

 static DerivErrInfo *
new_DerrMblock(Edaginfo *I, size_t len)
{
	DerivErrInfo *D;
	DerrMblock *M, **Mp;
	char *s;
	int nlco;
	size_t L, L1;

	len = len < 4096
		? 4096
		: (len + sizeof(real) - 1) & ~(sizeof(real) - 1);
	if (!(D = I->Derrs)) {
		if ((D = I->Derrs0)) {
			I->Derrs = D;
			M = D->curmb;
			if (M->len >= len)
				return D;
			}
		else {
			nlco = I->nlc_ + I->nlo_;
			L = sizeof(DerivErrInfo)
				+ nlco*(sizeof(int) + sizeof(DerrRecord*));
			L = (L + sizeof(real) - 1) & ~(sizeof(real) - 1);
			L1 = L + (sizeof(DerrMblock) - sizeof(real)) + len;
			D = (DerivErrInfo*)M1alloc_ASL(I, L1);
			memset(D, 0, L);
			I->Derrs = I->Derrs0 = D;
			D->R = (DerrRecord**)(D+1);
			D->busy = (int*)(D->R + nlco);
			M = (DerrMblock*)((char*)D + L);
			M->len = len;
			goto have_M;
			}
		}
	for(Mp = &D->freemb;; Mp = &M->next) {
		if (!(M = *Mp)) {
			M = (DerrMblock*)M1alloc_ASL(I, (sizeof(DerrMblock) - sizeof(real)) + len);
			M->len = len;
			break;
			}
		if (M->len >= len) {
			*Mp = M->next;
			break;
			}
		}
 have_M:
	M->next = D->curmb;
	D->curmb = M;
	D->mbnext = s = (char*)M->align;
	D->mblast = s + M->len;
	return D;
	}

 static DerrRecord *
getDR(ASL *asl)
{
	DerivErrInfo *D;
	DerrRecord *R;
	int k;
	size_t L;

	if ((k = co_index) < 0) {
		k = -(k + 1);
		if (k >= nlo)
			return 0;
		}
	else if (k >= nlc)
		return 0;
	L = (sizeof(DerrRecord) + sizeof(real)-1) & ~(sizeof(real)-1);
	if ((D = asl->i.Derrs)) {
		if (D->R[k])
			return 0;
		if (L <= D->mblast - D->mbnext)
			goto have_D;
		}
	D = new_DerrMblock(&asl->i, L);
 have_D:
	D->R[k] = R = (DerrRecord*)(D->mblast - L);
	D->mblast = (char*)R;
	D->busy[D->nbusy++] = k;
	return R;
	}

 void
deriv_errclear_ASL(Edaginfo *I)
{
	DerivErrInfo *D;
	DerrMblock *M, *M0, *M1;
	DerrRecord **R;
	char *s;
	int *b, *be;

	D = I->Derrs;
	I->Derrs = 0;
	R = D->R;
	for(b = D->busy, be = b + D->nbusy; b < be; ++b)
		R[*b] = 0;
	D->nbusy = 0;
	M0 = D->freemb;
	for(M = D->curmb; M; M0 = M, M = M1) {
		M1 = M->next;
		M->next = M0;
		}
	D->freemb = M0->next;
	M0->next = 0;
	D->curmb = M0;
	D->mbnext = s = (char*)M0->align;
	D->mblast = s + M0->len;
	}

#endif /*}} ASL_OLD_DERIV_CHECK*/

 void
introuble_ASL(ASL *asl, const char *who, real a, int jv)
{
	static const char fmt[] = "can't evaluate %s(%g).\n";
#ifndef ASL_OLD_DERIV_CHECK /*{*/
	DerrRecord *R;

	if (jv > 1 && !(want_deriv & 2)) {
		if ((R = getDR(asl))) {
			R->errprint = derrprint1;
			R->a = a;
			R->jv = jv;
			R->fmt = fmt;
			R->who = who;
			}
		return;
		}
#endif /*}*/
	jmp_check(err_jmp, jv);
	report_where(asl);
	Errprint(fmt, who, a);
	jmp_check(err_jmp1, jv);
	exit(1);
	}

 void
introuble2_ASL(ASL *asl, const char *who, real a, real b, int jv)
{
	static const char fmt[] = "can't evaluate %s(%g,%g).\n";
#ifndef ASL_OLD_DERIV_CHECK /*{*/
	DerrRecord *R;

	if (jv > 1 && !(want_deriv & 2)) {
		if ((R = getDR(asl))) {
			R->errprint = derrprint2;
			R->a = a;
			R->u.b = b;
			R->jv = jv;
			R->fmt = fmt;
			R->who = who;
			}
		return;
		}
#endif /*}*/
	jmp_check(err_jmp, jv);
	report_where(asl);
	Errprint(fmt, who, a, b);
	jmp_check(err_jmp1, jv);
	exit(1);
	}

 void
zero_div_ASL(ASL *asl, real L, const char *op)
{
	errno_set(EDOM);
	jmp_check(err_jmp, 1);
	report_where(asl);
	fprintf(Stderr, "can't compute %g%s0.\n", L, op);
	fflush(Stderr);
	jmp_check(err_jmp1, 1);
	exit(1);
	}

 void
fintrouble_ASL(ASL *asl, func_info *fi, const char *s, TMInfo *T)
{
	TMInfo *T1, *T1prev;
	int jv;
	static const char fmt[] = "Error in function %s:\n\t%s\n";

	jv = 1;
	switch(*s) {
	 case '\'':
		jv = 2;
		goto inc_s;
	 case '"':
		jv = 3;
 inc_s:
		++s;
	 }
#ifndef ASL_OLD_DERIV_CHECK /*{*/
	if (jv > 1 && !(want_deriv & 2)) {
		DerivErrInfo *D;
		DerrRecord *R;
		size_t L;

		if ((R = getDR(asl))) {
			D = asl->i.Derrs;
			L = strlen(s) + 1;
			if (L > D->mblast - D->mbnext)
				D = new_DerrMblock(&asl->i, L);
			memcpy(D->mbnext, s, L);
			R->u.s = D->mbnext;
			D->mbnext += L;
			R->errprint = derrprintf;
			R->jv = jv;
			R->fmt = fmt;
			R->who = fi->name;
			}
		return;
		}
#endif /*}*/
	jmp_check(err_jmp, jv);
	report_where(asl);
	fprintf(Stderr, fmt, fi->name, s);
	fflush(Stderr);
	for(T1 = T->u.prev; T1; T1 = T1prev) {
		T1prev = T1->u.prev;
		free(T1);
		}
	jmp_check(err_jmp1,jv);
	exit(1);
	}
