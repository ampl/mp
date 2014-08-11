#include "asl.h"

 static real *
ones(ASL *asl, int n)
{
	real *x, *x0, *xe;

	x = x0 = (real*)mem_ASL(asl, n*sizeof(real));
	xe = x + n;
	while(x < xe)
		*x++ = 1.;
	return x0;
	}

 static int
zcheck(ASL *asl, int i, real s, int n, fint *ierror, char *who)
{
#undef W0
#ifdef IEEE_MC68k
#define W0 0
#endif
#ifdef IEEE_8087
#define W0 1
#endif
#ifdef W0
	union { real r; unsigned int ui[2]; } u;
#define Inftest(z) (u.r = z, (u.ui[W0] & 0x7ff00000) == 0x7ff00000)
#else
#define Inftest(z) (z <= negInfinity || z >= Infinity)
#endif
	int rt;
	static const char *RN[6] = { "???", "f", "fg", "fgh", "pfg", "pfgh" };

	if ((n >= 0 && (i < 0 || i >= n))
	 || s == 0.
	 || Inftest(s)) {
			if (ierror && *ierror >= 0) {
				*ierror = 1;
				return 1;
				}
			fprintf(Stderr, "%s(", who);
			if (n >= 0)
				fprintf(Stderr, "%d, ", i);
			fprintf(Stderr, "%.g, nerror): bad argument\n", s);
 byebye:
			fflush(Stderr);
			if (err_jmp1)
				longjmp(err_jmp1->jb, 1);
			exit(1);
			}
	if (!Ograd) {
		if (ierror && *ierror >= 0) {
			*ierror = 1;
			return 1;
			}
		if ((rt = asl->i.ASLtype) < 1 || rt > 5)
			rt = 0;
		fprintf(Stderr, "%s called before %s_read().\n", who, RN[rt]);
		goto byebye;
		}
	if (ierror && *ierror >= 0)
		*ierror = 0;
	return 0;
	}

 static void
scaleadj(real s, int i, int m, real *scale, real *L, real *U, real *x)
{
	real u, v;

	scale += i;
	if (x)
		x[i] /= s;
	if (!U) {
		U = L + 1;
		i <<= 1;
		}
	L += i;
	U += i;
	*scale *= s;
	if (s > 0.) {
		if (*L > negInfinity) {
			if (m)
				*L *= s;
			else
				*L /= s;
			}
		if (*U < Infinity) {
			if (m)
				*U *= s;
			else
				*U /= s;
			}
		}
	else {
		u = -*L;
		v = -*U;
		if (u < Infinity) {
			if (m)
				u = *L * s;
			else
				u = *L / s;
			}
		if (v > negInfinity) {
			if (m)
				v = *U * s;
			else
				v = *U / s;
			}
		*L = v;
		*U = u;
		}
	}

 void
conscale_ASL(ASL *asl, int i, real s, fint *ierror)
{
	static char who[] = "conscale";

	if (!asl
	 || asl->i.ASLtype < ASL_read_fg
	 || asl->i.ASLtype > ASL_read_pfgh)
		badasl_ASL(asl, ASL_read_fg, who);
	if (zcheck(asl, i, s, n_con, ierror, who))
		return;
	if (s == 1.)
		return;
	if (!asl->i.cscale)
		asl->i.cscale = ones(asl, n_con);
	if (!asl->i.lscale)
		asl->i.lscale = asl->i.cscale;
	scaleadj(s, i, 1, asl->i.cscale, LUrhs, Urhsx, pi0);
	if (asl->i.lscale != asl->i.cscale)
		asl->i.lscale[i] *= s;
	}

 void
varscale_ASL(ASL *asl, int i, real s, fint *ierror)
{
	static char who[] = "varscale";

	if (!asl
	 || asl->i.ASLtype < ASL_read_fg
	 || asl->i.ASLtype > ASL_read_pfgh)
		badasl_ASL(asl, ASL_read_fg, who);
	if (zcheck(asl, i, s, n_var, ierror, who))
		return;
	if (!asl->i.vscale)
		asl->i.vscale = ones(asl, n_var);
	scaleadj(s, i, 0, asl->i.vscale, LUv, Uvx, X0);
	}

 void
lagscale_ASL(ASL *asl, real s, fint *ierror)
{
	static char who[] = "lagscale";
	real *c, *l, *le;
	size_t L;

	if (!asl
	 || (asl->i.ASLtype != ASL_read_pfgh
	  && asl->i.ASLtype != ASL_read_fgh))
		badasl_ASL(asl, ASL_read_pfgh, who);
	if (zcheck(asl, 0, s, -1, ierror, who))
		return;
	if (s == 1.)
		return;
	if (!asl->i.lscale)
		asl->i.lscale = ones(asl, n_con);
	else if (asl->i.lscale == asl->i.cscale) {
		L = n_con*sizeof(real);
		memcpy(asl->i.lscale = (real*)mem_ASL(asl,L), asl->i.cscale, L);
		}
	l = asl->i.lscale;
	le = l + n_con;
	if ((c = asl->i.cscale)) {
		while(l < le)
			*l++ = s * *c++;
		}
	else
		while(l < le)
			*l++ *= s;
	if ((c = pi0)) {
		le = c + n_con;
		s = 1. / s;
		while(c < le)
			*c++ *= s;
		}
	}
