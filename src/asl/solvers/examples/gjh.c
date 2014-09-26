/****************************************************************
Copyright (C) 2000 Lucent Technologies
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

/* "solver" gjh computes the gradient g, Jacobian J and Lagrangian  */
/* Hessian H at the current point and writes them to params in file */
/* stub.gjh.  The names "g", "J", "H" can be changed by assignments */
/* in $gjh_options to g, J, and H.  If keyword "sparse" appears in  */
/* $gjh_options, index sets SJ and SH of the structural nonzeros in */
/* J and H are also put in stub.gjh, along with params JS{SJ} and   */
/* HS{SH}, with (dense) J and H derived from JS and HS.		    */

#include "getstub.h"

 static char *gname, *Hname, *Jname;
 static int dense = 1;
 static I_Known dense_kw = { 1, &dense }, sparse_kw = { 0, &dense };

 static void
Sdecl(FILE *f, char *name, int m, int n)
{
	fprintf(f, "set S%s dimen 2; param %sS{S%s};\n",
		name, name, name);
	fprintf(f, "param %s{_i1 in 1..%d, _i2 in 1..%d}%s%s then %s%s",
			name, m, n, "\n\t:= if (_i1,_i2) in S",
			name, name, "S[_i1,_i2] else 0;\n");
	}

 static SufDecl
suftab[] = {
	{ "objweight", 0, ASL_Sufkind_obj |ASL_Sufkind_real, 0 }
	};

 static char *usage_msg[] = {
	"\tto write stub.gjh with params:",
	"\t\tg = gradient",
	"\t\tJ = Jacobian",
	"\t\tH = Hessian",
	"If \"sparse\" appears in $gjh_options, sparse variants JS{SJ} of J",
	"and HS{SH} of H are also defined.",
	"Assignments in $gjh_options can change these names.",
	"For multiple objectives, suffix objweight gives objective weights.",
	0};

 static keyword
keywds[] = {
	KW("dense", IK_val, &dense_kw, "dense-only: no HS{SH} or JS{SJ} (default)"),
	KW("g", C_val, &gname, "g=gradient param name"),
	KW("h", C_val, &Hname, "h=Hessian param name"),
	KW("j", C_val, &Jname, "j=Jacobian param name"),
	KW("sparse", IK_val, &sparse_kw, "supply sparse variants HS of H, JS of J")
	};

 static Option_Info
Oinfo = { "gjh", "gjh", "gjh_options", keywds, nkeywds,
	  ASL_OI_want_funcadd, 0, usage_msg };

 int
main(int argc, char **argv)
{
	ASL *asl;
	FILE *f;
	SufDesc *ow;
	cgrad *cg, **cgx;
	char *s;
	const char *sc;
	fint *hcs, *hrow, i1, i2, nh, nn;
	int i, j, m, n, no;
	real *g, *g1, t, *x, *w, *y;

	asl = ASL_alloc(ASL_read_pfgh);
	s = getstops(argv, &Oinfo);
	if (!s)
		return 1;
	suf_declare(suftab, sizeof(suftab)/sizeof(SufDecl));
	f = jac0dim(s,0);
	x = X0 = M1alloc(n_var*sizeof(real));
	want_xpi0 = 2;
	pfgh_read(f, ASL_findgroups);
	ow = suf_get("objweight", ASL_Sufkind_obj);
	no = 0;
	if (w = ow->u.r)
		no = -1;
	m = n_con;
	if (y = pi0)
		for(i = 0; i < m; i++)
			y[i] = -y[i];
	nh = sphsetup(no, w != 0, y != 0, 0);
	strcpy(stub_end, ".gjh");
	if (!(f = fopen(filename, "w"))) {
		printf("Cannot open \"%s\"\n", filename);
		return 1;
		}
	i = 3*strlen(filename) + 60;
	n = n_var;
	if ((nn = 2*n) < i)
		nn = i;
	if (!gname || !*gname)
		gname = "g";
	if (!Jname || !*Jname)
		Jname = "J";
	if (!Hname || !*Hname)
		Hname = "H";
	fprintf(f, "delete %s, %s, %sS, S%s, %s, %sS, S%s;\n",
		gname, Hname, Hname, Hname, Jname, Jname, Jname);
	fprintf(f, "param %s{1..%d} default 0;\n",
		gname, n);
	if (dense) {
		if (m > 0)
			fprintf(f,
				"param %s{1..%d, 1..%d} default 0;\n",
				Jname, m, n);
		else
			fprintf(f, "delete %s;\n", Jname);
		fprintf(f, "param %s{1..%d, 1..%d} default 0;\n",
			Hname, n, n);
		}
	else {
		if (m > 0)
			Sdecl(f, Jname, m, n);
		Sdecl(f, Hname, n, n);
		}
	if (nn < nh)
		nn = nh;
	if (nn < nzc)
		nn = nzc;
	g = (real*)M1alloc(nn*sizeof(real));
	xknown(x);
	g1 = 0;
	if (w) {
		g1 = g + n;
		memset(g, 0, n*sizeof(real));
		for(i = 0; i < n_obj; i++) {
			if (t = w[i]) {
				objgrd(i,x,g1,0);
				for(j = 0; j < n; j++)
					g[j] += t*g[j];
				}
			}
		}
	else if (n_obj > 0)
		objgrd(0, x, g1 = g, 0);
	fprintf(f, "data;\n");
	if (g1) {
		fprintf(f, "param %s :=\n", gname);
		for(i = 0; i < n; i++)
			if (g[i])
				fprintf(f, "\t%d\t%.g\n", i+1, g[i]);
		fprintf(f, ";\n");
		}
	if (nzc > 0) {
		jacval(x, g, 0);
		if (dense)
			fprintf(f, "param %s :=\n", Jname);
		else
			fprintf(f, "param :S%s: %sS :=\n", Jname, Jname);
		cgx = Cgrad;
		for(i = 1; i <= m; i++)
			if (cg = *cgx++) {
				if (dense) {
					fprintf(f, "[%d,*]\n", i);
					do fprintf(f, "\t%d\t%.g\n", cg->varno+1,
						g[cg->goff]);
					while(cg = cg->next);
					}
				else {
					do fprintf(f, "%d\t%d\t%.g\n", i,
						cg->varno+1, g[cg->goff]);
					while(cg = cg->next);
					}
				}
		fprintf(f, ";\n");
		}

	sc = nh ? "" : " ;";
	if (dense)
		fprintf(f, "param %s :=%s\n", Hname, sc);
	else
		fprintf(f, "param :S%s: %sS :=%s\n", Hname, Hname, sc);
	if (nh) {
		sphes(g, no, w, y);
		hcs = sputinfo->hcolstarts;
		hrow = sputinfo->hrownos;
		i1 = 0;
		for(j = 1; j <= n; j++) {
			i2 = hcs[j];
			if (i1 < i2) {
				if (dense) {
					fprintf(f, "[%d,*]\n", j);
					do fprintf(f, "\t%ld\t%.g\n",
						hrow[i1]+1, g[i1]);
				   	while(++i1 < i2);
					}
				else {
					do fprintf(f, "%d\t%ld\t%.g\n", j,
						hrow[i1]+1, g[i1]);
				   	while(++i1 < i2);
					}
				}
			}
		fprintf(f, ";\n");
		}
	fclose(f);
	sprintf(s = (char *)g,
	 "gjh: \"%s\" written.  Execute\n\n\tinclude \"%s\"\n\tremove \"%s\";\n",
		filename, filename, filename);
	write_sol(s, 0, 0, &Oinfo);
	ASL_free(&asl);
	return 0;
	}
