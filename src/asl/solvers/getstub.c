/****************************************************************
Copyright (C) 1997-1998 Lucent Technologies
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

#include "getstub.h"

#ifdef SYMANTEC
extern int _8087;
#endif

static char EmptyString[] = {0};

char *Lic_info_ASL = EmptyString, *Lic_info_add_ASL;

 static int
kw_width(keyword *kw, int n, keyword **pkwe)
{
	const char *s;
	int L = 1, m;
	keyword *kwe = kw;

	if (kw) {
		for(kwe += n; kw < kwe; kw++) {
			m = strlen(kw->name);
			if ((s = kw->desc) && *s == '=') {
				while(*++s > ' ');
				m += s - kw->desc;
				}
			if (L < m)
				L = m;
			}
		}
	*pkwe = kwe;
	return L + 2;
	}

 static void
tabexpand(int L, const char *name, const char *s)
{
	const char *t;
	int k, n;

	for(t = s;; ++t) {
		if (!*t) {
			printf("%-*s%s\n", L, name, s);
			return;
			}
		if (*t == '\n')
			break;
		}
	n = t - s;
	printf("%-*s%.*s\n", L, name, n, s);
	while(*(s = ++t)) {
		while(*t == '\t')
			++t;
		if ((n = t - s)) {
			n *= 8;
			if (n < L)
				n = L;
			for(k = 0; k < n; ++k)
				putchar(' ');
			s = t;
			}
		else if (*s > ' ')
			for(k = 0; k < L; ++k)
				putchar(' ');
		while(*t != '\n') {
			if (!*t++) {
				printf("%s\n", s);
				return;
				}
			}
		n = t - s;
		printf("%.*s\n", n, s);
		}
	}

 static void
shownames(Option_Info *oi)
{
	const char *s;
	keyword *v, *ve;
	int L, L1, L2, anl, te;

	if (oi) {
		te = oi->option_echo & ASL_OI_tabexpand;
		anl = oi->option_echo & ASL_OI_addnewline;
		L = kw_width(oi->keywds, oi->n_keywds, &ve);
		for(v = oi->keywds; v < ve; v++) {
			if ((s = v->desc)) {
				if (*s == '=') {
					while(*++s > ' ');
					L1 = s - v->desc;
					L2 = L - strlen(v->name);
					if (*s)
						s++;
					printf("%s%-*.*s%s\n", v->name, L2,
						L1, v->desc, s);
					}
				else if (te)
					tabexpand(L, v->name, s);
				else
					printf("%-*s%s\n", L, v->name, s);
				}
			else
				printf("%s\n", v->name);
			if (anl)
				putchar('\n');
			}
		}
	exit(0);
	}

 static void
ofix(char **o)
{
	/* get rid of "ix" and possibly "u" */
	char **o1, *s;

	for(o1 = o; (s = *o); o += 2) {
		switch(*s) {
		  case 'i':
		  case 'u':
			continue;
		  }
		o1[0] = o[0];
		o1[1] = o[1];
		o1 += 2;
		}
	*o1 = 0;
	}

 void
usage_noexit_ASL(Option_Info *oi, int rc)
{
	static const char *opts[] = {
		"-", "end of options",
		"=", "show name= possibilities",
		"?", "show usage",
#ifdef SYMANTEC
		"E", "force floating-point emulation (Symantec only)",
#endif
#ifndef NO_BOUNDSFILE_OPTION
		"bf","read boundsfile f",
#endif
		"e", "suppress echoing of assignments",
		"ix","import user-defined functions from x; -i? gives details",
#ifndef NO_BOUNDSFILE_OPTION
		"of", "write .sol file to file f",
#endif
		"s", "write .sol file (without -AMPL)",
		"u", "just show available user-defined functions",
		"v", "just show version",
		0};
	char **o;
	const char *s, *s1;
	keyword *kw, *kwe;
	int i, L, L1, L2;
	FILE *f = stdout;

	if (rc) {
		if (!Stderr)
			Stderr_init_ASL();
		f = Stderr;
		}
	kw = kwe = 0;
	o = 0;
	s = 0;
	L = 2;
	if (oi) {
		s = oi->sname;
		o = oi->usage;
		L = kw_width(kw = oi->options, oi->n_options, &kwe);
		}
	fprintf(f, "usage: %s [options] stub [-AMPL] [<assignment> ...]\n",
		s ? s : basename(progname));
	if (o)
		while((s = *o++))
			fprintf(f, "%s\n", s);
	fprintf(f, "\nOptions:\n");
	o = (char**)opts;
	if (!oi || !(oi->flags & ASL_OI_want_funcadd) || !ix_details_ASL[0])
		ofix(o);
	s = *o;
	for(;;) {
		if (kw < kwe)
			i = s ? strcmp(s, kw->name): 1;
		else if (s)
			i = -1;
		else
			break;
		if (i < 0)
			fprintf(f, "\t-%-*s{%s}\n", L, s, o[1]);
		else {
			if ((s1 = kw->desc))
				if (*s1 == '=') {
					while(*++s1 > ' ');
					L1 = s1 - kw->desc;
					L2 = L - strlen(kw->name);
					if (*s1)
						s1++;
					fprintf(f, "\t-%s%-*.*s{%s}\n",
						kw->name, L2, L1, kw->desc,s1);
					}
				else
					fprintf(f, "\t-%-*s{%s}\n", L,
						kw->name, s1);
			else
				fprintf(f, "\t-%s\n", kw->name);
			++kw;
			}
		if (i <= 0) {
			o += 2;
			s = *o;
			}
		}
	}

 void
usage_ASL(Option_Info *oi, int rc)
{
	usage_noexit_ASL(oi, rc);
	exit(rc);
	}

 char *
get_opt_ASL(Option_Info *oi, char *s)
{
	keyword *kw;
	char buf[256];
	char *b, *be, *s0, *s1;
	fint N;

	while(*s <= ' ' && *s) s++;
	if (!*s)
		return s;
	oi->nnl = 0;
	if (oi->option_echo & ASL_OI_defer_bsname) {
		printf("%s: ", oi->bsname);
		oi->option_echo &= ~ASL_OI_defer_bsname;
		oi->option_echo |= ASL_OI_echo | ASL_OI_echothis;
		}
	s0 = s;
	if ((kw = (keyword *)b_search_ASL(oi->keywds, (int)sizeof(keyword),
			oi->n_keywds, &s, &oi->eqsign))) {
		oi->option_echo = (oi->option_echo | ASL_OI_echothis)
				& ~ASL_OI_badvalue;
		s1 = (*kw->kf)(oi, kw, s);
		if (oi->option_echo & ASL_OI_badvalue) {
			fprintf(Stderr, "Bad value in \"%.*s\"\n",
				(int)(s1-s0), s0);
			oi->n_badopts++;
			while(*++s1 > ' ');
			}
		else if ((oi->option_echo & (ASL_OI_echo|ASL_OI_echothis))
				== (ASL_OI_echo|ASL_OI_echothis))
			printf("%.*s\n", s1-s0, s0);
		}
	else if (*s <= '9' && *s >= '0' && oi->feq) {
		N = strtol(s1 = s, &s1, 10);
		if (*s1 == '=')
			s1++;
		else if (*s1 > ' ')
			goto bad;
		while(*s1 <= ' ')
			if (!*s1++)
				goto bad;
		s = s1;
		while(*++s1 > ' ');
		if (!(oi->option_echo & ASL_OI_never_echo))
			printf("%.*s\n", s1-s0, s0);
		if ((*oi->feq)(&N, s, (fint)(s1-s)))
			oi->n_badopts++;
		}
	else if (oi->kwf) {
		b = buf;
		be = buf + sizeof(buf) - 2;
		while(*s > ' ' && *s != '=') {
			if ((*b = *s++) == '_') {
				if (!(oi->flags & ASL_OI_keep_underscores))
					*b = ' ';
				}
			if (b < be)
				b++;
			}
		*b++ = ' ';
		while(*s <= ' ' && *s)
			s++;
		if (*s != '=' || b >= be)
			goto bad;
		while(*++s <= ' ' && *s);
		while(*s > ' ') {
			*b = *s++;
			if (b < be)
				b++;
			}
		*b = 0;
		if (!(oi->option_echo & ASL_OI_never_echo))
			printf("%.*s\n", s-s0, s0);
		if ((*oi->kwf)(buf, (fint)(b-buf)))
			oi->n_badopts++;
		return s;
		}
	else {
 bad:
		for(s1 = s0; *s1 > ' ' && *s1 != '='; s1++);
		printf("Unknown keyword \"%.*s\"\n", s1-s0, s0);
		if (*s1 == '=')
			while(*++s1 > ' ');
		oi->n_badopts++;
		}
	return s1;
	}

 void
show_version_ASL(Option_Info *oi)
{
	const char *s;
	Const char *ver;
	int L;
	extern Const char *Version_Qualifier_ASL;

	if (!(s = oi->version)
	 && !(s = oi->bsname)
	 && !(s = progname))
		s = "???";	/* defend against bozos */
	L = strlen(s);
	while(L > 0 && s[L-1] == '\n')
		--L;
	if (!(ver = Version_Qualifier_ASL))
		ver = "";
	printf("%s%.*s%s", ver, L, s, oi->nnl ? "\n" : "");
	if (*sysdetails_ASL)
		printf(" (%s)", sysdetails_ASL);
	if (oi->driver_date > 0)
		printf(", driver(%ld)", oi->driver_date);
	printf(", ASL(%ld)\n", ASLdate_ASL);
	if (Lic_info_add_ASL)
		printf("%s\n", Lic_info_add_ASL);
	if (Lic_info_ASL && *Lic_info_ASL)
		printf("%s\n", Lic_info_ASL);
	}

 char *
Ver_val_ASL(Option_Info *oi, keyword *kw, char *v)
{
	char *s;
	int wantver;

	if (!v || *v < '0' || *v > '9')
		goto showver;
	wantver = (int)strtol(s = v, &v, 10);
	if (*(unsigned char*)v > ' ')
		return badval_ASL(oi,kw,s,v);
	if (wantver) {
 showver:
		if (oi->option_echo & ASL_OI_clopt) {	/* -v */
			show_version_ASL(oi);
			exit(0);
			}
		oi->flags |= ASL_OI_show_version;
		}
	else
		oi->flags &= ~ASL_OI_show_version;
	return v;
	}

 static void
ix_usage(VOID)
{
	const char **o = ix_details_ASL, *s;

	printf("-i options:\n");
	while((s = *o++))
		printf("\t%s\n", s);
	exit(0);
	}

 char *
getstub_ASL(ASL *asl, char ***pargv, Option_Info *oi)
{
	char *s, *s1;
	keyword *kw, *okw;
	int i, options;
	char **argv = *pargv;

	progname = *argv;
	if (!Stderr)
		Stderr_init_ASL();	/* set Stderr if necessary */
	if (!asl)
		badasl_ASL(asl,0,"getstub");
	amplflag = 0;
	options = 1;
	okw = 0;
	if (oi) {
		oi->nnl = 0;
		oi->asl = asl;
		okw = oi->options;
		if ((s = getenv("solver_msg"))) {
			i = (int)strtol(s, &s1, 10);
			if (s1 > s && !*s1 && i >= 0 && !(i & 1))
				oi->option_echo = ASL_OI_never_echo;
			}
		oi->option_echo = (oi->option_echo & ASL_OI_showname_bits)
			| ((oi->option_echo & ASL_OI_never_echo)
				? ASL_OI_never_echo : ASL_OI_clopt | ASL_OI_echo);
		oi->n_badopts = 0;
		}

	while((s = *++argv)) {
		if (*s == '-' && options) {
			s1 = s + 1;
			if (okw
			 && (kw = (keyword*)b_search_ASL(okw,
					(int)sizeof(keyword), oi->n_options,
					&s1, &oi->eqsign))) {
				if (!*s1 && argv[1]) {
					s = (*kw->kf)(oi, kw, argv[1]);
					if (s != argv[1])
						argv++;
					continue;
					}
				(*kw->kf)(oi, kw, s1);
				continue;
				}
			if (!s1[1])
			    switch(*s1) {
				case '=':
					shownames(oi);
				case '?':
					usage_ASL(oi,0);
				case '-':
					options = 0;
					continue;
#ifdef SYMANTEC
				case 'E':
					_8087 = 0;
					continue;
#endif
#ifndef NO_BOUNDSFILE_OPTION
				case 'b':
					if ((asl->i.boundsfile = *++argv))
						continue;
					break;
				case 'o':
					if ((asl->i.solfile = *++argv))
						continue;
					break;
#endif
				case 'e':
					if (oi)
						oi->option_echo &= ~ASL_OI_echo;
					continue;
				case 's':
					if (oi)
						oi->wantsol = 1;
					continue;
				case 'u':
					if (!oi || !(oi->flags && ASL_OI_want_funcadd))
						break;
					func_add(asl);
					show_funcs();
					exit(0);
				case 'v':
					if (oi)
						Ver_val_ASL(oi,0,0);
					continue;
				case 'i':
					if (ix_details_ASL[0]) {
						if ((s = argv[1])) {
							argv++;
							if (*s == '?' && !s[1])
								ix_usage();
							i_option_ASL = s;
							}
						continue;
						}
				}
			if (*s1 == 'i' && ix_details_ASL[0]) {
				if (s1[1] == '?' && !s1[2])
					ix_usage();
				i_option_ASL = s1 + 1;
				continue;
				}
#ifndef NO_BOUNDSFILE_OPTION
			if (*s1 == 'b' && s1[1]) {
				asl->i.boundsfile = s1 + 1;
				continue;
				}
			if (*s1 == 'o' && s1[1]) {
				asl->i.solfile =  s1 + 1;
				continue;
				}
#endif
			if (*s1 == '-') {
				if (!strcmp(++s1, "help")) {
					if (oi)
						usage_ASL(oi,0);
					}
				else if (!strcmp(s1, "version")) {
					if (oi) {
						Ver_val_ASL(oi,0,0);
						continue;
						}
					}
				}
			fprintf(Stderr, "%s: bad option %s\n", progname, s);
			usage_ASL(oi,1);
			}
		if (strchr(s,'='))
			break;
		if ((s1 = *++argv) && !strncmp(s1,"-AMPL",5)) {
			amplflag = 1;
			argv++;
			if (s1[5] == 'l') {
				switch(s1[6]) {
				  case 'n':
					oi->option_echo = ASL_OI_never_echo;
					goto no_echo;
				  case 's':
					oi->option_echo = ASL_OI_defer_bsname;
					goto no_echo;
				}}
			if (oi && oi->bsname
			 && !(oi->option_echo & ASL_OI_never_echo))
				need_nl = oi->nnl = printf("%s: ", oi->bsname);
			}
 no_echo:
		i = strlen(s) - 3;
		if (i > 0 && !strcmp(s+i,".nl"))
			s[i] = 0;
		break;
		}
	if (oi && oi->n_badopts)	/* possiby set by kw->kf */
		exit(1);
	*pargv = argv;
	return s;
	}

 int
getopts_ASL(ASL *asl, char **argv, Option_Info *oi)
{
	char *s;

	if (!Stderr)
		Stderr_init_ASL();
	if (!(oi->asl = asl))
		badasl_ASL(asl,0,"getopts");
	if (!oi->option_echo)
		oi->option_echo = ASL_OI_echo;
	oi->option_echo &= ASL_OI_echo | ASL_OI_never_echo | ASL_OI_defer_bsname;
	oi->n_badopts = 0;

	if (oi->opname && (s = getenv(oi->opname)))
		while(*s)
			s = get_opt_ASL(oi, s);

	while((s = *argv++))
		do s = get_opt_ASL(oi, s);
			while(*s);

	need_nl = oi->nnl;
	if (oi->flags & ASL_OI_show_version)
		show_version_ASL(oi);
	fflush(stdout);
	return oi->n_badopts;
	}

 char *
getstops_ASL(ASL *asl, char **argv, Option_Info *oi)
{
	char *s;

	if (!asl)
		badasl_ASL(asl,0,"getstops");
	if (!(s = getstub_ASL(asl, &argv, oi))) {
		fprintf(Stderr, "No stub!\n");
		usage_ASL(oi,1);
		}

	if (getopts_ASL(asl,argv,oi))
		exit(1);
	return s;
	}

 void
badopt_ASL(Option_Info *oi)
{
	oi->n_badopts++;
	oi->option_echo &= ~ASL_OI_echothis;
	}

 char *
badval_ASL(Option_Info *oi, keyword *kw, char *value, char *badc)
{
	char *s;
	int c, w, w1;

	fflush(stdout);
	for(s = badc; *s > ' '; s++);
	w = (int)(strlen(kw->name) + 2 + (badc-value));
	w1 = (int)(s - value);
	fprintf(Stderr, "\n%s%s%.*s\n%*s\nBad character ",
		kw->name, oi->eqsign, w1, value, w, "*");
	c = *(unsigned char *)badc;
	fprintf(Stderr, c >= ' ' && c < 0x7f ? "'%c'" : "'\\x%x'", c);
	fprintf(Stderr, " in numeric string \"%.*s\".\n", w1, value);
	fflush(Stderr);
	badopt_ASL(oi);
	return s;
	}
/* 20000703: basename --> basename_ASL in asl.h */
