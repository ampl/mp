/*******************************************************************
Copyright (C) 2017 AMPL Optimization, Inc.; written by David M. Gay.

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that the copyright notice and this permission notice and warranty
disclaimer appear in supporting documentation.

The author and AMPL Optimization, Inc. disclaim all warranties with
regard to this software, including all implied warranties of
merchantability and fitness.  In no event shall the author be liable
for any special, indirect or consequential damages or any damages
whatsoever resulting from loss of use, data or profits, whether in an
action of contract, negligence or other tortious action, arising out
of or in connection with the use or performance of this software.
*******************************************************************/

#include "nlp.h"

 fint
qpcheck_ASL(ASL *a, fint **rowqp, fint **colqp, double **delsqp)
{

	ASL_fg *asl;
	fint rv;
	int akind, i;

	if ((akind = a->i.ASLtype) != ASL_read_pfgh)
		akind = ASL_read_fg;
	ASL_CHECK(a, akind, "qpcheck");
	asl = (ASL_fg*)a;
	i = obj_no;
	if (i < 0 || i >= n_obj)
		return 0;
	rv = nqpcheck(i, rowqp, colqp, delsqp);
	if (rv < 0) {
		if (rv == -2L)
			fprintf(Stderr,
			 "Quadratic objective involves division by 0.\n");
		else
			fprintf(Stderr,
			 "Sorry, %s can't handle nonlinearities.\n",
				progname ? progname : "");
		exit(1);
		}
	return rv;
	}
