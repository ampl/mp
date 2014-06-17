/****************************************************************
Copyright (C) 2013 AMPL Optimization LLC; written by David M. Gay.

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

/* variant of qpcheck.c */

#include "nlp.h"

 ssize_t
qpcheckZ_ASL(ASL *asl, fint **rowqp, size_t **colqp, double **delsqp)
{

	int i;
	ssize_t rv;

	ASL_CHECK(asl, ASL_read_fg, "qpcheck");
	i = obj_no;
	if (i < 0 || i >= n_obj)
		return 0;
	rv = nqpcheckZ_ASL(asl, i, rowqp, colqp, delsqp);
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
