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

 struct
Objrep {
	/* For minimization (or maximization, with suitable inequality */
	/* sense adjustments) of Objective = c0 + c1*v */
	/* with v = _var[ivo] defined by constraint ico of the form */
	/* (_con[ico].body = c2*v + f(x)) >= rhs  with c2 > 0   or */
	/* (_con[ico].body = c2*v + f(x)) <= rhs  with c2 < 0, */
	/* We change the Objective to */
	/* Objective = c0 + c1*(rhs - f(x))/c2 = c0a + c12*f(x) */
	/* with c12 = -c1 / c2  and  c0a = c0 - c12*rhs. */
	/* Objective gradient = c12 * grad(f(x)). */
	/* We remove constraint _con[ico] and variable v; writesol */
	/* calls obj_adj_xy_ASL to compute v = (Objectve - c0)/c1 */
	/* and dual variable value -c12 for the removed constraint. */

	int ico;	/* index of constraint that gives the objective */
	int ivo;	/* index of objective variable */
	int nxval;	/* X generation for f */
	void (*opify)(ASL*);	/* whether obj_adj_xy may need to call qp_opify */
	real c0, c0a, c1, c12, f;
	cgrad *cg;	/* Original gradient if modified by mqpcheck(). */
	cgrad *cg0;	/* Copy of cg for use in obj_adj_xy_ASL.  For solvers */
			/* like minos and snopt that separately evaluate the */
			/* linear and nonlinear parts of objectives, cg and cg0 */
			/* may differ. */
	};
