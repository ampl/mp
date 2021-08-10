# choi.mod

param M integer >= 1;	# subjects or consumers
param N integer >= 1;	# brand names

set subjects := 1 .. M;
set brands := 1 .. N;
set ingred;

param chi >= 0;		# randomness increases as chi decreases
param K;		# constant term in probability function
			# representing "no purchase" option
param x {brands, ingred};	# amount of ingredients
param y {subjects, ingred};	# preferences
param v {subjects};		# importance weights
param w0 {subjects};		# unscaled w values
param w{i in subjects} := -chi * w0[i];	# constant in DU eqn
param b {subjects};		# constant in DU eqn
param c {brands};			# average cost for producing j
param DU {i in subjects, j in brands} :=
	-chi * (v[i] * sum {k in ingred} (x[j,k]-y[i,k])^2 + b[i]);

param p_lo {j in brands} default c[j];
param p_up {j in brands} default Infinity;

var p {j in brands} := c[j] + .01;

s.t. mprofit {j in brands}:	/* marginal profit for brand j */
	p_lo[j] <= p[j] <= p_up[j]
     complements
	(-1/M) * sum {i in subjects}
		( exp(w[i]*p[j]+DU[i,j])
	 	/ ( K + sum {jj in brands} exp(w[i]*p[jj]+DU[i,jj]) )
		* (1 + (p[j]-c[j]) * w[i] *
		   ( K + sum {jj in brands:jj != j} exp(w[i]*p[jj]+DU[i,jj]) )
		   / ( K + sum {jj in brands} exp(w[i]*p[jj]+DU[i,jj]) )
		  ) );
