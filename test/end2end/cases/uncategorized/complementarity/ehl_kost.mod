#  Reference: Michael M. Kostreva,
#  "Elasto-Hydrodynamic Lubrication: a Non-linear Complementarity
#  Problem", Int. Journal for Num. Methods in Fluids (4), 377-397 (1984).

#  The lubricant film gap and the pressure between two lubricated
#  elastic cylinders in line contact are calculated.  When the
#  pressure is positive, Reynolds's equation must be satisfied;
#  when the pressure is 0, the surfaces must diverge.

#  The load (in pounds) is represented by alpha.  The speed
#  (in rpm) of the cylinders is represented by lambda.  A more
#  elaborate version of this model would involve a table of
#  (alpha, lambda) pairs.

#  This variant involves no "defined variables".

param N integer >= 1, default 100; # grid 0..N;
                                   # half-points go 1 .. N
param pi := 4 * atan(1);
param xa := -3;
param xf > xa, := 2;
param dx := (xf - xa) / N;
param alpha := 2.832;
param lambda := 6.057;

param w {i in 0..N} := if i in {0,N} then 0.5 else 1;

var k := 1.6;
var p {i in 1..N} >= 0, := max(0, 1 - abs((xa +1 + i*dx)/2));

s.t. psum: ( 1 - dx*2/pi * sum {i in 1..N} w[i]*p[i] ) == 0;

s.t. reynolds {i in 1..N}: p[i] >= 0 complements
	(lambda / dx)
	  * (
		  ((xa + (i+.5)*dx)^2 + k + 1
		   + 1/pi * (sum {l in 0..N}
				w[l] * (l-i-.5)*dx
				     * log(abs(l-i-.5)*dx)
				     * ((if l < N then p[l+1])
					- (if l > 1 then p[l-1]) )
				)
			)
		- ((xa + (i-.5)*dx)^2 + k + 1
		   + 1/pi * (sum {l in 0..N}
				w[l] * (l-i+.5)*dx
				     * log(abs(l-i+.5)*dx)
				     * ((if l < N then p[l+1])
					- (if l > 1 then p[l-1]) )
				)
			)
		)


      - (1/dx^2) * (
		  ((xa + (i+.5)*dx)^2 + k + 1
		   + 1/pi * (sum {l in 0..N}
				w[l] * (l-i-.5)*dx
				     * log(abs(l-i-.5)*dx)
				     * ((if l < N then p[l+1])
					- (if l > 1 then p[l-1]) )
				)
			)^3
		  * ((if i < N then p[i+1])-p[i]) /
			exp(alpha*((if i < N then p[i+1])+p[i])*.5)
		- ((xa + (i-.5)*dx)^2 + k + 1
		   + 1/pi * (sum {l in 0..N}
				w[l] * (l-i+.5)*dx
				     * log(abs(l-i+.5)*dx)
				     * ((if l < N then p[l+1])
					- (if l > 1 then p[l-1]) )
				)
			)^3 * (p[i] - (if i > 1 then p[i-1])) /
				exp(alpha*(p[i]+(if i > 1 then p[i-1]))*.5)
		)
	>= 0;
