include gsl.ampl;

param l {1..4};
param u {1..4};

var x {j in 1..4} >= l[j], <= u[j] := 1;

param a := 0.0001;
param b := 1;
param d := 1;
param n := 24;

minimize obj:
  ( a*n - (b*(exp(x[1])-1) - x[3])*x[4]/(exp(x[1]) - 1 + x[4]) )/x[1];

subject to constr1:
    x[3] - 2*gsl_cdf_ugaussian_P(-x[2]) = 0;
subject to constr2:
    x[4] = gsl_cdf_ugaussian_P(-x[2] + d*sqrt(n)) +
           gsl_cdf_ugaussian_P(-x[2] - d*sqrt(n));

data;

param l := 
  1  0.0001
  2  0
  3  0
  4  0
  ;

param u := 
  1  100
  2  100
  3    2
  4    2
  ;

#printf "optimal solution as starting point \n";
#let x[1] := 0.06785874;
#let x[2] := 3.6461717;
#let x[3] := 0.00026617;
#let x[4] := 0.8948622;

display obj;

solve;

display x, obj;

print 'Best known objective value:', -0.920425026;
