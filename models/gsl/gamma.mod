# Find minimum of the gamma function for x > 0.

include gsl.ampl;

var x >= 1e-5;
minimize obj: gsl_sf_gamma(x);
solve;
print x;
