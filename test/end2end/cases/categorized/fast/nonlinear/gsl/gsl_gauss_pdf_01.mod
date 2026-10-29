#####################################
#
#####################################

function gsl_ran_gaussian_pdf;


var x >=-11 <=5;

var obj_var = gsl_ran_gaussian_pdf(x+2, 1.8);

minimize Obj: obj_var;

s.t. F2F1: obj_var >= gsl_ran_gaussian_pdf(x+6, 1.6);

s.t. F2F3: obj_var >= gsl_ran_gaussian_pdf(x-1.36, 2.03);