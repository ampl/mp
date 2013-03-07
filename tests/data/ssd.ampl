function ssd_uniform;
load ../solvers/ssdsolver/ssd.dll;

var x >= 42;
minimize o: x;
subject to ssd: ssd_uniform(x, 0);
