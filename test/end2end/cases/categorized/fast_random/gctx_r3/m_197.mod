var x{1..5} integer >= -1 <= 1;

minimize o: -2*x[1] + 5*x[2] + 0*x[3] - x[4] + 0*x[5];

s.t. c1: (x[5] >= 0) <==> (x[5] >= 0);
s.t. c2: (x[2] <= x[5]) and ((x[5] >= 0) or not (x[5] >= 0));
