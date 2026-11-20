var x{1..4} integer >= -2 <= 0;

minimize o: -x[1] + 0*x[2] - 3*x[3] + 2*x[4];

s.t. c1: (x[4] > -2) <==> (x[1] + x[1] != -1);
