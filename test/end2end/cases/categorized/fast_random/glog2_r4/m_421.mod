var x{1..4} integer >= 0 <= 3;

minimize o: 3*x[1] + 3*x[2] + 3*x[3] + 0*x[4];

s.t. c1: (x[4] <= 4) <==> (x[1] + x[1] > -1);
