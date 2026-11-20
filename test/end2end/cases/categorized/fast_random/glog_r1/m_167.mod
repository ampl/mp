var x{1..5} integer >= -2 <= 2;

minimize o: -4*x[1] - x[2] + 5*x[3] - x[4] + 2*x[5];

s.t. c1: (x[2] < 0) <==> (x[1] + x[3] >= 1);
