var x{1..6} integer >= 0 <= 1;

minimize o: -5*x[1] + 0*x[2] - 5*x[3] + 2*x[4] - x[5] + 4*x[6];

s.t. c1: max(x[3], x[5]) - min(x[3], x[1]) <= 1;
s.t. c2: (x[4] != 1) <==> (x[5] + x[6] == 1);
s.t. c3: max(x[2], x[4]) - min(x[3], x[1]) < 0;
