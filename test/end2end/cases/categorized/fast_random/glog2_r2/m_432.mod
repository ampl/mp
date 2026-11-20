var x{1..6} integer >= 0 <= 4;

minimize o: -2*x[1] - x[2] + 5*x[3] - 5*x[4] + 2*x[5] + 0*x[6];

s.t. c1: (if x[1] > 3 then x[1] else x[1]) >= 0;
s.t. c2: (x[5] > 3) <==> (x[1] + x[5] == 1);
