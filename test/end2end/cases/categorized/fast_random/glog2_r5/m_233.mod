var x{1..6} integer >= 0 <= 4;

minimize o: x[1] + 2*x[2] - 5*x[3] + 0*x[4] - 3*x[5] - 4*x[6];

s.t. c1: (if x[6] != 0 then x[6] else x[6]) == 1;
s.t. c2: x[5] == 0;
