var x{1..7} integer >= 0 <= 1;

minimize o: -2*x[1] - 5*x[2] + 0*x[3] - 5*x[4] - 3*x[5] - 5*x[6] + 0*x[7];

s.t. c1: (x[5] != 0) or (x[7] != 0) or (x[4] == 0);
s.t. c2: (if x[4] != 0 then x[3] else x[5]) < 1;
