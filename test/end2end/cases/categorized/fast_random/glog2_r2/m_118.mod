var x{1..7} integer >= 0 <= 1;

minimize o: 5*x[1] + 3*x[2] - 3*x[3] - 3*x[4] + 0*x[5] + 2*x[6] + 3*x[7];

s.t. c1: (if x[3] <= 0 then x[3] else x[1]) != 0;
