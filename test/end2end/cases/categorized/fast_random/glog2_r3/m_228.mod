var x{1..7} integer >= 0 <= 1;

minimize o: 2*x[1] + 5*x[2] - 2*x[3] - 4*x[4] - x[5] - 3*x[6] - 4*x[7];

s.t. c1: (if x[6] < 1 then x[3] else x[3]) > 0;
