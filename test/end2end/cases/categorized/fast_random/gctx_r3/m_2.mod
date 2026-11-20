var x{1..5} integer >= 0 <= 2;

minimize o: -3*x[1] + 4*x[2] - 5*x[3] + 4*x[4] + 2*x[5];

s.t. c1: (min(x[2], x[2]) != 1) ==> ((min(x[2], x[2]) != 1) or (x[2] <= 2));
s.t. c2: (min(x[2], x[2]) != 1) or (min(x[4], x[5]) <= 0);
