var x{1..5} integer >= 0 <= 3;

minimize o: -x[1] + 2*x[2] - 2*x[3] - 5*x[4] + 5*x[5];

s.t. c1: not ((x[4] != 0) and (x[3] == 2));
