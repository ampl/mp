var x{1..3} integer >= 0 <= 2;

minimize o: -2*x[1] + 5*x[2] - 5*x[3];

s.t. c1: (x[3] > 1) and (x[2] < 2) ==> (x[3] >= 1);
