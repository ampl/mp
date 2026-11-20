var x{1..3} integer >= 0 <= 2;

minimize o: -4*x[1] + 3*x[2] - 2*x[3];

s.t. c1: (x[2] < 2) and (x[1] >= 0) ==> (x[1] >= 2);
