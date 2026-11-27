var x{1..3} integer >= -2 <= 0;

minimize o: x[1] + 2*x[2] - 2*x[3];

s.t. c1: (x[2] <= 1) ==> (x[1] < 0);
