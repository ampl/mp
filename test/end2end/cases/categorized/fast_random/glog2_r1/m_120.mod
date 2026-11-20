var x{1..3} integer >= 0 <= 4;

minimize o: -4*x[1] - 4*x[2] - x[3];

s.t. c1: (x[3] <= 0) ==> (x[1] < 3);
