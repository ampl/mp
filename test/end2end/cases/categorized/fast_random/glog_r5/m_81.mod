var x{1..4} integer >= -2 <= 0;

minimize o: 3*x[1] + x[2] - 4*x[3] + 4*x[4];

s.t. c1: (x[2] < -2) ==> (x[4] == 0);
