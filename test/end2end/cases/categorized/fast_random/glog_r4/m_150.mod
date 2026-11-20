var x{1..5} integer >= -2 <= 2;

minimize o: x[1] + 4*x[2] + 2*x[3] - 4*x[4] + x[5];

s.t. c1: (x[1] != -2) and (x[3] <= -2) ==> (x[1] >= -2);
