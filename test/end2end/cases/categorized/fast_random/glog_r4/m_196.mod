var x{1..3} integer >= -2 <= 0;

minimize o: -3*x[1] - 5*x[2] + 2*x[3];

s.t. c1: (x[3] == 0) or (x[2] > 0) or (x[3] == -2);
s.t. c2: (x[3] != -2) ==> (x[3] < 0);
