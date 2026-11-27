var x{1..3} integer >= -2 <= 2;

minimize o: 5*x[1] - 4*x[2] - 3*x[3];

s.t. c1: (x[3] == 1) and (x[3] > -2) ==> (x[2] >= 1);
s.t. c2: not ((x[2] > -2) and (x[2] < -2));
s.t. c3: (x[3] == 0) ==> (x[3] <= 0) else (x[3] > 0);
