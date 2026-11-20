var x{1..3} integer >= -2 <= 0;

maximize o: -5*x[1] + 2*x[2] + 4*x[3];

s.t. c1: (x[2] == -1) or (x[1] > 0) or (x[3] < -2);
s.t. c2: (x[1] <= -1) ==> (x[3] >= -2) else (x[1] > -1);
