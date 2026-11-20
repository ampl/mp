var x{1..6} integer >= -2 <= 1;

maximize o: 3*x[1] - 4*x[2] + x[3] - 5*x[4] - 2*x[5] + 5*x[6];

s.t. c1: (x[1] < -2) ==> (x[6] >= -2) else (x[4] < -2);
