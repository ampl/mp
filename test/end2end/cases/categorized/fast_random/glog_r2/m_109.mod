var x{1..6} integer >= 0 <= 4;

maximize o: 0*x[1] + 2*x[2] + 4*x[3] + x[4] + 5*x[5] - 3*x[6];

s.t. c1: (x[6] <= 0) ==> (x[6] < 3);
