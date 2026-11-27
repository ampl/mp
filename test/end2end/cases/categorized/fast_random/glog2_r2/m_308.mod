var x{1..6} integer >= 0 <= 1;

maximize o: 0*x[1] + x[2] - x[3] + 6*x[4] + 0*x[5] - 5*x[6];

s.t. c1: (x[6] >= 0) ==> (x[3] >= 1) else (x[1] > 0);
