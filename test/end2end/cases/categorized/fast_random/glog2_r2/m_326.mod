var x{1..6} integer >= 0 <= 4;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5] + 0*x[6];

s.t. c1: (x[6] > 0) ==> (x[1] > 4) else (x[4] > 4);
