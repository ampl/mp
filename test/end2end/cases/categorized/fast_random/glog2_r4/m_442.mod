var x{1..6} integer >= 0 <= 4;

maximize o: 3*x[1] + 6*x[2] - 5*x[3] - x[4] - 5*x[5] - 4*x[6];

s.t. c1: (x[1] >= 3) ==> (x[3] > 2) else (x[6] < 4);
