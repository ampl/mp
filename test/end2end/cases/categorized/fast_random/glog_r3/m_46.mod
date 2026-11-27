var x{1..6} integer >= 0 <= 4;

maximize o: -4*x[1] - 5*x[2] - 4*x[3] - 2*x[4] + x[5] - 5*x[6];

s.t. c1: (x[6] >= 2) ==> (x[5] == 2) else (x[3] != 3);
