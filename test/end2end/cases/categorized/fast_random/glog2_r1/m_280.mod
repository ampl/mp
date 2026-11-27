var x{1..6} integer >= 0 <= 2;

maximize o: -2*x[1] - 5*x[2] + 5*x[3] + 4*x[4] - x[5] - x[6];

s.t. c1: (x[5] != 0) ==> (x[1] != 2) else (x[3] == 2);
