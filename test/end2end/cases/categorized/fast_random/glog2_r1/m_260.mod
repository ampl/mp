var x{1..4} integer >= 0 <= 4;

maximize o: 4*x[1] - 5*x[2] - 3*x[3] + x[4];

s.t. c1: (x[4] > 3) and (x[4] != 2) ==> (x[1] >= 3);
