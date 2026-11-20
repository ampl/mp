var x{1..5} integer >= 0 <= 1;

maximize o: -3*x[1] + 0*x[2] + x[3] + 3*x[4] - 5*x[5];

s.t. c1: (x[4] != 1) and (x[4] <= 1) ==> (x[4] >= 1);
