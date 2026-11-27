var x{1..4} integer >= 0 <= 1;

maximize o: 2*x[1] - 3*x[2] + 2*x[3] + 3*x[4];

s.t. c1: (x[3] > 1) ==> (x[4] > 0) else (x[1] == 1);
