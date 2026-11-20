var x{1..4} integer >= 0 <= 1;

maximize o: -2*x[1] + 2*x[2] - 4*x[3] - 5*x[4];

s.t. c1: (x[3] < 0) ==> (x[3] > 1) else (x[1] <= 1);
