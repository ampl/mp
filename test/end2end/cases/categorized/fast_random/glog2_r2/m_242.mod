var x{1..5} integer >= -2 <= 1;

maximize o: 3*x[1] - 2*x[2] + 6*x[3] - 2*x[4] + 3*x[5];

s.t. c1: (x[3] >= -2) ==> (x[2] <= 0) else (x[1] < 0);
