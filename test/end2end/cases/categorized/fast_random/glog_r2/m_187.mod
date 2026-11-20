var x{1..5} integer >= 0 <= 3;

maximize o: 6*x[1] + 6*x[2] + 3*x[3] + 2*x[4] - 3*x[5];

s.t. c1: (x[1] <= 1) ==> (x[5] < 1);
