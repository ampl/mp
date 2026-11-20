var x{1..4} integer >= -2 <= 1;

maximize o: 2*x[1] - 2*x[2] + x[3] + 3*x[4];

s.t. c1: (x[4] != -1) ==> (x[2] > -1);
