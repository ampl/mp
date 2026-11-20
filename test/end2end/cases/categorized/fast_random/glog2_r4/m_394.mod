var x{1..4} integer >= -3 <= 1;

maximize o: 3*x[1] + 3*x[2] + 0*x[3] - 4*x[4];

s.t. c1: (x[2] != 1) ==> (x[1] != -1);
