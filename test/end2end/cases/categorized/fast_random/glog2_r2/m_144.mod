var x{1..5} integer >= -3 <= 2;

maximize o: 4*x[1] + 6*x[2] - x[3] - 3*x[4] + 3*x[5];

s.t. c1: (x[1] != 2) ==> (x[2] != -3) else (x[4] > -2);
