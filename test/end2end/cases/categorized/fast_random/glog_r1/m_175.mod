var x{1..5} integer >= 0 <= 3;

maximize o: -3*x[1] + x[2] - 4*x[3] - 4*x[4] - x[5];

s.t. c1: (x[5] > 3) ==> (x[5] <= 3) else (x[2] == 3);
