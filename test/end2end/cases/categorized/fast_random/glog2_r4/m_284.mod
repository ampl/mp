var x{1..5} integer >= -2 <= 3;

maximize o: x[1] - 4*x[2] + 4*x[3] + 6*x[4] + 4*x[5];

s.t. c1: (x[4] != 3) ==> (x[4] == -1) else (x[3] < 3);
