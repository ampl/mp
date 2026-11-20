var x{1..5} integer >= 0 <= 3;

maximize o: -2*x[1] + 3*x[2] + 2*x[3] + x[4] - 2*x[5];

s.t. c1: (x[4] == 1) ==> (x[5] != 3) else (x[2] > 3);
s.t. c2: max(x[3], x[4]) - min(x[3], x[1]) > 3;
s.t. c3: (x[4] == 1) ==> (x[5] != 3) else (x[2] > 3);
s.t. c4: max(x[3], x[4]) - min(x[3], x[1]) > 3;
