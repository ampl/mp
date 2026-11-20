var x{1..4} integer >= 0 <= 4;

maximize o: -5*x[1] - 5*x[2] + 2*x[3] + x[4];

s.t. c1: max(x[1], x[3]) - min(x[3], x[1]) == 1;
s.t. c2: (x[4] < 0) ==> (x[1] == 1) else (x[1] <= 4);
