var x{1..3} integer >= 0 <= 2;

maximize o: 4*x[1] + 4*x[2] + 4*x[3];

s.t. c1: (x[2] <= 0) ==> (x[2] >= 0) else (x[2] <= 0);
s.t. c2: max(x[3], x[3]) - min(x[3], x[3]) >= 1;
