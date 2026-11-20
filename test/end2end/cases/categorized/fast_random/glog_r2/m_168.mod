var x{1..3} integer >= -2 <= 2;

maximize o: 5*x[1] + 5*x[2] + 4*x[3];

s.t. c1: max(x[3], x[1]) - min(x[2], x[3]) > 1;
s.t. c2: abs(x[3] - x[1]) <= 2;
s.t. c3: (x[3] <= -1) and (x[3] > -1) ==> (x[3] >= -1);
