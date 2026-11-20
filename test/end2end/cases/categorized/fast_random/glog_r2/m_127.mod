var x{1..4} integer >= -2 <= 2;

maximize o: x[1] + 4*x[2] + 3*x[3] - 3*x[4];

s.t. c1: alldiff(x[4], x[2], x[1]);
s.t. c2: max(x[3], x[4]) - min(x[4], x[4]) >= 2;
s.t. c3: (x[2] == 2) ==> (x[3] > 0);
