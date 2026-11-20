var x{1..4} integer >= -2 <= 0;

maximize o: 2*x[1] - 3*x[2] - 5*x[3] - 3*x[4];

s.t. c1: (x[3] != -2) ==> (x[1] == 0);
s.t. c2: max(x[2], x[4]) - min(x[2], x[1]) == -2;
