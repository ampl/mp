var x{1..4} integer >= -2 <= 2;

maximize o: 2*x[1] + 3*x[2] + x[3] - x[4];

s.t. c1: alldiff(x[1], x[3]);
s.t. c2: max(x[1], x[1]) - min(x[2], x[4]) <= -3;
