var x{1..4} integer >= -2 <= 1;

maximize o: 4*x[1] + 4*x[2] + 4*x[3] - x[4];

s.t. c1: alldiff(x[3], x[1]);
s.t. c2: abs(x[4] - x[2]) != 0;
