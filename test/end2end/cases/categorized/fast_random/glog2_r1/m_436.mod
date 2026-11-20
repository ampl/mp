var x{1..5} integer >= 0 <= 2;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5];

s.t. c1: abs(x[4] - x[4]) < 2;
s.t. c2: alldiff(x[1], x[4], x[2], x[3]);
