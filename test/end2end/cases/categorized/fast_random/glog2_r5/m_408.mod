var x{1..5} integer >= 0 <= 1;

maximize o: 2*x[1] + 0*x[2] - 4*x[3] + 4*x[4] + x[5];

s.t. c1: alldiff(x[3], x[5], x[1]);
s.t. c2: abs(x[3] - x[4]) < 1;
