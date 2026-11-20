var x{1..6} integer >= 0 <= 2;

maximize o: -5*x[1] - 2*x[2] + x[3] + 3*x[4] + x[5] - 5*x[6];

s.t. c1: abs(x[4] - x[2]) < 2;
s.t. c2: alldiff(x[1], x[5], x[4], x[6]);
s.t. c3: abs(x[6] - x[3]) != 1;
