var x{1..6} integer >= -2 <= 2;

minimize o: -2*x[1] + x[2] - 2*x[3] + 2*x[4] - 2*x[5] + 5*x[6];

s.t. c1: abs(x[1] - x[4]) >= 1;
s.t. c2: alldiff(x[6], x[2]);
s.t. c3: abs(x[1] - x[4]) >= 1;
s.t. c4: alldiff(x[6], x[2]);
