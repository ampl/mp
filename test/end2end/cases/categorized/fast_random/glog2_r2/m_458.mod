var x{1..7} integer >= 0 <= 1;

maximize o: -2*x[1] - 5*x[2] - 4*x[3] + 6*x[4] + 3*x[5] + 5*x[6] + 3*x[7];

s.t. c1: abs(x[2] - x[7]) > 1;
s.t. c2: alldiff(x[4], x[2]);
