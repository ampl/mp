var x{1..4} integer >= 0 <= 4;

maximize o: -333*x[1] + 17*x[2] + 1000*x[3] + 3000*x[4];

s.t. c1: abs(x[4] - x[3]) >= 4;
s.t. c2: alldiff(x[3], x[2], x[4]);
